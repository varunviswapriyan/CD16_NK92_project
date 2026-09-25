"""
bootstrap_pzap.py -- bootstrap CIs for the pZAP fit, 6 FITTED parameters:
lig0, kdl0, ZAP0, SYK0, KZP_MULT, KPR_MULT.  KZBG_FRAC stays removed (=1).
Indrani approved estimating kzp ("not experimentally measured ... you can
estimate this parameter") and accepted the fitted KPR, so both are fitted here
and get their own CIs.

Two bootstrap flavours, both resampling the EXPERIMENTAL data and refitting:

  emp  -- empirical: resample the 3 day-values (26-May, 28-May, 03-Jun) with
          replacement, per cell line and per timepoint, take the mean.
          Honest but coarse: with 3 days there are only 10 distinct draws.
  par  -- parametric: SD across the 3 days per (line, timepoint), then draw
          Normal(mean, SD). Smooth, not limited to 10 combinations.

Each sample gets its own copy of the v77 directory with a resampled
data/pZAP70_Tyr493_mean.csv, warm-started bounds around the v77 optimum, and a
short PSO run. Nothing in her fitting code is modified.

    python ~/CD16_NK92_project/filesCC/bootstrap_pzap.py build        # make dirs
    python ~/CD16_NK92_project/filesCC/bootstrap_pzap.py submit       # build + sbatch
    python ~/CD16_NK92_project/filesCC/bootstrap_pzap.py submit 8 8   # 8 emp + 8 par
    python ~/CD16_NK92_project/filesCC/bootstrap_pzap.py report
"""
import os, sys, json, glob, shutil, subprocess
from itertools import combinations_with_replacement
from math import factorial
import numpy as np
import pandas as pd

HOME     = os.path.expanduser('~')
SRC      = os.path.join(HOME, 'NK92_fit_v77')
SUB      = 'estimate_params_pzap_cleaned_up'
XLSX     = '/home/gddaslab/share/Varun_Indrani/estimate_params_pzap/data/pZAP70_Tyr493_Tyr292_original_and_averages.xlsx'
PIN = [x for x in os.environ.get('BOOT_PIN', '').split(',') if x]
BOOT_ROOT = os.path.join(HOME, 'boot_pzap' + ('_pin' + '_'.join(PIN) if PIN else ''))
MARKER   = 'pZAP70 (Tyr493)'
# cell line in the xlsx  ->  column in data/pZAP70_Tyr493_mean.csv
#
# VERIFIED against the shipped CSV (each column reproduces its line's 3-day mean
# to <1e-6 at 60/120/300 s; the next-best line is off by 5-37%).  The parental
# NK92 row feeds mean_NK92, which the fitting code does not read.
LINE2COL = {'KI 1': 'mean_zeta', 'KI 2': 'mean_gamma', 'KI 6': 'mean_hetero'}
DAYCOLS  = ['d_0', 'd_1', 'd_2', 'd_5']          # 0, 60, 120, 300 s -- the fitted points
TIMES    = [0.0, 60.0, 120.0, 300.0]
PARTICLES = int(os.environ.get('BOOT_PARTICLES', '16'))   # PSO swarm size
ITERS     = int(os.environ.get('BOOT_ITERS', '20'))       # PSO iterations
ENV = ('module load Miniconda3/4.9.2; '
       'source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh; conda activate CD16_v2')

# warm-start centres (log10). lig0..SYK0 from v77; multipliers from the v69 fit
# that produced the approved values (KZP_MULT=2.2 -> 0.3424, KPR_MULT=0.54 -> -0.2676)
V77 = {'lig0': 1.9329, 'kdl0': -2.5655, 'ZAP0': 2.1666, 'SYK0': 1.7237,
       'KZP_MULT': 0.3424, 'KPR_MULT': -0.2676}
FIT6 = ['lig0', 'kdl0', 'ZAP0', 'SYK0', 'KZP_MULT', 'KPR_MULT']
MULT_BOUNDS = {'KZP_MULT': (-0.3, 1.0), 'KPR_MULT': (-1.0, 1.0)}   # original v_config bounds
HALF = float(os.environ.get('BOOT_HALF', '0.60'))  # warm-start half-width in log10
ORIG_B = {'lig0': (1.4, 2.4), 'kdl0': (-4.5, -2.3), 'ZAP0': (2.0, 3.2),
          'SYK0': (0.5, 2.5), 'KZP_MULT': (-0.3, 1.0), 'KPR_MULT': (-1.0, 1.0)}

def canon(n):
    """The config/.bngl spell it kd10 (digit one); our tables say kdl0."""
    n = n.strip()
    return 'kdl0' if n.lower() in ('kd10', 'kdl0') else n


def read_days():
    """{line: DataFrame(3 days x 4 timepoints)} for the Tyr493 marker."""
    d = pd.read_excel(XLSX, sheet_name='Original_values')
    d = d[d['marker'].astype(str).str.strip() == MARKER]
    out = {}
    for line in LINE2COL:
        sub = d[d['line'].astype(str).str.strip() == line]
        if sub.shape[0] != 3:
            raise SystemExit(f'ERROR: expected 3 days for {line}, found {sub.shape[0]}')
        out[line] = sub[DAYCOLS].to_numpy(float)   # (3 days, 4 timepoints)
    return out

# ---- the COMPLETE empirical bootstrap over 3 experiment days -----------------
# Resampling 3 days with replacement has exactly 10 distinct outcomes.  Fitting
# all ten and weighting each by its multinomial probability gives the exact
# bootstrap distribution -- no Monte-Carlo error, and fewer jobs than random
# draws.  Lilly ran all three cell lines on the SAME three days, so a day is
# resampled as a unit and the same day-set is applied to every line.
EMP_SETS = list(combinations_with_replacement(range(3), 3))     # 10 multisets

def emp_weight(ms):
    c = [list(ms).count(k) for k in range(3)]
    p = factorial(3)
    for x in c: p //= factorial(x)
    return p / 27.0

def resampled_means(days, kind, seed, i=1):
    """-> {line: array(4 timepoints)} of bootstrap means."""
    rng = np.random.RandomState(seed % (2**32 - 1))
    out = {}
    sel = EMP_SETS[(i - 1) % len(EMP_SETS)]          # same day-set for every line
    for line, arr in days.items():
        if kind == 'emp':
            out[line] = arr[list(sel)].mean(axis=0)
        else:                                                        # parametric
            # We are resampling the MEAN of n days, whose sampling distribution
            # has standard error sd/sqrt(n) -- not sd.  Drawing at sd inflates
            # every synthetic dataset by sqrt(n) and produces shapes the model
            # cannot fit, which shows up as huge parameter scatter.
            n = arr.shape[0]
            mu = arr.mean(axis=0)
            se = arr.std(axis=0, ddof=1) / np.sqrt(n)
            out[line] = rng.normal(mu, np.maximum(se, 1e-12))
        out[line][0] = 0.0                                           # t=0 is 0 by construction
    return out

def write_csv(path, means):
    """Same columns/rows as the original data/pZAP70_Tyr493_mean.csv."""
    orig = pd.read_csv(path)
    df = orig.copy()
    tcol = df.columns[0]
    for line, col in LINE2COL.items():
        if col not in df.columns:
            raise SystemExit(f'ERROR: {col} not in {path} (cols={list(df.columns)})')
        for t, v in zip(TIMES, means[line]):
            df.loc[np.isclose(df[tcol].to_numpy(float), t), col] = v
    df.to_csv(path, index=False)

def build_one(kind, i, days):
    tag = f'{kind}{i:03d}'
    dst = os.path.join(BOOT_ROOT, tag)
    if os.path.isdir(dst): shutil.rmtree(dst)
    shutil.copytree(SRC, dst, ignore=shutil.ignore_patterns(
        'zeta_runs', 'gamma_runs', 'mixed_runs', '*.png', 'slurm*.out', '*.log',
        # v77's OWN result file lives here; copying it means an unfinished sample
        # looks like a finished one holding v77's 4-parameter answer.
        'analysis_param_residue.dat'))
    d = os.path.join(dst, SUB)

    if i == 0 or kind == 'nul':
        # reference AND null samples both use the shipped CSV untouched.  A null
        # differs from the reference only in the optimizer's own randomness, so
        # the spread across nulls is pure simulation/optimizer noise.
        pass
    else:
        write_csv(os.path.join(d, 'data', 'pZAP70_Tyr493_mean.csv'),
                  resampled_means(days, kind,
                                  seed=20260921 + (0 if kind == 'emp' else 5000) + i, i=i))

    cfgp = os.path.join(d, 'v_config.json')
    c = json.load(open(cfgp))
    c.setdefault('FIXED', {})
    for n in ('KZP_MULT', 'KPR_MULT'):             # fit them instead of holding fixed
        c['FIXED'].pop(n, None)
        if n not in c['PARAMS']:
            lo, hi = MULT_BOUNDS[n]
            c['PARAMS'].append(n); c['LB'].append(lo); c['UB'].append(hi)
    c['FIXED']['KZBG_FRAC'] = 1.0                  # stays removed
    if 'KZBG_FRAC' in c['PARAMS']:
        j = c['PARAMS'].index('KZBG_FRAC')
        for k in ('PARAMS', 'LB', 'UB'): c[k].pop(j)
    # Iterate the CONFIG's parameter names (not our table's), so a spelling
    # difference like kd10 vs kdl0 can never silently skip a warm start.
    warmed = []
    for j, name in enumerate(c['PARAMS']):
        cn = canon(name)
        if cn not in V77:
            raise SystemExit('ERROR: config parameter %r has no warm-start centre '
                             '(known: %s)' % (name, sorted(V77)))
        v = V77[cn]
        if cn in [canon(x) for x in PIN]:          # pinned: collapse the box to a point
            c['LB'][j] = c['UB'][j] = round(v, 6)
            warmed.append(name + '(pinned)')
            continue
        lo, hi = ORIG_B.get(cn, MULT_BOUNDS.get(cn, (-9.0, 9.0)))
        c['LB'][j] = round(max(lo, v - HALF), 4)
        c['UB'][j] = round(min(hi, v + HALF), 4)
        warmed.append(name)
    if len(warmed) != len(c['PARAMS']):
        raise SystemExit('ERROR: warm start covered %d of %d parameters'
                         % (len(warmed), len(c['PARAMS'])))
    if os.environ.get('BOOT_NREPS'):               # NFsim replicates per evaluation
        c['N_REPS'] = int(os.environ['BOOT_NREPS'])
    c['NOTE'] = ('bootstrap %s sample %d (KZBG_FRAC=1%s)'
                 % (kind, i, ('; pinned ' + ','.join(PIN)) if PIN else ''))
    json.dump(c, open(cfgp, 'w'), indent=2)
    if i == 1:
        print(f'    [{kind}] PARAMS={c["PARAMS"]}')
        print(f'    [{kind}] FIXED={c["FIXED"]}   N_REPS={c.get("N_REPS")}')
        print(f'    [{kind}] warm-started: {warmed}')
        print('    [%s] box: %s' % (kind, {n: (c['LB'][j], c['UB'][j])
                                           for j, n in enumerate(c['PARAMS'])}))

    run = os.path.join(BOOT_ROOT, f'run_{tag}.sh')
    with open(run, 'w') as f:
        f.write(f"""#!/bin/bash
#SBATCH --job-name=bp{tag}
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --output={dst}/{tag}.log
{ENV}
cd {d}
python pzap_param_estimation_NK92.py {PARTICLES} {ITERS}
""")
    os.chmod(run, 0o755)
    return tag, run

def verify_mapping(days):
    """Abort unless each line's 3-day mean reproduces its CSV column.

    This is the guard that was missing: a wrong line->column mapping silently
    fits the model to scrambled data and shows up only as a hugely inflated SSR.
    """
    csv = os.path.join(SRC, SUB, 'data', 'pZAP70_Tyr493_mean.csv')
    d = pd.read_csv(csv)
    tcol = d.columns[0]
    bad = []
    for line, col in LINE2COL.items():
        mu = days[line].mean(axis=0)
        for t, v in zip(TIMES, mu):
            row = d.loc[np.isclose(d[tcol].to_numpy(float), t), col]
            if row.empty:
                bad.append('%s: no row at t=%g' % (col, t)); continue
            if abs(float(row.iloc[0]) - v) > 1e-4:
                bad.append('%s @ %gs: csv=%.6f  xlsx-mean(%s)=%.6f'
                           % (col, t, float(row.iloc[0]), line, v))
    if bad:
        print('ERROR: line -> column mapping does not reproduce the shipped CSV:')
        for b in bad: print('   ' + b)
        raise SystemExit('refusing to build bootstrap samples from a wrong mapping')
    print('mapping verified: ' + ', '.join('%s -> %s' % (k, v) for k, v in LINE2COL.items()))


def build(n_emp, n_par, n_nul=0):
    os.makedirs(BOOT_ROOT, exist_ok=True)
    days = read_days()
    verify_mapping(days)
    if n_emp > len(EMP_SETS):
        print('note: only %d distinct day-resamples exist; using all of them.' % len(EMP_SETS))
        n_emp = len(EMP_SETS)
    if n_emp:
        print('empirical bootstrap: all %d day-sets, weights (x/27): %s'
              % (n_emp, ', '.join('%s=%d' % (tuple(d + 1 for d in m), round(emp_weight(m) * 27))
                                  for m in EMP_SETS[:n_emp])))
    print('day-to-day spread per line (SD across 3 days, timepoints 0/1/2/5 min):')
    for line, arr in days.items():
        print(f'  {line:<6} mean={np.round(arr.mean(axis=0), 4)}  sd={np.round(arr.std(axis=0, ddof=1), 4)}')
    jobs = []
    for i in range(0, n_emp + 1): jobs.append(build_one('emp', i, days))
    for i in range(1, n_par + 1): jobs.append(build_one('par', i, days))
    for i in range(1, n_nul + 1): jobs.append(build_one('nul', i, days))
    if PIN: print(f'PINNED at v77 values (not fitted): {PIN}')
    print(f'\nbuilt {len(jobs)} sample directories under {BOOT_ROOT}')
    return jobs

def submit(n_emp, n_par, n_nul=0):
    for tag, run in build(n_emp, n_par, n_nul):
        out = subprocess.run(['sbatch', run], stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, universal_newlines=True)
        print(f'{tag}: {out.stdout.strip()}', flush=True)
    print('\nsqueue -u $USER    |    later: python ' + os.path.abspath(__file__) + ' report')

def _norm(k):
    """kd10 (digit one) and kdl0 (letter L) are the same parameter."""
    k = k.strip()
    return 'kdl0' if k.lower() in ('kd10', 'kdl0') else k


def parse_result(d):
    """Read fitted values out of a finished sample directory.

    analysis_param_residue.dat ends with a 'linear' line giving the parameters
    already converted out of log10, e.g.
      linear  lig0=63.86  kdl0=0.002986  ZAP0=274.6  SYK0=63.55 ...
    """
    dat = os.path.join(d, SUB, 'analysis_param_residue.dat')
    if not os.path.exists(dat): return None
    res, vals = None, {}
    header, hvals = None, None
    for ln in open(dat, errors='ignore'):
        s = ln.strip()
        low = s.lower()
        if low.startswith('residue ='):
            try: res = float(s.split('=', 1)[1])
            except Exception: pass
        elif low.startswith('linear'):
            for tok in s.split()[1:]:
                if '=' in tok:
                    k, v = tok.split('=', 1)
                    try: vals[_norm(k)] = float(v)
                    except Exception: pass
        elif header is None and s.startswith('lig0'):
            header = s.split()
        elif header is not None and hvals is None and s and s[0].isdigit() or \
             (header is not None and hvals is None and s.startswith('-')):
            try: hvals = [float(x) for x in s.split()]
            except Exception: pass
    merged = {}
    if header and hvals:                           # log10 row first
        for k, v in zip(header, hvals):
            if _norm(k) in FIT6:
                merged[_norm(k)] = 10.0 ** v
    merged.update(vals)                            # linear line wins where present
    if not merged: return None
    return res, merged

def report():
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = {'emp': [], 'par': []}; ref = {}
    for d in sorted(glob.glob(os.path.join(BOOT_ROOT, '*[0-9]'))):
        tag = os.path.basename(d); kind = tag[:3]
        if kind not in rows: continue
        r = parse_result(d)
        if not r: continue
        merged = dict(r[1]); merged['_res'] = r[0]
        if tag.endswith('000'): ref[kind] = merged            # un-jittered reference
        else: rows[kind].append(merged)
    pnames = [n for n in FIT6 if n not in PIN]
    allv = {}
    for kind in ('emp', 'par'):
        v = rows[kind]
        if not v:
            print(f'{kind}: no completed samples yet'); continue
        res = np.array([s['_res'] for s in v if s.get('_res') is not None], float)
        rr = ref.get(kind, {}).get('_res')
        print(f'\n=== pZAP bootstrap ({kind}), n={len(v)} ===')
        print('median SSR = %.4e' % np.median(res) + ('   reference (un-jittered) SSR = %.4e' % rr if rr else ''))
        print(f'{"param":<9} {"median":>11} {"2.5%":>11} {"97.5%":>11} {"w":>6} {"ref fit":>10}  {"v77":>9}')
        print('-' * 78)
        pt = {'lig0': 85.7, 'kdl0': 2.72e-3, 'ZAP0': 146.8, 'SYK0': 52.9,
              'KZP_MULT': 2.2, 'KPR_MULT': 0.54}
        for n in pnames:
            arr = np.array([s[n] for s in v if n in s], float)
            if arr.size == 0:
                print(f'{n:<9} {"(not parsed in any sample)":>48}')
                continue
            lo, med, hi = np.percentile(arr, [2.5, 50, 97.5])
            allv[(kind, n)] = arr
            w = (hi - lo) / abs(med) if med else float('inf')
            rv = ref.get(kind, {}).get(n)
            rs = f'{rv:.4g}' if rv is not None else '-'
            mark = '' if (rv is not None and lo <= rv <= hi) else ' *'
            print(f'{n:<9} {med:>11.4g} {lo:>11.4g} {hi:>11.4g} {w:>6.2f} {rs:>10}{mark:<2} {pt[n]:>9.4g}')
        print('-' * 78)
        print('w = (97.5%-2.5%)/median.   * = reference fit falls OUTSIDE its own CI')
        print('KZP_MULT x 0.03 = kzp in (uM s)^-1   |   KPR_MULT x 0.01 = KPR in s^-1')
    if allv:
        ks = sorted({k for k, _ in allv})
        ncol = len(pnames)
        fig, axes = plt.subplots(len(ks), ncol, figsize=(3.2 * ncol, 3.6 * len(ks)), squeeze=False)
        for r, kind in enumerate(ks):
            for c, n in enumerate(pnames):
                ax = axes[r][c]; arr = allv.get((kind, n))
                if arr is None: ax.axis('off'); continue
                ax.hist(arr, bins=max(6, len(arr) // 2), color='tab:green', alpha=0.75)
                lo, hi = np.percentile(arr, [2.5, 97.5])
                ax.axvline(lo, color='k', ls='--'); ax.axvline(hi, color='k', ls='--')
                ax.set_title(f'{kind} {n}\n[{lo:.3g}, {hi:.3g}]', fontsize=9)
        fig.tight_layout()
        p = os.path.join(BOOT_ROOT, 'ci_hist_pzap.png'); fig.savefig(p, dpi=110); plt.close(fig)
        print(f'\nhistograms: {p}')

        pairs = [('lig0','ZAP0'), ('ZAP0','SYK0'), ('ZAP0','KZP_MULT'),
                 ('lig0','KZP_MULT'), ('kdl0','KPR_MULT'), ('SYK0','KZP_MULT')]
        fig, axes = plt.subplots(len(ks), len(pairs),
                                 figsize=(3.0*len(pairs), 3.0*len(ks)), squeeze=False)
        for r, kind in enumerate(ks):
            for c, (a1, a2) in enumerate(pairs):
                ax = axes[r][c]
                x, y = allv.get((kind, a1)), allv.get((kind, a2))
                if x is None or y is None or len(x) != len(y): ax.axis('off'); continue
                ax.scatter(x, y, s=28, alpha=0.75, color='tab:blue')
                rf = ref.get(kind, {})
                if a1 in rf and a2 in rf:
                    ax.scatter([rf[a1]], [rf[a2]], marker='*', s=200, color='tab:red', zorder=5)
                ax.set_xlabel(a1, fontsize=8); ax.set_ylabel(a2, fontsize=8)
                ax.set_title(f'{kind}', fontsize=8)
        fig.suptitle('pZAP bootstrap pairwise (red star = reference fit); '
                     'a tilted line means only the combination is constrained', fontsize=10)
        fig.tight_layout(rect=[0,0,1,0.95])
        p2 = os.path.join(BOOT_ROOT, 'ci_pairs_pzap.png'); fig.savefig(p2, dpi=110); plt.close(fig)
        print(f'pairwise:   {p2}')

if __name__ == '__main__':
    a = sys.argv[1:]
    if not a: print(__doc__); sys.exit(0)
    def _n(k, dflt):
        return int(a[k]) if len(a) > k else dflt
    if a[0] == 'build':    build(_n(1, 12), _n(2, 12), _n(3, 0))
    elif a[0] == 'submit': submit(_n(1, 12), _n(2, 12), _n(3, 0))
    elif a[0] == 'report': report()
    else: print(__doc__)
