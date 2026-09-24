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
import numpy as np
import pandas as pd

HOME     = os.path.expanduser('~')
SRC      = os.path.join(HOME, 'NK92_fit_v77')
SUB      = 'estimate_params_pzap_cleaned_up'
XLSX     = '/home/gddaslab/share/Varun_Indrani/estimate_params_pzap/data/pZAP70_Tyr493_Tyr292_original_and_averages.xlsx'
BOOT_ROOT = os.path.join(HOME, 'boot_pzap')
MARKER   = 'pZAP70 (Tyr493)'
# cell line in the xlsx  ->  column in data/pZAP70_Tyr493_mean.csv
LINE2COL = {'NK92': 'mean_zeta', 'KI 1': 'mean_gamma', 'KI 2': 'mean_hetero'}
DAYCOLS  = ['d_0', 'd_1', 'd_2', 'd_5']          # 0, 60, 120, 300 s -- the fitted points
TIMES    = [0.0, 60.0, 120.0, 300.0]
PARTICLES, ITERS = 16, 20                        # warm-started: short run is enough
ENV = ('module load Miniconda3/4.9.2; '
       'source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh; conda activate CD16_v2')

# warm-start centres (log10). lig0..SYK0 from v77; multipliers from the v69 fit
# that produced the approved values (KZP_MULT=2.2 -> 0.3424, KPR_MULT=0.54 -> -0.2676)
V77 = {'lig0': 1.9329, 'kdl0': -2.5655, 'ZAP0': 2.1666, 'SYK0': 1.7237,
       'KZP_MULT': 0.3424, 'KPR_MULT': -0.2676}
FIT6 = ['lig0', 'kdl0', 'ZAP0', 'SYK0', 'KZP_MULT', 'KPR_MULT']
MULT_BOUNDS = {'KZP_MULT': (-0.3, 1.0), 'KPR_MULT': (-1.0, 1.0)}   # original v_config bounds
HALF = 0.30                                       # warm-start half-width in log10

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

def resampled_means(days, kind, seed):
    """-> {line: array(4 timepoints)} of bootstrap means."""
    rng = np.random.RandomState(seed % (2**32 - 1))
    out = {}
    for line, arr in days.items():
        if kind == 'emp':
            idx = rng.randint(0, arr.shape[0], size=arr.shape[0])    # resample days
            out[line] = arr[idx].mean(axis=0)
        else:                                                        # parametric
            mu = arr.mean(axis=0)
            sd = arr.std(axis=0, ddof=1)
            out[line] = rng.normal(mu, np.maximum(sd, 1e-12))
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
        'zeta_runs', 'gamma_runs', 'mixed_runs', '*.png', 'slurm*.out', '*.log'))
    d = os.path.join(dst, SUB)

    write_csv(os.path.join(d, 'data', 'pZAP70_Tyr493_mean.csv'),
              resampled_means(days, kind, seed=20260921 + (0 if kind == 'emp' else 5000) + i))

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
    for n, v in V77.items():                       # warm start: narrow box, clipped to orig bounds
        if n in c['PARAMS']:
            j = c['PARAMS'].index(n)
            lo, hi = MULT_BOUNDS.get(n, (-9.0, 9.0))
            c['LB'][j] = round(max(lo, v - HALF), 4)
            c['UB'][j] = round(min(hi, v + HALF), 4)
    c['NOTE'] = f'bootstrap {kind} sample {i} (6 fitted params, KZBG_FRAC=1)'
    json.dump(c, open(cfgp, 'w'), indent=2)
    if i == 1:
        print(f'    [{kind}] PARAMS={c["PARAMS"]}')
        print(f'    [{kind}] FIXED={c["FIXED"]}')

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

def build(n_emp, n_par):
    os.makedirs(BOOT_ROOT, exist_ok=True)
    days = read_days()
    print('day-to-day spread per line (SD across 3 days, timepoints 0/1/2/5 min):')
    for line, arr in days.items():
        print(f'  {line:<6} mean={np.round(arr.mean(axis=0), 4)}  sd={np.round(arr.std(axis=0, ddof=1), 4)}')
    jobs = []
    for i in range(1, n_emp + 1): jobs.append(build_one('emp', i, days))
    for i in range(1, n_par + 1): jobs.append(build_one('par', i, days))
    print(f'\nbuilt {len(jobs)} sample directories under {BOOT_ROOT}')
    return jobs

def submit(n_emp, n_par):
    for tag, run in build(n_emp, n_par):
        out = subprocess.run(['sbatch', run], stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, universal_newlines=True)
        print(f'{tag}: {out.stdout.strip()}', flush=True)
    print('\nsqueue -u $USER    |    later: python ' + os.path.abspath(__file__) + ' report')

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
                    try: vals[k.strip()] = float(v)
                    except Exception: pass
        elif header is None and s.startswith('lig0'):
            header = s.split()
        elif header is not None and hvals is None and s and s[0].isdigit() or \
             (header is not None and hvals is None and s.startswith('-')):
            try: hvals = [float(x) for x in s.split()]
            except Exception: pass
    if not vals and header and hvals:              # fall back to the log10 row
        for k, v in zip(header, hvals):
            if k in FIT6:
                vals[k] = 10.0 ** v
    if not vals: return None
    return res, vals

def report():
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    rows = {'emp': [], 'par': []}
    for d in sorted(glob.glob(os.path.join(BOOT_ROOT, '*[0-9]'))):
        tag = os.path.basename(d); kind = tag[:3]
        if kind not in rows: continue
        r = parse_result(d)
        if r:
            merged = dict(r[1]); merged['_res'] = r[0]
            rows[kind].append(merged)
    pnames = FIT6
    allv = {}
    for kind in ('emp', 'par'):
        v = rows[kind]
        if not v:
            print(f'{kind}: no completed samples yet'); continue
        print(f'\n=== pZAP bootstrap ({kind}), n={len(v)} ===')
        print(f'{"param":<6} {"median":>12} {"2.5%":>12} {"97.5%":>12}   v77 point est')
        print('-' * 70)
        pt = {'lig0': 85.7, 'kdl0': 2.72e-3, 'ZAP0': 146.8, 'SYK0': 52.9,
              'KZP_MULT': 2.2, 'KPR_MULT': 0.54}
        for n in pnames:
            arr = np.array([s[n] for s in v if n in s], float)
            if arr.size == 0: continue
            lo, med, hi = np.percentile(arr, [2.5, 50, 97.5])
            allv[(kind, n)] = arr
            print(f'{n:<6} {med:>12.4g} {lo:>12.4g} {hi:>12.4g}   {pt[n]:.4g}')
        print('-' * 70)
        print('KZP_MULT x 0.03 = kzp in (uM s)^-1   |   KPR_MULT x 0.01 = KPR in s^-1')
    if allv:
        ks = sorted({k for k, _ in allv})
        fig, axes = plt.subplots(len(ks), 4, figsize=(16, 4 * len(ks)), squeeze=False)
        for r, kind in enumerate(ks):
            for c, n in enumerate(pnames):
                ax = axes[r][c]; arr = allv.get((kind, n))
                if arr is None: ax.axis('off'); continue
                ax.hist(arr, bins=max(6, len(arr) // 2), color='tab:green', alpha=0.75)
                lo, hi = np.percentile(arr, [2.5, 97.5])
                ax.axvline(lo, color='k', ls='--'); ax.axvline(hi, color='k', ls='--')
                ax.set_title(f'{kind} {n}\n[{lo:.3g}, {hi:.3g}]', fontsize=9)
        fig.tight_layout()
        p = os.path.join(BOOT_ROOT, 'ci_hist_pzap.png'); fig.savefig(p, dpi=110)
        print(f'\nhistograms: {p}')

if __name__ == '__main__':
    a = sys.argv[1:]
    if not a: print(__doc__); sys.exit(0)
    if a[0] == 'build':  build(int(a[1]) if len(a) > 1 else 12, int(a[2]) if len(a) > 2 else 12)
    elif a[0] == 'submit': submit(int(a[1]) if len(a) > 1 else 12, int(a[2]) if len(a) > 2 else 12)
    elif a[0] == 'report': report()
    else: print(__doc__)
