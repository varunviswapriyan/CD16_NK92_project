"""
bootstrap_ca2.py -- identifiability campaign for the Ca fit.

The n01 bootstrap gave comparable SSR across every sample while C1 ranged over
three orders of magnitude: the objective has a flat valley because F(pZAP) has
BOTH a Hill term and a linear term, so k3/k4 (and C1) trade off against each
other. Narrow CIs require removing that redundancy, not a different statistic.

This runs the same parametric bootstrap under several reduced parameterisations
and reports, for each, the 95% CI *and its relative width* per parameter next to
the median SSR -- so you can pick the variant that keeps the fit while making
the parameters identifiable.

Variants (all shared across adaptors, per-ITAM input, her ODE):
  full      C1 C2 g k3 k4 S        n01 as sent (baseline, expect wide CIs)
  k3fix     C1 C2 g k4 S           Hill threshold fixed at the n01 value
  k4fix     C1 C2 g k3 S           linear coefficient fixed
  hillonly  C1 C2 g k3 S           k4 = 0: drive is the Hill term alone
  linonly   C1 C2 g k4 S           Hill term removed: drive is linear alone
  Sfix      C1 C2 g k3 k4          S fixed (her 5 params, scale calibrated)
  min4      C1 C2 g k4             k3 and S both fixed
  min4h     C1 C2 g k3             k4 = 0 and S fixed
  lin1      C1 C2 g S              linear drive, k4 = 1 (C1 carries the gain)
  lin1S     C1 C2 g                same, S also fixed

    cd ~/Ca_fit_c02
    python .../bootstrap_ca2.py submit              # all variants, 12 samples each
    python .../bootstrap_ca2.py submit 15 k3fix hillonly linonly
    python .../bootstrap_ca2.py report              # CI + relative-width table
    python .../bootstrap_ca2.py band k3fix          # prediction envelope figure
    python .../bootstrap_ca2.py k3fix 7             # run one sample locally
"""
import os, sys, json, time, glob, subprocess
import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.optimize import minimize
from multiprocessing import Pool
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PZAP_PATH = 'optimized_model_pzap/model_output_pzap_v77.csv'
CA_PATH   = 'ca_data/Ca_NK92.csv'
ROOT      = 'out_boot_ca2'
ITAMS     = {'zeta': 6.0, 'gamma': 2.0, 'hetero': 4.0}
CONDS     = ['zeta', 'gamma', 'hetero']
DISPLAY   = {'zeta': 'CD3z', 'gamma': 'FceRIg', 'hetero': 'Hetero'}
COLORS    = {'zeta': 'tab:blue', 'gamma': 'tab:green', 'hetero': 'tab:orange'}
T0, VE, Z = 30.0, 25.0, 602.0

# n01 point estimate -- the values used when a parameter is held fixed
N01 = {'C1': 1267.0, 'C2': 7.752, 'g': 0.004959, 'k3': 12.81, 'k4': 0.1039, 'S': 31.6}
BASEFIX = {'nn': 9.0, 'be': 0.0, 'b': 0.111, 'k1': 0.7, 'k2': 0.7, 'h0': 1.0}
LOG_B = {'C1': [-2, 5], 'C2': [-3, 2], 'g': [-5, 0], 'k3': [-4, 2], 'k4': [-4, 2], 'S': [0, 3]}

VARIANTS = {
    'full':     dict(fit=['C1','C2','g','k3','k4','S'], fix={},                      hill=True,  lin=True),
    'k3fix':    dict(fit=['C1','C2','g','k4','S'],      fix={'k3': N01['k3']},       hill=True,  lin=True),
    'k4fix':    dict(fit=['C1','C2','g','k3','S'],      fix={'k4': N01['k4']},       hill=True,  lin=True),
    'hillonly': dict(fit=['C1','C2','g','k3','S'],      fix={'k4': 0.0},             hill=True,  lin=False),
    'linonly':  dict(fit=['C1','C2','g','k4','S'],      fix={'k3': N01['k3']},       hill=False, lin=True),
    'Sfix':     dict(fit=['C1','C2','g','k3','k4'],     fix={'S': N01['S']},         hill=True,  lin=True),
    'min4':     dict(fit=['C1','C2','g','k4'],          fix={'k3': N01['k3'], 'S': N01['S']}, hill=True, lin=True),
    'min4h':    dict(fit=['C1','C2','g','k3'],          fix={'k4': 0.0, 'S': N01['S']},       hill=True, lin=False),
    # --- identifiable re-parameterisations of the linear-drive model ---
    # F = k4*z enters as C1*h*k4*z, so C1 and k4 are only ever multiplied.
    # Fix k4 = 1 and let C1 carry the coupling: same model, same curves, no valley.
    'lin1':     dict(fit=['C1','C2','g','S'],           fix={'k4': 1.0},             hill=False, lin=True),
    'lin1S':    dict(fit=['C1','C2','g'],               fix={'k4': 1.0, 'S': 31.5},  hill=False, lin=True),
}

N_STARTS, PARTICLES, ITERS = 2, 30, 90
ENV = ('module load Miniconda3/4.9.2; '
       'source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh; conda activate CD16_v2')

def pick(df, c):
    for x in c:
        if x in df.columns: return x
    raise KeyError(f'missing {c}; have {list(df.columns)}')

def load_base():
    ca = pd.read_csv(CA_PATH); pz = pd.read_csv(PZAP_PATH)
    tc = pick(ca, ['time_seconds','time']); tp = pick(pz, ['time','time_seconds'])
    t_abs = ca[tc].to_numpy(float); m = t_abs >= T0
    base = {}
    for k in CONDS:
        base[k] = dict(t=t_abs[m],
                       mean=ca[pick(ca,[f'mean_{k}'])].to_numpy(float)[m],
                       se=ca[f'SE_{k}'].to_numpy(float)[m],
                       in_t=pz[tp].to_numpy(float),
                       in_raw=pz[pick(pz,[f'mean_pZAP_{k}', f'mean_{k}'])].to_numpy(float)/(VE*Z),
                       itam=ITAMS[k])
    tmax = min(d['t'].max() for d in base.values())
    for d in base.values():
        keep = d['t'] <= tmax
        d['t'], d['mean'], d['se'] = d['t'][keep], d['mean'][keep], d['se'][keep]
    return base

def resample(base, seed, jitter=True):
    rng = np.random.RandomState(seed % (2**32 - 1))
    out = {}
    for k, d in base.items():
        y = rng.normal(d['mean'], np.maximum(d['se'], 1e-9)) if jitter else d['mean'].copy()
        out[k] = dict(d, y=y)
    return out

def rhs(y, t, z, p, cfg):
    ca, h = y
    F = 0.0
    if cfg['hill']: F += z**p['nn'] / (z**p['nn'] + p['k3']**p['nn'])
    if cfg['lin']:  F += p['k4'] * z
    cs = ca / p['S']
    fb = ((p['b']*p['k1']) + cs) / (p['k1'] + cs)
    s = p['k2']*p['k2']
    return [(p['C1']*h*F)*fb - p['g']*ca + p['be'],
            (p['C2']*F)*((s/(s + cs**2)) - h)]

def solve(d, p, cfg):
    t = d['t']; zc = d['in_raw']/d['itam']
    z = np.interp(t, d['in_t'], zc, left=zc[0], right=zc[-1])
    y0 = [float(d['y'][0]), float(p['h0'])]
    ca = np.empty_like(t); ca[0] = y0[0]
    for i in range(1, len(t)):
        y0 = odeint(rhs, y0, [t[i-1], t[i]], args=(z[i], p, cfg))[1]
        ca[i] = y0[0]
    return ca

DATA, CFG, PN = None, None, []

def _params(x):
    p = dict(BASEFIX); p.update(N01); p.update(CFG['fix'])
    for n, xi in zip(PN, x): p[n] = 10.0**xi
    return p

def cost(x):
    p = _params(x); tot = 0.0
    for k, d in DATA.items():
        try: mdl = solve(d, p, CFG)
        except Exception: return 1e30
        if not np.all(np.isfinite(mdl)): return 1e30
        tot += float(np.sum((d['y'] - mdl)**2))
    return tot

def cost_batch(X):
    with Pool(processes=min(len(X), os.cpu_count() or 1)) as pool:
        return np.asarray(pool.map(cost, [np.asarray(r) for r in X]), float)

def fit(seed):
    from pyswarms.single.global_best import GlobalBestPSO
    lb = np.array([LOG_B[n][0] for n in PN], float)
    ub = np.array([LOG_B[n][1] for n in PN], float)
    bx, bf = None, np.inf
    for s in range(N_STARTS):
        np.random.seed((seed*1000 + s) % (2**31 - 1))
        opt = GlobalBestPSO(n_particles=PARTICLES, dimensions=len(PN),
                            options={'c1':1.5,'c2':1.5,'w':0.6}, bounds=(lb,ub))
        f, x = opt.optimize(cost_batch, iters=ITERS); x = np.asarray(x)
        r1 = minimize(cost, x0=x, method='L-BFGS-B', bounds=list(zip(lb,ub)),
                      options={'maxiter':3000,'ftol':1e-13})
        r2 = minimize(cost, x0=(r1.x if r1.fun < f else x), method='Nelder-Mead',
                      options={'maxiter':4000,'xatol':1e-9,'fatol':1e-9})
        c = min([(f,x),(r1.fun,r1.x),(r2.fun,r2.x)], key=lambda q: q[0])
        xf = np.clip(c[1], lb, ub); ff = cost(xf)
        if ff < bf: bx, bf = xf, ff
    return bx, bf

def run_sample(variant, i):
    global DATA, CFG, PN
    CFG = VARIANTS[variant]; PN = CFG['fit']
    out = os.path.join(ROOT, variant); os.makedirs(out, exist_ok=True)
    t0 = time.time()
    print(f'=== {variant} sample {i}: fit={PN} fix={CFG["fix"]} ===', flush=True)
    try:
        DATA = resample(load_base(), seed=20260922 + 97*i, jitter=(i > 0))
        x, f = fit(i if i > 0 else 1)
        p = _params(x)
        res = dict(ok=True, variant=variant, sample=i, ssr=float(f),
                   params={n: float(p[n]) for n in PN},
                   fixed={k: float(v) for k, v in CFG['fix'].items()},
                   minutes=round((time.time()-t0)/60, 1))
        print('  done SSR={:.4e} '.format(f) + ' '.join(f'{n}={p[n]:.4g}' for n in PN), flush=True)
    except Exception as e:
        res = dict(ok=False, variant=variant, sample=i, error=repr(e))
        print(f'  FAILED: {e!r}', flush=True)
    json.dump(res, open(os.path.join(out, f'boot_{i:03d}.json'), 'w'), indent=2)

def submit(n, variants):
    os.makedirs('logs_boot_ca2', exist_ok=True)
    script = os.path.abspath(__file__)
    for v in variants:
        for i in range(0, n + 1):          # sample 0 = un-jittered reference fit
            cmd = ['sbatch', '-J', f'{v[:5]}{i:02d}', '-N1','-n1','-c','32','--mem=16G',
                   '-t','3:00:00','-o', f'logs_boot_ca2/{v}_{i:03d}.log',
                   '--wrap', f'{ENV}; cd {os.getcwd()}; python {script} {v} {i}']
            r = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               universal_newlines=True)
            print(f'{v} {i}: {r.stdout.strip()}', flush=True)
    print(f'\n{len(variants)}x{n+1} jobs submitted.  squeue -u $USER')
    print(f'later: python {script} report')

def collect(variant):
    rows = [json.load(open(f)) for f in sorted(glob.glob(os.path.join(ROOT, variant, 'boot_*.json')))]
    ok = [r for r in rows if r.get('ok')]
    ref = [r for r in ok if r['sample'] == 0]
    boot = [r for r in ok if r['sample'] > 0]
    if boot:
        ss = np.array([r['ssr'] for r in boot], float)
        boot = [r for r, s in zip(boot, ss) if s < 5*np.median(ss)]
    return ref[0] if ref else None, boot, len(rows)

def report():
    print(f'{"variant":<9} {"n":>3} {"med SSR":>10} {"ref SSR":>10}   parameter 95% CIs (relative width)')
    print('-'*110)
    summary = {}
    for v in VARIANTS:
        ref, boot, ntot = collect(v)
        if not boot:
            if ntot: print(f'{v:<9} {"-":>3} {"-":>10} {"-":>10}   ({ntot} files, none converged yet)')
            continue
        ss = np.array([r['ssr'] for r in boot], float)
        parts, widths = [], []
        for n in VARIANTS[v]['fit']:
            a = np.array([r['params'][n] for r in boot], float)
            lo, med, hi = np.percentile(a, [2.5, 50, 97.5])
            w = (hi - lo)/abs(med) if med else np.inf
            widths.append(w)
            parts.append(f'{n}={med:.3g}[{lo:.3g},{hi:.3g}] w={w:.1f}')
        summary[v] = dict(n=len(boot), med_ssr=float(np.median(ss)),
                          ref_ssr=float(ref['ssr']) if ref else None,
                          median_rel_width=float(np.median(widths)),
                          params={n: dict(zip(('lo','med','hi'),
                                  map(float, np.percentile([r['params'][n] for r in boot],[2.5,50,97.5]))))
                                  for n in VARIANTS[v]['fit']})
        rs = f'{ref["ssr"]:.3e}' if ref else '-'
        print(f'{v:<9} {len(boot):>3} {np.median(ss):>10.3e} {rs:>10}   ' + '  '.join(parts))
    print('-'*110)
    print('w = (97.5% - 2.5%) / median.  w < 1 means the parameter is well determined;')
    print('w > 10 means it is not identifiable on its own.\n')
    if summary:
        print(f'{"variant":<9} {"#fit":>4} {"ref SSR":>10} {"median w":>9}   verdict')
        print('-'*62)
        for v, s in sorted(summary.items(), key=lambda kv: kv[1]['median_rel_width']):
            rs = f'{s["ref_ssr"]:.3e}' if s['ref_ssr'] else '-'
            verdict = 'tight CIs' if s['median_rel_width'] < 1 else (
                      'usable' if s['median_rel_width'] < 3 else 'still degenerate')
            print(f'{v:<9} {len(VARIANTS[v]["fit"]):>4} {rs:>10} {s["median_rel_width"]:>9.2f}   {verdict}')
        json.dump(summary, open(os.path.join(ROOT, 'ci_summary.json'), 'w'), indent=2)
        print(f'\nwritten: {ROOT}/ci_summary.json')
        print('pick the variant with the lowest median w whose ref SSR is still acceptable,')
        print('then: python ' + os.path.abspath(__file__) + ' band <variant>')

def band(variant):
    """Prediction envelope: simulate every bootstrap parameter set over the data."""
    global CFG, PN
    CFG = VARIANTS[variant]; PN = CFG['fit']
    ref, boot, _ = collect(variant)
    if not boot: print(f'no results for {variant}'); return
    base = load_base()
    data = {k: dict(d, y=d['mean']) for k, d in base.items()}
    fig, ax = plt.subplots(figsize=(11, 7))
    for k in CONDS:
        d = data[k]
        curves = []
        for r in boot:
            p = dict(BASEFIX); p.update(N01); p.update(CFG['fix']); p.update(r['params'])
            try: curves.append(solve(d, p, CFG))
            except Exception: pass
        if not curves: continue
        C = np.vstack(curves)
        lo, hi = np.percentile(C, [2.5, 97.5], axis=0)
        ax.fill_between(d['t'], lo, hi, color=COLORS[k], alpha=0.25, lw=0)
        ax.plot(d['t'], np.median(C, axis=0), color=COLORS[k], lw=2.5, label=f'{DISPLAY[k]} model')
        ax.errorbar(d['t'], d['mean'], yerr=d['se'], fmt='.', ms=3, color=COLORS[k],
                    alpha=0.5, elinewidth=0.7, capsize=1.5, label=f'{DISPLAY[k]} exp')
    ax.set_xlabel('Time (s)'); ax.set_ylabel('Ca signal')
    ax.set_title(f'{variant}: model prediction with 95% bootstrap band (n={len(boot)})')
    ax.legend(ncol=2, fontsize=9)
    p = os.path.join(ROOT, f'band_{variant}.png')
    fig.tight_layout(); fig.savefig(p, dpi=110); plt.close(fig)
    print(f'wrote {p}')

if __name__ == '__main__':
    a = sys.argv[1:]
    if not a: print(__doc__); sys.exit(0)
    if a[0] == 'submit':
        n = int(a[1]) if len(a) > 1 and a[1].isdigit() else 12
        named = [x for x in a[1:] if not x.isdigit()]
        bad = [x for x in named if x not in VARIANTS]
        if bad:
            print('ERROR: unknown variant(s): ' + ', '.join(bad))
            print('known: ' + ' '.join(VARIANTS))
            print('(is this copy of the script up to date?)')
            sys.exit(1)
        submit(n, named or list(VARIANTS))
    elif a[0] == 'report': report()
    elif a[0] == 'band':   band(a[1] if len(a) > 1 else 'k3fix')
    elif a[0] in VARIANTS: run_sample(a[0], int(a[1]) if len(a) > 1 else 1)
    else: print(__doc__)
