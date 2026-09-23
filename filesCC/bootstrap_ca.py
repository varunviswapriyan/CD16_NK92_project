"""
bootstrap_ca.py -- parametric bootstrap CIs for the Ca fit (n01: her 5 shared
params + unit-scale S, per-ITAM input, v77 pZAP curves).

Each sample: draw a synthetic dataset by sampling Normal(mean, SE) at every
timepoint for all three adaptors, refit n01, record the 6 fitted parameters.
Spread across samples -> 95% CI per parameter (2.5th / 97.5th percentile).

    cd ~/Ca_fit_c02
    python ~/CD16_NK92_project/filesCC/bootstrap_ca.py submit      # 25 jobs
    python ~/CD16_NK92_project/filesCC/bootstrap_ca.py submit 40   # 40 jobs
    python ~/CD16_NK92_project/filesCC/bootstrap_ca.py report      # CI table + histograms
    python ~/CD16_NK92_project/filesCC/bootstrap_ca.py 7           # run sample 7 locally
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
OUT_DIR   = 'out_boot_ca'
ITAMS     = {'zeta': 6.0, 'gamma': 2.0, 'hetero': 4.0}
CONDS     = ['zeta', 'gamma', 'hetero']
T0, VE, Z = 30.0, 25.0, 602.0

# n01 = her 5 shared params + S (fluorescence -> uM scale)
PNAMES = ['C1', 'C2', 'g', 'k3', 'k4', 'S']
LOG_B  = {'C1': [-2, 5], 'C2': [-3, 2], 'g': [-5, 0], 'k3': [-4, 2], 'k4': [-4, 2], 'S': [0, 3]}
FIXED  = {'nn': 9.0, 'be': 0.0, 'b': 0.111, 'k1': 0.7, 'k2': 0.7, 'h0': 1.0, 'al': 1.0, 'tau': 0.0}

N_STARTS, PARTICLES, ITERS = 2, 30, 90
ENV = ('module load Miniconda3/4.9.2; '
       'source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh; conda activate CD16_v2')

def pick(df, cands):
    for c in cands:
        if c in df.columns: return c
    raise KeyError(f'missing {cands}; have {list(df.columns)}')

def load_base():
    """Original mean+SE data and the v77 pZAP inputs."""
    ca = pd.read_csv(CA_PATH); pz = pd.read_csv(PZAP_PATH)
    tc = pick(ca, ['time_seconds', 'time']); tp = pick(pz, ['time', 'time_seconds'])
    base = {}
    t_abs = ca[tc].to_numpy(float); m = t_abs >= T0
    for k in CONDS:
        base[k] = dict(
            t=t_abs[m],
            mean=ca[pick(ca, [f'mean_{k}'])].to_numpy(float)[m],
            se=ca[f'SE_{k}'].to_numpy(float)[m],
            in_t=pz[tp].to_numpy(float),
            in_raw=pz[pick(pz, [f'mean_pZAP_{k}', f'mean_{k}'])].to_numpy(float) / (VE * Z),
            itam=ITAMS[k])
    tmax = min(d['t'].max() for d in base.values())
    for d in base.values():
        keep = d['t'] <= tmax
        d['t'], d['mean'], d['se'] = d['t'][keep], d['mean'][keep], d['se'][keep]
    return base

def resample(base, seed):
    """Parametric bootstrap: y ~ Normal(mean, SE) at each timepoint."""
    rng = np.random.default_rng(seed)
    out = {}
    for k, d in base.items():
        y = rng.normal(d['mean'], np.maximum(d['se'], 1e-9))
        out[k] = dict(d, y=y)
    return out

# ---------------- her ODE (Ca scaled by S inside the feedback terms) ----------------
def rhs(y, t, z, p):
    ca, h = y
    nn = p['nn']
    F = (z ** nn / (z ** nn + p['k3'] ** nn)) + (p['k4'] * z)
    cs = ca / p['S']
    fb = ((p['b'] * p['k1']) + cs) / (p['k1'] + cs)
    s = p['k2'] * p['k2']
    dca = (p['C1'] * h * F) * fb - (p['g'] * ca) + p['be']
    dh = (p['C2'] * F) * ((s / (s + cs ** 2)) - h)
    return [dca, dh]

def solve(d, p):
    t = d['t']
    zc = d['in_raw'] / d['itam']          # per-ITAM input
    z = np.interp(t, d['in_t'], zc, left=zc[0], right=zc[-1])
    y0 = [float(d['y'][0]), float(p['h0'])]
    ca = np.empty_like(t); ca[0] = y0[0]
    for i in range(1, len(t)):
        y0 = odeint(rhs, y0, [t[i - 1], t[i]], args=(z[i], p))[1]
        ca[i] = y0[0]
    return ca

DATA = None
def _params(x):
    p = dict(FIXED)
    for n, xi in zip(PNAMES, x): p[n] = 10.0 ** xi
    return p

def cost(x):
    p = _params(x); tot = 0.0
    for k, d in DATA.items():
        try: mdl = solve(d, p)
        except Exception: return 1e30
        if not np.all(np.isfinite(mdl)): return 1e30
        tot += float(np.sum((d['y'] - mdl) ** 2))
    return tot

def cost_batch(X):
    with Pool(processes=min(len(X), os.cpu_count() or 1)) as pool:
        return np.asarray(pool.map(cost, [np.asarray(r) for r in X]), float)

def fit(seed):
    from pyswarms.single.global_best import GlobalBestPSO
    lb = np.array([LOG_B[n][0] for n in PNAMES], float)
    ub = np.array([LOG_B[n][1] for n in PNAMES], float)
    best_x, best_f = None, np.inf
    for s in range(N_STARTS):
        np.random.seed(seed * 1000 + s)
        opt = GlobalBestPSO(n_particles=PARTICLES, dimensions=len(PNAMES),
                            options={'c1': 1.5, 'c2': 1.5, 'w': 0.6}, bounds=(lb, ub))
        f, x = opt.optimize(cost_batch, iters=ITERS); x = np.asarray(x)
        r1 = minimize(cost, x0=x, method='L-BFGS-B', bounds=list(zip(lb, ub)),
                      options={'maxiter': 3000, 'ftol': 1e-13})
        r2 = minimize(cost, x0=(r1.x if r1.fun < f else x), method='Nelder-Mead',
                      options={'maxiter': 4000, 'xatol': 1e-9, 'fatol': 1e-9})
        cand = min([(f, x), (r1.fun, r1.x), (r2.fun, r2.x)], key=lambda q: q[0])
        xf = np.clip(cand[1], lb, ub); ff = cost(xf)
        if ff < best_f: best_x, best_f = xf, ff
    return best_x, best_f

def run_sample(i):
    global DATA
    os.makedirs(OUT_DIR, exist_ok=True)
    t0 = time.time()
    print(f'=== Ca bootstrap sample {i} ===', flush=True)
    try:
        DATA = resample(load_base(), seed=20260921 + i)
        x, f = fit(i)
        p = _params(x)
        res = dict(ok=True, sample=i, ssr=float(f),
                   params={n: float(p[n]) for n in PNAMES},
                   minutes=round((time.time() - t0) / 60, 1))
        print(f'  done SSR={f:.4e} ' + ' '.join(f'{n}={p[n]:.4g}' for n in PNAMES), flush=True)
    except Exception as e:
        res = dict(ok=False, sample=i, error=repr(e))
        print(f'  FAILED: {e!r}', flush=True)
    with open(os.path.join(OUT_DIR, f'boot_{i:03d}.json'), 'w') as fh:
        json.dump(res, fh, indent=2)

def submit(n):
    os.makedirs('logs_boot_ca', exist_ok=True)
    script = os.path.abspath(__file__)
    for i in range(1, n + 1):
        cmd = ['sbatch', '-J', f'bca{i:03d}', '-N1', '-n1', '-c', '32', '--mem=16G',
               '-t', '3:00:00', '-o', f'logs_boot_ca/boot_{i:03d}.log',
               '--wrap', f'{ENV}; cd {os.getcwd()}; python {script} {i}']
        out = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             universal_newlines=True)
        print(f'sample {i}: {out.stdout.strip()}', flush=True)
    print(f'\n{n} Ca bootstrap jobs submitted.  squeue -u $USER')
    print(f'When done: python {script} report')

def report():
    files = sorted(glob.glob(os.path.join(OUT_DIR, 'boot_*.json')))
    rows = [json.load(open(f)) for f in files]
    ok = [r for r in rows if r.get('ok')]
    if not ok:
        print('no completed samples yet'); return
    # drop fits that clearly failed to converge (SSR far above the median)
    ssr = np.array([r['ssr'] for r in ok], float)
    keep = [r for r, s in zip(ok, ssr) if s < 5 * np.median(ssr)]
    print(f'{len(rows)} files, {len(ok)} converged, {len(keep)} kept after outlier trim\n')
    print(f'{"param":<6} {"median":>12} {"2.5%":>12} {"97.5%":>12}   point est (n01)')
    print('-' * 72)
    POINT = {'C1': None, 'C2': None, 'g': None, 'k3': None, 'k4': None, 'S': None}
    try:
        pe = json.load(open('out_filesD/results_n01_S.json'))['params']
        POINT.update({k: pe.get(k) for k in POINT})
    except Exception:
        pass
    summary = {}
    for n in PNAMES:
        v = np.array([r['params'][n] for r in keep], float)
        lo, med, hi = np.percentile(v, [2.5, 50, 97.5])
        summary[n] = dict(median=float(med), ci_lo=float(lo), ci_hi=float(hi))
        pt = POINT.get(n)
        pts = f'{pt:.4g}' if isinstance(pt, (int, float)) else '-'
        print(f'{n:<6} {med:>12.4g} {lo:>12.4g} {hi:>12.4g}   {pts}')
    print('-' * 72)
    ss = np.array([r['ssr'] for r in keep], float)
    print(f'SSR across samples: median {np.median(ss):.3e}  range [{ss.min():.3e}, {ss.max():.3e}]')
    with open(os.path.join(OUT_DIR, 'ci_summary.json'), 'w') as fh:
        json.dump(dict(n_samples=len(keep), params=summary), fh, indent=2)

    fig, axes = plt.subplots(2, 3, figsize=(15, 8))
    for ax, n in zip(axes.ravel(), PNAMES):
        v = np.array([r['params'][n] for r in keep], float)
        ax.hist(v, bins=max(8, len(v) // 3), color='tab:blue', alpha=0.75)
        lo, hi = np.percentile(v, [2.5, 97.5])
        ax.axvline(lo, color='k', ls='--', lw=1); ax.axvline(hi, color='k', ls='--', lw=1)
        pt = POINT.get(n)
        if isinstance(pt, (int, float)): ax.axvline(pt, color='tab:red', lw=2)
        ax.set_title(f'{n}  [{lo:.3g}, {hi:.3g}]', fontsize=10)
    fig.suptitle(f'Ca (n01) bootstrap parameter distributions, n={len(keep)}'
                 '   dashed = 95% CI, red = point estimate')
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    path = os.path.join(OUT_DIR, 'ci_hist_ca.png'); fig.savefig(path, dpi=110); plt.close(fig)
    print(f'histograms: {path}')

if __name__ == '__main__':
    a = sys.argv[1:]
    if not a: print(__doc__); sys.exit(0)
    if a[0] == 'submit': submit(int(a[1]) if len(a) > 1 else 25)
    elif a[0] == 'report': report()
    else: run_sample(int(a[0]))
