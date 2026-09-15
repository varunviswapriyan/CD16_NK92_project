"""
filesC.py -- Ca fitting, shared parameters, four modes, fully self-contained.

ONE COMMAND runs everything as four PARALLEL slurm jobs:

    cd ~/Ca_fit_c02
    python ~/CD16_NK92_project/filesCC/filesC.py submit

That submits raw / se / syk / syk_se as separate 32-core jobs (each ~15-25 min,
all running at once). The pSYK input needed by the syk modes is generated
AUTOMATICALLY from the v69 BioNetGen .gdat outputs -- no manual export step.

Modes (can also be run directly, e.g. `python filesC.py raw`):
  raw     shared 5-param fit, plain SSR              (Indrani's requested config)
  se      shared 5-param fit, SE-weighted (chi-sq)
  syk     TWO-KINASE model (her JI structure): Ca driven by pZAP AND pSYK,
          all params shared (6 params: C1, C1S, C2, g, k3, k4)
  syk_se  two-kinase + SE-weighted

Each mode: 4 PSO starts (48 particles x 150 iters) + L-BFGS-B + Nelder-Mead,
best kept; plots mean +/- SE; writes its own results_{mode}.json (no clashes).
Summary reports raw SSR (vs c02 = 8.402e4, Indrani = 8.9e4) and the gamma-peak
gap in SE units -- the number that says whether the fit is actually good.
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

# ---------------- config ----------------
PZAP_PATH = 'optimized_model_pzap/model_output_pzap.csv'
PSYK_PATH = 'optimized_model_pzap/model_output_psyk.csv'
CA_PATH   = 'ca_data/Ca_NK92.csv'
V69_DIR   = os.path.expanduser('~/NK92_fit_v69/estimate_params_pzap_cleaned_up')
RUNS_MAP  = {'zeta': 'zeta_runs', 'gamma': 'gamma_runs', 'hetero': 'mixed_runs'}
K_RECENT  = 10          # average the K most recent analysis dirs (best-fit sims)
T0        = 30.0
VE, Z     = 25.0, 602.0
OUT_DIR   = 'out_filesC'

LOG_BOUNDS = {'C1': [-2.0, 5.0], 'C1S': [-2.0, 5.0], 'C2': [-3.0, 2.0],
              'g': [-5.0, 0.0], 'k3': [-4.0, 2.0], 'k4': [-4.0, 2.0]}
FIXED = {'nn': 9.0, 'be': 0.0}

N_STARTS, PARTICLES, ITERS = 4, 48, 150
BASE_SEED = 20260915
C02_REF, INDRANI_REF = 8.402e4, 8.9e4

CONDS   = ['zeta', 'gamma', 'hetero']
DISPLAY = {'zeta': 'CD3z', 'gamma': 'FceRIg', 'hetero': 'Hetero'}
COLORS  = {'zeta': 'tab:blue', 'gamma': 'tab:green', 'hetero': 'tab:orange'}

MODES = {
    'raw':    dict(two_kinase=False, weighted=False, label='shared, raw SSR'),
    'se':     dict(two_kinase=False, weighted=True,  label='shared, SE-weighted'),
    'syk':    dict(two_kinase=True,  weighted=False, label='two-kinase (pZAP+pSYK), shared, raw SSR'),
    'syk_se': dict(two_kinase=True,  weighted=True,  label='two-kinase (pZAP+pSYK), shared, SE-weighted'),
}

ENV_SETUP = ('module load Miniconda3/4.9.2; '
             'source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh; '
             'conda activate CD16_v2')

# ---------------- auto-generate pSYK export from v69 .gdat files ----------------
def _read_gdat_psyk(path):
    try:
        with open(path) as f:
            header = f.readline()
        names = header.lstrip('#').split()
        arr = np.loadtxt(path, comments='#')
        if arr.ndim != 2 or arr.shape[1] != len(names): return None, None
        cols = {n: arr[:, i] for i, n in enumerate(names)}
        t = cols.get('time')
        if t is None: return None, None
        if 'pSYK_total' in cols:
            return t, cols['pSYK_total']
        if 'pSYK_bound' in cols and 'pSYK_free' in cols:
            return t, cols['pSYK_bound'] + cols['pSYK_free']
        hit = next((n for n in names if 'pSYK' in n), None)
        return (t, cols[hit]) if hit else (None, None)
    except Exception:
        return None, None


def export_psyk():
    print(f'Generating {PSYK_PATH} from v69 .gdat files...', flush=True)
    series, tgrid = {}, None
    for cond, sub in RUNS_MAP.items():
        adirs = sorted(glob.glob(os.path.join(V69_DIR, sub, 'analysis*')),
                       key=os.path.getmtime)[-K_RECENT:]
        traces = []
        for dd in adirs:
            for g in glob.glob(os.path.join(dd, '**', '*.gdat'), recursive=True):
                t, v = _read_gdat_psyk(g)
                if t is not None: traces.append((t, v))
        if not traces:
            raise FileNotFoundError(
                f'no usable .gdat with a pSYK observable under {V69_DIR}/{sub}/analysis*')
        t0 = traces[0][0]
        vals = [np.interp(t0, t, v) for t, v in traces]
        series[cond] = np.mean(vals, axis=0); tgrid = t0
        print(f'  {cond}: averaged {len(traces)} gdat trace(s) from {len(adirs)} analysis dir(s)', flush=True)
    os.makedirs(os.path.dirname(PSYK_PATH), exist_ok=True)
    pd.DataFrame({'time': tgrid,
                  **{f'mean_pSYK_{c}': series[c] for c in CONDS}}).to_csv(PSYK_PATH, index=False)
    print(f'  wrote {PSYK_PATH}', flush=True)

# ---------------- ODEs ----------------
def calcium(y, t, pzap, C1, C2, g, k3, k4, nn, be):
    # verbatim from fit_ca.py
    b, k1, k2 = 0.111, 0.7, 0.7
    s = k2 * k2
    ca, h = y
    F = (pzap ** nn / (pzap ** nn + k3 ** nn)) + (k4 * pzap)
    dca = (C1 * h * F) * (((b * k1) + ca) / (k1 + ca)) - (g * ca) + be
    dh = (C2 * F) * ((s / (s + ca ** 2)) - h)
    return [dca, dh]


def calcium2(y, t, pzap, psyk, C1, C1S, C2, g, k3, k4, nn, be):
    # two-kinase extension, same structure (her JI model: ZAP + SYK arms)
    b, k1, k2 = 0.111, 0.7, 0.7
    s = k2 * k2
    ca, h = y
    FZ = (pzap ** nn / (pzap ** nn + k3 ** nn)) + (k4 * pzap)
    FS = (psyk ** nn / (psyk ** nn + k3 ** nn)) + (k4 * psyk)
    dca = (h * (C1 * FZ + C1S * FS)) * (((b * k1) + ca) / (k1 + ca)) - (g * ca) + be
    dh = (C2 * (FZ + FS)) * ((s / (s + ca ** 2)) - h)
    return [dca, dh]


def solve(t, pz, ps, p, ca0, two_kinase, h0=1.0):
    o = np.argsort(t); t = np.asarray(t, float)[o]
    pz = np.asarray(pz, float)[o]
    ps = np.asarray(ps, float)[o] if ps is not None else None
    y0 = [float(ca0), float(h0)]
    ca = np.empty_like(t); ca[0] = y0[0]
    for i in range(1, len(t)):
        if two_kinase:
            sol = odeint(calcium2, y0, [t[i-1], t[i]],
                         args=(pz[i], ps[i], p['C1'], p['C1S'], p['C2'], p['g'],
                               p['k3'], p['k4'], p['nn'], p['be']))
        else:
            sol = odeint(calcium, y0, [t[i-1], t[i]],
                         args=(pz[i], p['C1'], p['C2'], p['g'], p['k3'], p['k4'],
                               p['nn'], p['be']))
        y0 = sol[1]; ca[i] = y0[0]
    return ca

# ---------------- data ----------------
def pick(df, cands):
    for c in cands:
        if c in df.columns: return c
    raise KeyError(f'missing {cands}; have {list(df.columns)}')


def load(need_syk):
    if need_syk and not os.path.exists(PSYK_PATH):
        export_psyk()
    ca_df = pd.read_csv(CA_PATH); pz_df = pd.read_csv(PZAP_PATH)
    ps_df = pd.read_csv(PSYK_PATH) if need_syk else None
    tc = pick(ca_df, ['time_seconds', 'time']); tp = pick(pz_df, ['time', 'time_seconds'])
    data = {}
    for k in CONDS:
        t_abs = ca_df[tc].to_numpy(float)
        sig = ca_df[pick(ca_df, [f'mean_{k}'])].to_numpy(float)
        se = ca_df[f'SE_{k}'].to_numpy(float) if f'SE_{k}' in ca_df.columns else np.ones_like(sig)
        m = t_abs >= T0
        pz_t = pz_df[tp].to_numpy(float)
        pz = pz_df[pick(pz_df, [f'mean_pZAP_{k}', f'mean_{k}'])].to_numpy(float) / (VE * Z)
        d = dict(t=t_abs[m], y=sig[m], se=se[m], pz_t=pz_t, pz=pz)
        if ps_df is not None:
            tps = pick(ps_df, ['time', 'time_seconds'])
            d['ps_t'] = ps_df[tps].to_numpy(float)
            d['ps'] = ps_df[pick(ps_df, [f'mean_pSYK_{k}', f'mean_{k}'])].to_numpy(float) / (VE * Z)
        data[k] = d
    tmax = min(d['t'].max() for d in data.values())
    for k, d in data.items():
        keep = d['t'] <= tmax
        d['t'], d['y'], d['se'] = d['t'][keep], d['y'][keep], d['se'][keep]
        d['pz_at_t'] = np.interp(d['t'], d['pz_t'], d['pz'], left=d['pz'][0], right=d['pz'][-1])
        if 'ps' in d:
            d['ps_at_t'] = np.interp(d['t'], d['ps_t'], d['ps'], left=d['ps'][0], right=d['ps'][-1])
    return data

# ---------------- objective ----------------
DATA = None
G = {'two_kinase': False, 'weighted': False, 'pnames': []}


def _params(x):
    p = dict(FIXED)
    for name, xi in zip(G['pnames'], x): p[name] = 10.0 ** xi
    return p


def _sims(x):
    p = _params(x)
    out = {}
    for k, d in DATA.items():
        try:
            mdl = solve(d['t'], d['pz_at_t'], d.get('ps_at_t'), p, d['y'][0], G['two_kinase'])
        except Exception:
            return None
        if not np.all(np.isfinite(mdl)): return None
        out[k] = mdl
    return out


def cost(x):
    sims = _sims(x)
    if sims is None: return 1e30
    tot = 0.0
    for k, d in DATA.items():
        r2 = (d['y'] - sims[k]) ** 2
        if G['weighted']:
            tot += float(np.sum(r2 / np.maximum(d['se'], 1e-6) ** 2))
        else:
            tot += float(np.sum(r2))
    return tot


def raw_ssr(x):
    sims = _sims(x)
    if sims is None: return float('inf')
    return float(sum(np.sum((d['y'] - sims[k]) ** 2) for k, d in DATA.items()))


def cost_batch(X):
    with Pool(processes=min(len(X), os.cpu_count() or 1)) as pool:
        return np.asarray(pool.map(cost, [np.asarray(r) for r in X]), float)

# ---------------- fit + plot ----------------
def fit(mode_cfg, mode):
    from pyswarms.single.global_best import GlobalBestPSO
    G['two_kinase'] = mode_cfg['two_kinase']
    G['weighted'] = mode_cfg['weighted']
    G['pnames'] = (['C1', 'C1S', 'C2', 'g', 'k3', 'k4'] if mode_cfg['two_kinase']
                   else ['C1', 'C2', 'g', 'k3', 'k4'])
    lb = np.array([LOG_BOUNDS[n][0] for n in G['pnames']])
    ub = np.array([LOG_BOUNDS[n][1] for n in G['pnames']])
    best_x, best_f = None, np.inf
    for s in range(N_STARTS):
        np.random.seed(BASE_SEED + 1000 * s + abs(hash(mode)) % 997)
        opt = GlobalBestPSO(n_particles=PARTICLES, dimensions=len(G['pnames']),
                            options={'c1': 1.5, 'c2': 1.5, 'w': 0.6}, bounds=(lb, ub))
        f, x = opt.optimize(cost_batch, iters=ITERS)
        x = np.asarray(x)
        r1 = minimize(cost, x0=x, method='L-BFGS-B', bounds=list(zip(lb, ub)),
                      options={'maxiter': 3000, 'ftol': 1e-13})
        x1 = r1.x if r1.fun < f else x
        r2 = minimize(cost, x0=x1, method='Nelder-Mead',
                      options={'maxiter': 5000, 'xatol': 1e-9, 'fatol': 1e-9})
        cand = min([(f, x), (r1.fun, r1.x), (r2.fun, r2.x)], key=lambda z: z[0])
        xf = np.clip(cand[1], lb, ub); ff = cost(xf)
        if ff < best_f: best_x, best_f = xf, ff
        print(f'  [{mode}] start {s+1}/{N_STARTS}: obj={ff:.6e} (best {best_f:.6e})', flush=True)
    return best_x, best_f


def plot(mode, mode_cfg, x, data):
    p = _params(x)
    ssr = raw_ssr(x)
    fig, (aL, aR) = plt.subplots(1, 2, figsize=(19.2, 7.2))
    for k in CONDS:
        d = data[k]
        aL.plot(d['pz_t'], d['pz'], color=COLORS[k], label=f'{DISPLAY[k]} pZAP')
        if 'ps' in d and mode_cfg['two_kinase']:
            aL.plot(d['ps_t'], d['ps'], '--', color=COLORS[k], alpha=0.7,
                    label=f'{DISPLAY[k]} pSYK')
    aL.set_xlabel('Time (s)'); aL.set_ylabel('kinase input (uM)')
    aL.set_title('Inputs to Ca ODE'); aL.legend(fontsize=8)
    for k in CONDS:
        d = data[k]
        aR.errorbar(d['t'], d['y'], yerr=d['se'], fmt='.', ms=4, color=COLORS[k],
                    alpha=0.6, elinewidth=0.8, capsize=2,
                    label=f'{DISPLAY[k]} exp (mean±SE)')
        mdl = solve(d['t'], d['pz_at_t'], d.get('ps_at_t'), p, d['y'][0],
                    mode_cfg['two_kinase'])
        aR.plot(d['t'], mdl, '-', lw=2.5, color=COLORS[k], label=f'{DISPLAY[k]} model')
    aR.set_xlabel('Time (s)'); aR.set_ylabel('Ca signal')
    aR.set_title(f'{mode_cfg["label"]}   total SSR={ssr:.3e}')
    aR.legend(ncol=2, fontsize=8)
    ptxt = '  '.join(f'{n}={p[n]:.3g}' for n in G['pnames'])
    fig.suptitle(f'{mode}: {mode_cfg["label"]}   {ptxt}', fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    path = os.path.join(OUT_DIR, f'plot_ca_{mode}.png')
    fig.savefig(path, dpi=100); plt.close(fig)
    return path, ssr, p

# ---------------- self-submitting parallel jobs ----------------
def submit():
    script = os.path.abspath(__file__)
    os.makedirs('logs', exist_ok=True)
    # generate the pSYK export ONCE up front so parallel jobs don't race
    if not os.path.exists(PSYK_PATH):
        export_psyk()
    for m in MODES:
        cmd = ['sbatch', '-J', f'fc_{m}', '-N1', '-n1', '-c', '32', '--mem=16G',
               '-t', '4:00:00', '-o', f'logs/fc_{m}.log',
               '--wrap', f'{ENV_SETUP}; cd {os.getcwd()}; python {script} {m}']
        out = subprocess.run(cmd, capture_output=True, text=True)
        print(f'{m}: {out.stdout.strip() or out.stderr.strip()}', flush=True)
    print('\nAll four submitted. Check:  squeue -u $USER')
    print('Progress:                  tail -f logs/fc_raw.log')
    print('When done:                 grep "DONE\\|FAILED" logs/fc_*.log')


def main():
    global DATA
    args = sys.argv[1:]
    if 'submit' in args:
        submit(); return
    os.makedirs(OUT_DIR, exist_ok=True)
    req = [a for a in args if a in MODES] or ['raw', 'se']
    need_syk = any(MODES[m]['two_kinase'] for m in req)
    print(f'modes: {req}', flush=True)
    DATA = load(need_syk)
    summary = {}
    for mode in req:
        cfg = MODES[mode]
        print(f'\n=== {mode}: {cfg["label"]} ===', flush=True)
        t0 = time.time()
        try:
            x, f = fit(cfg, mode)
            path, ssr, p = plot(mode, cfg, x, DATA)
            d = DATA['gamma']
            mdl = solve(d['t'], d['pz_at_t'], d.get('ps_at_t'), p, d['y'][0], cfg['two_kinase'])
            i = int(np.argmax(d['y']))
            gap_se = (d['y'][i] - mdl[i]) / max(d['se'][i], 1e-9)
            res = dict(label=cfg['label'], ok=True, objective=float(f),
                       raw_ssr=float(ssr), params={n: float(p[n]) for n in G['pnames']},
                       gamma_peak_gap_in_SE=round(float(gap_se), 2), plot=path,
                       minutes=round((time.time() - t0) / 60, 1))
            print(f'  DONE raw SSR={ssr:.4e}  gamma peak gap={gap_se:+.1f} SE', flush=True)
        except Exception as e:
            res = dict(label=cfg['label'], ok=False, error=repr(e))
            print(f'  FAILED: {e!r}', flush=True)
        summary[mode] = res
        with open(os.path.join(OUT_DIR, f'results_{mode}.json'), 'w') as fjs:
            json.dump(res, fjs, indent=2)

    print('\n' + '=' * 72)
    print(f'{"mode":<8} {"raw SSR":>12} {"gamma gap":>11}   label')
    print('-' * 72)
    for m, r in summary.items():
        if r.get('ok'):
            print(f'{m:<8} {r["raw_ssr"]:>12.4e} {r["gamma_peak_gap_in_SE"]:>8.1f} SE   {r["label"]}')
        else:
            print(f'{m:<8} {"--":>12} {"--":>11}   FAILED {r["error"]}')
    print('-' * 72)
    print(f'{"c02":<8} {C02_REF:>12.4e}              previous best shared')
    print(f'{"Indrani":<8} {INDRANI_REF:>12.4e}              her logged shared fit')


if __name__ == '__main__':
    main()
