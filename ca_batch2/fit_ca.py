#!/usr/bin/env python
"""
fit_ca.py  -  Ca ODE parameter fit from pZAP input (NK92), configurable version of
Indrani's estimate_ca_params.py. Same ODE, same data handling, same output CSV format.

Differences:
  * bounds are configurable (hers pinned C1 at its lower bound of 10)
  * any subset of [C1, C2, g, k3, k4, nn, be] can be fitted; the rest are held FIXED
      nn = Hill coefficient (hers: 9), be = Ca baseline influx (hers: 0)
  * optimizer: "pso" (pyswarms global search + L-BFGS-B polish) or "lbfgs" (multi-start)
  * optional SE-weighted SSR (chi-square) using the SE_* columns of the Ca data
  * reports which fitted parameters ended within 2% of a bound

Usage:  python fit_ca.py ca_config.json
Config keys (all optional except pzap):
  pzap, ca, t0, params, fixed, bounds, optimizer, particles, iters, starts, seed, weighted, note
"""
import sys, os, json, time
import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.optimize import minimize
from multiprocessing import Pool

ALL = ['C1', 'C2', 'g', 'k3', 'k4', 'nn', 'be']
DEFAULT_FIXED = {'nn': 9.0, 'be': 0.0}
# log10 bounds. Indrani's: C1 (1,4)  C2 (-2,1)  g (-5,-1)  k3 (-3,1)  k4 (-3,1)
DEFAULT_BOUNDS = {'C1': [-2.0, 5.0], 'C2': [-3.0, 2.0], 'g': [-5.0, 0.0], 'k3': [-4.0, 2.0],
                  'k4': [-4.0, 2.0], 'nn': [0.0, 1.1], 'be': [-3.0, 1.0]}

cfg = json.load(open(sys.argv[1])) if len(sys.argv) > 1 else {}
PZAP_PATH = cfg.get('pzap', 'optimized_model_pzap/model_output_pzap.csv')
CA_PATH   = cfg.get('ca', 'ca_data/Ca_NK92.csv')
T0        = float(cfg.get('t0', 30.0))
PARAMS    = cfg.get('params', ['C1', 'C2', 'g', 'k3', 'k4'])
FIXED     = dict(DEFAULT_FIXED); FIXED.update(cfg.get('fixed', {}))
BOUNDS    = dict(DEFAULT_BOUNDS); BOUNDS.update(cfg.get('bounds', {}))
OPT       = cfg.get('optimizer', 'pso')
PARTICLES = int(cfg.get('particles', 32))
ITERS     = int(cfg.get('iters', 100))
STARTS    = int(cfg.get('starts', 48))
SEED      = int(cfg.get('seed', 7))
WEIGHTED  = bool(cfg.get('weighted', False))
NOTE      = cfg.get('note', '')
VE, Z = 25.0, 602.0


# ---------------- ODE (identical to Indrani's calcium(), with nn and be exposed) ----------------
def calcium(y, t, pzap, C1, C2, g, k3, k4, nn, be):
    b, k1, k2 = 0.111, 0.7, 0.7
    s = k2 * k2
    ca, h = y
    F = (pzap ** nn / (pzap ** nn + k3 ** nn)) + (k4 * pzap)
    dca = (C1 * h * F) * (((b * k1) + ca) / (k1 + ca)) - (g * ca) + be
    dh = (C2 * F) * ((s / (s + ca ** 2)) - h)
    return [dca, dh]


def solve_ca(t, pzap, p, ca0, h0=1.0):
    o = np.argsort(t); t = np.asarray(t, float)[o]; pzap = np.asarray(pzap, float)[o]
    y0 = [float(ca0), float(h0)]
    ca = np.empty_like(t); ca[0] = y0[0]
    for i in range(1, len(t)):
        sol = odeint(calcium, y0, [t[i - 1], t[i]], args=(pzap[i], p['C1'], p['C2'], p['g'], p['k3'], p['k4'], p['nn'], p['be']))
        y0 = sol[1]; ca[i] = y0[0]
    return ca


# ---------------- data (mirrors load_data) ----------------
def pick(df, cands):
    for c in cands:
        if c in df.columns: return c
    raise KeyError(f'missing {cands}; have {list(df.columns)}')


def load():
    ca_df = pd.read_csv(CA_PATH); pz_df = pd.read_csv(PZAP_PATH)
    tc = pick(ca_df, ['time_seconds', 'time']); tp = pick(pz_df, ['time', 'time_seconds'])
    data = {}
    for k in ['zeta', 'gamma', 'hetero']:
        t_abs = ca_df[tc].to_numpy(float); sig = ca_df[pick(ca_df, [f'mean_{k}'])].to_numpy(float)
        se = ca_df[f'SE_{k}'].to_numpy(float) if f'SE_{k}' in ca_df.columns else np.ones_like(sig)
        m = t_abs >= T0
        pz_t = pz_df[tp].to_numpy(float); pz_raw = pz_df[pick(pz_df, [f'mean_pZAP_{k}', f'mean_{k}'])].to_numpy(float)
        data[k] = dict(t=t_abs[m], y=sig[m], se=se[m], pz_t=pz_t, pz=pz_raw / (VE * Z))
    tmax = min(d['t'].max() for d in data.values())
    for k, d in data.items():
        keep = d['t'] <= tmax
        d['t'], d['y'], d['se'] = d['t'][keep], d['y'][keep], d['se'][keep]
        d['pz_at_t'] = np.interp(d['t'], d['pz_t'], d['pz'], left=d['pz'][0], right=d['pz'][-1])
        d['w'] = 1.0 / np.maximum(d['se'], 1e-6) ** 2 if WEIGHTED else np.ones_like(d['y'])
    return data


DATA = None


def model_curve(d, p):
    return solve_ca(d['t'], d['pz_at_t'], p, ca0=d['y'][0])


def params_from_x(x):
    p = dict(FIXED); p.update({n: 10.0 ** v for n, v in zip(PARAMS, x)}); return p


def cost(x):
    p = params_from_x(x)
    tot = 0.0
    for k, d in DATA.items():
        try:
            mdl = model_curve(d, p)
        except Exception:
            return 1e30
        if not np.all(np.isfinite(mdl)): return 1e30
        tot += float(np.sum(d['w'] * (d['y'] - mdl) ** 2))
    return tot


def cost_batch(X):
    with Pool(processes=min(len(X), os.cpu_count() or 1)) as pool:
        return np.asarray(pool.map(cost, [np.asarray(r) for r in X]), float)


def main():
    global DATA
    t_start = time.time()
    DATA = load()
    lb = np.array([BOUNDS[n][0] for n in PARAMS]); ub = np.array([BOUNDS[n][1] for n in PARAMS])
    print(f'=== fit_ca  {NOTE}')
    print(f'pzap={PZAP_PATH}  t0={T0}  weighted={WEIGHTED}  optimizer={OPT}')
    print(f'fitted={PARAMS}  fixed={FIXED}')
    print('log10 bounds: ' + '  '.join(f'{n}[{BOUNDS[n][0]},{BOUNDS[n][1]}]' for n in PARAMS))
    rng = np.random.default_rng(SEED)

    if OPT == 'pso':
        from pyswarms.single.global_best import GlobalBestPSO
        np.random.seed(SEED)
        opt = GlobalBestPSO(n_particles=PARTICLES, dimensions=len(PARAMS),
                            options={'c1': 1.5, 'c2': 1.5, 'w': 0.5}, bounds=(lb, ub))
        best_cost, best_x = opt.optimize(cost_batch, iters=ITERS)
        print(f'PSO best cost {best_cost:.6e} at {best_x}')
        res = minimize(cost, x0=best_x, method='L-BFGS-B', bounds=list(zip(lb, ub)), options={'maxiter': 2000, 'ftol': 1e-12})
        if res.fun > best_cost: res.x, res.fun = best_x, best_cost
        print(f'after L-BFGS-B polish: {res.fun:.6e}')
    else:
        best = None
        for s in range(STARTS):
            x0 = rng.uniform(lb, ub)
            r = minimize(cost, x0=x0, method='L-BFGS-B', bounds=list(zip(lb, ub)), options={'maxiter': 2000, 'ftol': 1e-12})
            if best is None or r.fun < best.fun: best = r
            print(f'start {s + 1}/{STARTS}: {r.fun:.4e}  (best {best.fun:.4e})', flush=True)
        res = best

    p = params_from_x(res.x)
    print(f'runtime: {(time.time() - t_start) / 60:.1f} min')
    print(f'Best total SSR = {res.fun:.6e}')
    print('Best parameters (linear):')
    for n, v in zip(PARAMS, res.x):
        span = BOUNDS[n][1] - BOUNDS[n][0]
        flag = '   <-- ON BOUND' if (v - BOUNDS[n][0] < 0.02 * span or BOUNDS[n][1] - v < 0.02 * span) else ''
        print(f'  {n} = {10 ** v:.8e}{flag}')
    for n in FIXED:
        if n not in PARAMS: print(f'  {n} = {FIXED[n]:g}  (fixed)')
    per = {}; rows = []; unw = 0.0
    for k, d in DATA.items():
        mdl = model_curve(d, p)
        per[k] = float(np.sum((d['y'] - mdl) ** 2)); unw += per[k]
        for ta, e, m in zip(d['t'], d['y'], mdl):
            rows.append({'curve': k, 'time_abs_s': float(ta), 'time_scaled_s': float(ta - T0), 'exp_ca': float(e), 'model_ca': float(m), 'residual': float(e - m)})
    print('Per-curve SSR (unweighted):')
    for k in per: print(f'  {k}: {per[k]:.6e}')
    print(f'Total unweighted SSR = {unw:.6e}')
    pd.DataFrame(rows).to_csv('Ca_model.csv', index=False)
    json.dump({'note': NOTE, 'pzap': PZAP_PATH, 'params': PARAMS, 'fixed': FIXED, 'x_log10': list(map(float, res.x)),
               'linear': {n: float(10 ** v) for n, v in zip(PARAMS, res.x)}, 'ssr_objective': float(res.fun),
               'ssr_unweighted': unw, 'per_curve': per, 'weighted': WEIGHTED, 't0': T0}, open('ca_result.json', 'w'), indent=2)
    print('Wrote Ca_model.csv, ca_result.json')


if __name__ == '__main__':
    main()
