"""
filesCC.py -- Ca variant sweep c10-c20 + auto-MAIN c16
======================================================
Built directly on fit_ca.py (Ca_fit_c02): SAME ODE (Indrani's calcium(),
nn=9 be=0 fixed), SAME segment-wise odeint solver, SAME data loader and
pZAP scaling (VE*Z = 25*602), SAME log10 bounds, SAME PSO(pyswarms)+L-BFGS-B.
Only the variant machinery is new: which params are replicated per condition,
fixed, or how residuals are weighted. Model rules untouched.

PER-CONDITION REPLICATION family (can separate FceRIg from Hetero):
  c10  per-condition C1   (coupling gain)
  c11  per-condition g    (decay/downstream gain)
  c12  per-condition Ca0  (baseline offset control)
  c19  per-condition k3   (Hill threshold = input sensitivity)
  c20  per-condition k4   (linear input coupling)

SHARED-PARAM family (objective / identifiability tweaks):
  c13  k4 fixed @ 0.0618 (c08 value), raw SSR
  c14  condition-balanced SSR
  c15  peak-window weighted (3x, t_abs 40-130 s)
  c17  log-residual SSR (relative error)
  c18  combo: k4 fixed + balanced + peak-window

AUTO-MAIN:
  c16  winner of {c10,c11,c19,c20} + k4 fixed + balanced SSR, extra starts

RUN FROM ~/Ca_fit_c02 so the relative data paths resolve (same as fit_ca.py):
    cd ~/Ca_fit_c02
    python ~/CD16_NK92_project/filesCC/filesCC.py            # all variants
    python ~/CD16_NK92_project/filesCC/filesCC.py c13 c10    # subset
Outputs -> ~/Ca_fit_c02/out_filesCC/  (plot_ca_cXX.png, results.json)
Raw SSR is reported for EVERY variant -> apples-to-apples with c02 = 8.402e4.
"""

import sys, os, json, time
import numpy as np
import pandas as pd
from scipy.integrate import odeint
from scipy.optimize import minimize
from multiprocessing import Pool

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---------------- config (mirrors fit_ca.py defaults) ----------------
PZAP_PATH = 'optimized_model_pzap/model_output_pzap.csv'
CA_PATH   = 'ca_data/Ca_NK92.csv'
T0        = 30.0
VE, Z     = 25.0, 602.0
OUT_DIR   = 'out_filesCC'

ALL_PARAMS = ['C1', 'C2', 'g', 'k3', 'k4']
FIXED_ALWAYS = {'nn': 9.0, 'be': 0.0}
# log10 bounds, identical to fit_ca.py DEFAULT_BOUNDS
BOUNDS = {'C1': [-2.0, 5.0], 'C2': [-3.0, 2.0], 'g': [-5.0, 0.0],
          'k3': [-4.0, 2.0], 'k4': [-4.0, 2.0]}

K4_FIXED   = 0.0618          # from c08
PEAK_WIN   = (40.0, 130.0)   # absolute seconds
PEAK_W     = 3.0
PARTICLES  = 32
ITERS      = 100
SEED       = 7
C02_REF    = 8.402e4

CONDS   = ['zeta', 'gamma', 'hetero']
DISPLAY = {'zeta': 'CD3z', 'gamma': 'FceRIg', 'hetero': 'Hetero'}
COLORS  = {'zeta': 'tab:blue', 'gamma': 'tab:green', 'hetero': 'tab:orange'}

# ---------------- ODE: verbatim from fit_ca.py ----------------
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
        sol = odeint(calcium, y0, [t[i - 1], t[i]],
                     args=(pzap[i], p['C1'], p['C2'], p['g'], p['k3'], p['k4'],
                           p['nn'], p['be']))
        y0 = sol[1]; ca[i] = y0[0]
    return ca

# ---------------- data: verbatim logic from fit_ca.py ----------------
def pick(df, cands):
    for c in cands:
        if c in df.columns: return c
    raise KeyError(f'missing {cands}; have {list(df.columns)}')


def load():
    ca_df = pd.read_csv(CA_PATH); pz_df = pd.read_csv(PZAP_PATH)
    tc = pick(ca_df, ['time_seconds', 'time']); tp = pick(pz_df, ['time', 'time_seconds'])
    data = {}
    for k in CONDS:
        t_abs = ca_df[tc].to_numpy(float); sig = ca_df[pick(ca_df, [f'mean_{k}'])].to_numpy(float)
        se = ca_df[f'SE_{k}'].to_numpy(float) if f'SE_{k}' in ca_df.columns else np.ones_like(sig)
        m = t_abs >= T0
        pz_t = pz_df[tp].to_numpy(float)
        pz_raw = pz_df[pick(pz_df, [f'mean_pZAP_{k}', f'mean_{k}'])].to_numpy(float)
        data[k] = dict(t=t_abs[m], y=sig[m], se=se[m], pz_t=pz_t, pz=pz_raw / (VE * Z))
    tmax = min(d['t'].max() for d in data.values())
    for k, d in data.items():
        keep = d['t'] <= tmax
        d['t'], d['y'], d['se'] = d['t'][keep], d['y'][keep], d['se'][keep]
        d['pz_at_t'] = np.interp(d['t'], d['pz_t'], d['pz'], left=d['pz'][0], right=d['pz'][-1])
    return data

# ---------------- variant machinery ----------------
# (vector entries: shared params once, per-condition params replicated;
#  all model params searched in log10 like fit_ca.py; Ca0 searched linear)
VARIANTS = {
    'c10': dict(label='per-condition C1 (coupling gain)',        per_cond=['C1'], fixed={}, obj='raw'),
    'c11': dict(label='per-condition g (decay gain)',            per_cond=['g'],  fixed={}, obj='raw'),
    'c12': dict(label='per-condition Ca0 (offset control)',      per_cond=[],     fixed={}, obj='raw', ca0_fit=True),
    'c19': dict(label='per-condition k3 (Hill threshold)',       per_cond=['k3'], fixed={}, obj='raw'),
    'c20': dict(label='per-condition k4 (linear coupling)',      per_cond=['k4'], fixed={}, obj='raw'),
    'c13': dict(label=f'shared, k4 fixed @ {K4_FIXED}',          per_cond=[], fixed={'k4': K4_FIXED}, obj='raw'),
    'c14': dict(label='shared, condition-balanced SSR',          per_cond=[], fixed={}, obj='balanced'),
    'c15': dict(label='shared, peak-weighted (3x, 40-130s)',     per_cond=[], fixed={}, obj='peak'),
    'c17': dict(label='shared, log-residual SSR',                per_cond=[], fixed={}, obj='logres'),
    'c18': dict(label='shared combo: k4 fix + balanced + peak',  per_cond=[], fixed={'k4': K4_FIXED}, obj='balanced+peak'),
}
PER_COND_FAM = ['c10', 'c11', 'c12', 'c19', 'c20']
SHARED_FAM   = ['c13', 'c14', 'c15', 'c17', 'c18']
TIER1_POOL   = {'c10': 'C1', 'c11': 'g', 'c19': 'k3', 'c20': 'k4'}


def build_spec(v, data):
    entries = []  # (name, base, cond|None, lo, hi, is_log)
    for n in ALL_PARAMS:
        if n in v['fixed']: continue
        lo, hi = BOUNDS[n]
        if n in v['per_cond']:
            for c in CONDS: entries.append((f'{n}_{c}', n, c, lo, hi, True))
        else:
            entries.append((n, n, None, lo, hi, True))
    if v.get('ca0_fit'):
        for c in CONDS:
            y0 = float(data[c]['y'][0])
            entries.append((f'Ca0_{c}', 'Ca0', c, 0.3 * y0, 3.0 * y0, False))
    return entries


def params_from_x(x, entries, fixed):
    per = {c: dict(FIXED_ALWAYS, **fixed) for c in CONDS}
    for xi, (nm, base, cond, lo, hi, is_log) in zip(x, entries):
        val = 10.0 ** xi if is_log else xi
        if cond is None:
            for c in CONDS: per[c][base] = val
        else:
            per[cond][base] = val
    return per


def describe(x, entries, fixed):
    out = {nm: float(10.0 ** xi if is_log else xi)
           for xi, (nm, _, _, _, _, is_log) in zip(x, entries)}
    out.update({k: float(v) for k, v in fixed.items()})
    return out

# ---- globals for Pool workers (fork carries them, like fit_ca.py) ----
DATA = None
G_ENTRIES = None
G_FIXED = None
G_OBJ = 'raw'


def _sim_all(x):
    per = params_from_x(x, G_ENTRIES, G_FIXED)
    out = {}
    for k, d in DATA.items():
        p = per[k]
        ca0 = p.get('Ca0', d['y'][0])
        try:
            mdl = solve_ca(d['t'], d['pz_at_t'], p, ca0)
        except Exception:
            return None
        if not np.all(np.isfinite(mdl)): return None
        out[k] = mdl
    return out


def cost(x):
    sims = _sim_all(x)
    if sims is None: return 1e30
    tot = 0.0
    for k, d in DATA.items():
        mdl = sims[k]; y = d['y']; t = d['t']
        if G_OBJ == 'logres':
            eps = 1e-9
            tot += float(np.sum((np.log(np.clip(mdl, eps, None)) - np.log(np.clip(y, eps, None))) ** 2))
            continue
        r2 = (y - mdl) ** 2
        if G_OBJ == 'balanced':
            tot += float(np.sum(r2)) / float(np.var(y))
        elif G_OBJ == 'peak':
            w = np.where((t >= PEAK_WIN[0]) & (t <= PEAK_WIN[1]), PEAK_W, 1.0)
            tot += float(np.sum(w * r2))
        elif G_OBJ == 'balanced+peak':
            w = np.where((t >= PEAK_WIN[0]) & (t <= PEAK_WIN[1]), PEAK_W, 1.0)
            tot += float(np.sum(w * r2)) / float(np.var(y))
        else:
            tot += float(np.sum(r2))
    return tot


def raw_ssr(x):
    sims = _sim_all(x)
    if sims is None: return float('inf')
    return float(sum(np.sum((d['y'] - sims[k]) ** 2) for k, d in DATA.items()))


def cost_batch(X):
    with Pool(processes=min(len(X), os.cpu_count() or 1)) as pool:
        return np.asarray(pool.map(cost, [np.asarray(r) for r in X]), float)

# ---------------- fit one variant: PSO + L-BFGS-B (same recipe) ----------------
def fit_variant(cid, v, data, n_pso=1):
    global G_ENTRIES, G_FIXED, G_OBJ
    entries = build_spec(v, data)
    G_ENTRIES, G_FIXED, G_OBJ = entries, v['fixed'], v['obj']
    lb = np.array([e[3] for e in entries]); ub = np.array([e[4] for e in entries])

    from pyswarms.single.global_best import GlobalBestPSO
    best_x, best_f = None, np.inf
    for s in range(n_pso):
        np.random.seed(SEED + s)
        opt = GlobalBestPSO(n_particles=PARTICLES, dimensions=len(entries),
                            options={'c1': 1.5, 'c2': 1.5, 'w': 0.5}, bounds=(lb, ub))
        f, x = opt.optimize(cost_batch, iters=ITERS)
        res = minimize(cost, x0=np.asarray(x), method='L-BFGS-B',
                       bounds=list(zip(lb, ub)), options={'maxiter': 2000, 'ftol': 1e-12})
        xf, ff = (res.x, res.fun) if res.fun < f else (np.asarray(x), f)
        if ff < best_f: best_x, best_f = np.clip(xf, lb, ub), ff
        print(f'  [{cid}] start {s + 1}/{n_pso}: obj={ff:.6e} (best {best_f:.6e})', flush=True)
    return best_x, best_f, raw_ssr(best_x), entries


def plot_variant(cid, v, x, entries, ssr, data):
    per = params_from_x(x, entries, v['fixed'])
    fig, (aL, aR) = plt.subplots(1, 2, figsize=(19.2, 7.2))
    for k in CONDS:
        d = data[k]
        aL.plot(d['pz_t'], d['pz'], color=COLORS[k], label=DISPLAY[k])
    aL.set_xlabel('Time (s)'); aL.set_ylabel('pZAP input (uM)')
    aL.set_title('pZAP input to Ca ODE'); aL.legend()
    for k in CONDS:
        d = data[k]; p = per[k]
        aR.plot(d['t'], d['y'], '.', ms=3, color=COLORS[k], alpha=0.6, label=f'{DISPLAY[k]} exp')
        mdl = solve_ca(d['t'], d['pz_at_t'], p, p.get('Ca0', d['y'][0]))
        aR.plot(d['t'], mdl, '-', lw=2.5, color=COLORS[k], label=f'{DISPLAY[k]} model')
    aR.set_xlabel('Time (s)'); aR.set_ylabel('Ca signal')
    aR.set_title(f'Ca fit   total SSR={ssr:.3e}'); aR.legend(ncol=2, fontsize=8)
    ptxt = '  '.join(f'{k}={v_:.3g}' for k, v_ in describe(x, entries, v['fixed']).items())
    fig.suptitle(f'{cid}: {v["label"]}   {ptxt}', fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    path = os.path.join(OUT_DIR, f'plot_ca_{cid}.png')
    fig.savefig(path, dpi=100); plt.close(fig)
    return path


def main():
    global DATA
    os.makedirs(OUT_DIR, exist_ok=True)
    requested = [a for a in sys.argv[1:] if a.startswith('c')] or None
    print(f'=== filesCC sweep ===\npzap={PZAP_PATH}  ca={CA_PATH}  t0={T0}', flush=True)
    DATA = load()
    results = {}

    def run(cid, v, n_pso=1):
        print(f'\n=== {cid}: {v["label"]} ===', flush=True)
        t0 = time.time()
        try:
            x, obj, ssr, entries = fit_variant(cid, v, DATA, n_pso)
            plot = plot_variant(cid, v, x, entries, ssr, DATA)
            results[cid] = dict(label=v['label'], ok=True, objective=float(obj),
                                raw_ssr=float(ssr), params=describe(x, entries, v['fixed']),
                                x_log10=[float(t) for t in x], plot=plot,
                                minutes=round((time.time() - t0) / 60, 1))
            print(f'  DONE raw SSR = {ssr:.4e}  ({results[cid]["minutes"]} min)', flush=True)
        except Exception as e:
            results[cid] = dict(label=v['label'], ok=False, error=repr(e))
            print(f'  FAILED: {e!r} -- continuing', flush=True)
        with open(os.path.join(OUT_DIR, 'results.json'), 'w') as f:
            json.dump(results, f, indent=2)

    order = ['c13', 'c14', 'c15', 'c17', 'c18', 'c10', 'c11', 'c19', 'c20', 'c12']
    for cid in order:
        if requested and cid not in requested: continue
        run(cid, VARIANTS[cid])

    if requested is None or 'c16' in requested:
        tier1 = {k: results[k]['raw_ssr'] for k in TIER1_POOL
                 if k in results and results[k].get('ok')}
        wp = TIER1_POOL[min(tier1, key=tier1.get)] if tier1 else 'C1'
        if tier1: print(f'\nTier-1 winner param: {wp}', flush=True)
        run('c16', dict(label=f'MAIN: per-cond {wp} + k4 fixed + balanced',
                        per_cond=[wp], fixed={'k4': K4_FIXED}, obj='balanced'),
            n_pso=2)

    print('\n' + '=' * 76)
    print(f'{"id":<5} {"family":<9} {"raw SSR":>12}   label')
    print('-' * 76)
    fam = lambda c: 'per-cond' if c in PER_COND_FAM else ('shared' if c in SHARED_FAM else 'MAIN')
    for cid in sorted(results):
        r = results[cid]
        if r.get('ok'):
            print(f'{cid:<5} {fam(cid):<9} {r["raw_ssr"]:>12.4e}   {r["label"]}')
        else:
            print(f'{cid:<5} {fam(cid):<9} {"--":>12}   FAILED {r["label"]}: {r["error"]}')
    print('-' * 76)
    print(f'{"c02":<5} {"(ref)":<9} {C02_REF:>12.4e}   previous MAIN')
    print(f'\nOutputs in {OUT_DIR}/ (plots + results.json)')


if __name__ == '__main__':
    main()
