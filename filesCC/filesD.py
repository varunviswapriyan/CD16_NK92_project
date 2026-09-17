"""
filesD.py -- overnight Ca campaign: ALL PARAMETERS SHARED, differences between
adaptors ONLY via effective ITAM number (input = pZAP / ITAM count), v77 inputs.
Goal: fix the FceRIg peak amplitude and CD3z rise timing that m01 misses.

    cd ~/Ca_fit_c02
    python ~/CD16_NK92_project/filesCC/filesD.py submit    # all modes, parallel
    python ~/CD16_NK92_project/filesCC/filesD.py report    # table when done

KEY IDEA (Tier 1): her ODE's Ca-dependent feedback terms use k1=k2=0.7 uM
(Xenopus) but the data is a fluorescence signal ~50-200, so both feedback
terms are numerically DEAD (positive feedback pinned at 1, gate target at 0).
Fitting a shared unit-scale S (Ca_uM = signal/S) re-activates HER OWN feedback
biology with one shared parameter. Several modes test this.

Modes (all shared; base = her 5 params C1,C2,g,k3,k4 + per-ITAM input):
  n00_base          m01 reference (nothing extra)
  n01_S             + unit scale S                              [Tier 1]
  n02_k1k2          + k1,k2 refit in signal units               [Tier 1]
  n03_S_nn          + S + Hill coefficient nn fitted            [Tier 1]
  n04_S_k1k2_b      + S + k1,k2,b (all Xenopus constants free)  [Tier 1]
  n05_nn            + nn fitted (currently fixed at 9)          [Tier 2 timing]
  n06_nn_h0         + nn + initial gate h0
  n07_tau           + shared input delay tau                    [Tier 2 timing]
  n08_casc          + IP3-like first-order stage (rate kip)     [Tier 3 cascade]
  n09_casc_S        + cascade + S
  n10_sum_S         input = (pZAP+pSYK)/ITAM, + S               [Tier 4 inputs]
  n11_bound_S       input = bound kinase/ITAM, + S
  n12_alpha_S       input = pZAP/ITAM^alpha (fitted), + S
  n13_S_seW         n01 with SE-weighted objective              [Tier 5 objective]
  n14_S_peakW       n01 with peak-window weighting
  n15_S_nn_casc     + S + nn + cascade                          [Tier 6 combos]
  n16_all           + S + nn + k1,k2 + cascade + alpha (diagnostic ceiling)
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
PZAP_PATH = os.environ.get('PZAP_OVERRIDE', 'optimized_model_pzap/model_output_pzap_v77.csv')
OBS_PATHS = {'pSYK_total': 'optimized_model_pzap/model_output_psyk_v77.csv',
             'bound_ZAP_SYK': 'optimized_model_pzap/model_output_boundzs_v77.csv'}
CA_PATH   = 'ca_data/Ca_NK92.csv'
V_DIR     = os.path.expanduser('~/NK92_fit_v77/estimate_params_pzap_cleaned_up')
RUNS_MAP  = {'zeta': 'zeta_runs', 'gamma': 'gamma_runs', 'hetero': 'mixed_runs'}
K_RECENT  = 12
T0, VE, Z = 30.0, 25.0, 602.0
OUT_DIR   = os.environ.get('OUT_OVERRIDE', 'out_filesD')
ITAMS     = {'zeta': 6.0, 'gamma': 2.0, 'hetero': 4.0}
PEAK_WIN, PEAK_W = (40.0, 130.0), 3.0

LOG_B = {'C1': [-2, 5], 'C2': [-3, 2], 'g': [-5, 0], 'k3': [-4, 2], 'k4': [-4, 2],
         'S': [0, 3], 'k1': [-2, 3], 'k2': [-2, 3], 'kip': [-3, 1]}
LIN_B = {'nn': [1.0, 12.0], 'b': [0.0, 1.0], 'h0': [0.05, 1.0], 'tau': [0.0, 40.0], 'al': [0.0, 2.0]}
DEFAULTS = {'nn': 9.0, 'be': 0.0, 'b': 0.111, 'k1': 0.7, 'k2': 0.7, 'S': 1.0,
            'h0': 1.0, 'tau': 0.0, 'al': 1.0, 'kip': None}

N_STARTS, PARTICLES, ITERS = 3, 40, 120
BASE_SEED = 777
REF = {'m01 (v77 in)': 7.598e4, 'shared baseline': 8.538e4}

CONDS   = ['zeta', 'gamma', 'hetero']
DISPLAY = {'zeta': 'CD3z', 'gamma': 'FceRIg', 'hetero': 'Hetero'}
COLORS  = {'zeta': 'tab:blue', 'gamma': 'tab:green', 'hetero': 'tab:orange'}
ENV = ('module load Miniconda3/4.9.2; '
       'source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh; conda activate CD16_v2')

BASE = ['C1', 'C2', 'g', 'k3', 'k4']
def M(extra=(), src='pzap', casc=False, obj='raw', alpha=False):
    return dict(extra=list(extra), src=src, casc=casc, obj=obj, alpha=alpha)

MODES = {
    'n00_base':      M(),
    'n01_S':         M(['S']),
    'n02_k1k2':      M(['k1', 'k2']),
    'n03_S_nn':      M(['S', 'nn']),
    'n04_S_k1k2_b':  M(['S', 'k1', 'k2', 'b']),
    'n05_nn':        M(['nn']),
    'n06_nn_h0':     M(['nn', 'h0']),
    'n07_tau':       M(['tau']),
    'n08_casc':      M(['kip'], casc=True),
    'n09_casc_S':    M(['kip', 'S'], casc=True),
    'n10_sum_S':     M(['S'], src='sum'),
    'n11_bound_S':   M(['S'], src='bound_ZAP_SYK'),
    'n12_alpha_S':   M(['al', 'S'], alpha=True),
    'n13_S_seW':     M(['S'], obj='se'),
    'n14_S_peakW':   M(['S'], obj='peak'),
    'n15_S_nn_casc': M(['S', 'nn', 'kip'], casc=True),
    'n16_all':       M(['S', 'nn', 'k1', 'k2', 'kip', 'al'], casc=True, alpha=True),
}

# ---------------- observable export (from v77 gdat) ----------------
def _read_gdat(path, obs):
    try:
        names = open(path).readline().lstrip('#').split()
        a = np.loadtxt(path, comments='#')
        if a.ndim != 2 or a.shape[1] != len(names): return None, None
        c = {n: a[:, i] for i, n in enumerate(names)}
        low = {n.lower(): n for n in names}
        t = c.get(low.get('time'))
        if t is None: return None, None
        if obs.lower() in low: return t, c[low[obs.lower()]]
        return None, None
    except Exception:
        return None, None

def export_obs(obs, out):
    if os.path.exists(out): return
    print(f'exporting {obs} -> {out}', flush=True)
    series, tg = {}, None
    for cond, sub in RUNS_MAP.items():
        ad = sorted(glob.glob(os.path.join(V_DIR, sub, 'analysis*')), key=os.path.getmtime)[-K_RECENT:]
        tr = []
        for a in ad:
            for g in glob.glob(os.path.join(a, '**', '*.gdat'), recursive=True):
                t, v = _read_gdat(g, obs)
                if t is not None: tr.append((t, v))
        if not tr: raise FileNotFoundError(f'no gdat with {obs} under {sub}')
        t0 = tr[0][0]; series[cond] = np.mean([np.interp(t0, t, v) for t, v in tr], axis=0); tg = t0
    os.makedirs(os.path.dirname(out), exist_ok=True)
    pd.DataFrame({'time': tg, **{f'mean_{c}': series[c] for c in CONDS}}).to_csv(out, index=False)

# ---------------- data ----------------
def pick(df, cands):
    for c in cands:
        if c in df.columns: return c
    raise KeyError(f'missing {cands}; have {list(df.columns)}')

def load(cfg):
    ca_df = pd.read_csv(CA_PATH); pz_df = pd.read_csv(PZAP_PATH)
    src_df = None
    if cfg['src'] == 'sum':
        export_obs('pSYK_total', OBS_PATHS['pSYK_total']); src_df = pd.read_csv(OBS_PATHS['pSYK_total'])
    elif cfg['src'] in OBS_PATHS:
        export_obs(cfg['src'], OBS_PATHS[cfg['src']]); src_df = pd.read_csv(OBS_PATHS[cfg['src']])
    tc = pick(ca_df, ['time_seconds', 'time']); tp = pick(pz_df, ['time', 'time_seconds'])
    data = {}
    for k in CONDS:
        t_abs = ca_df[tc].to_numpy(float)
        sig = ca_df[pick(ca_df, [f'mean_{k}'])].to_numpy(float)
        se = ca_df[f'SE_{k}'].to_numpy(float) if f'SE_{k}' in ca_df.columns else np.ones_like(sig)
        m = t_abs >= T0
        pz_t = pz_df[tp].to_numpy(float)
        pz = pz_df[pick(pz_df, [f'mean_pZAP_{k}', f'mean_{k}'])].to_numpy(float) / (VE * Z)
        if src_df is not None:
            ts = pick(src_df, ['time', 'time_seconds'])
            other = np.interp(pz_t, src_df[ts].to_numpy(float),
                              src_df[pick(src_df, [f'mean_{k}', f'mean_pSYK_{k}'])].to_numpy(float)) / (VE * Z)
            z = pz + other if cfg['src'] == 'sum' else other
        else:
            z = pz
        data[k] = dict(t=t_abs[m], y=sig[m], se=se[m], in_t=pz_t, in_raw=z, itam=ITAMS[k])
    tmax = min(d['t'].max() for d in data.values())
    for k, d in data.items():
        keep = d['t'] <= tmax
        d['t'], d['y'], d['se'] = d['t'][keep], d['y'][keep], d['se'][keep]
    return data

# ---------------- ODE (her equation + optional unit scale / cascade) ----------------
def rhs(y, t, z, p, casc):
    if casc:
        ca, h, ip = y; zi = ip; dip = p['kip'] * (z - ip)
    else:
        ca, h = y; zi = z
    nn = p['nn']
    F = (zi ** nn / (zi ** nn + p['k3'] ** nn)) + (p['k4'] * zi)
    cs = ca / p['S']                         # signal -> uM for the feedback terms
    fb = ((p['b'] * p['k1']) + cs) / (p['k1'] + cs)
    s = p['k2'] * p['k2']
    dca = (p['C1'] * h * F) * fb - (p['g'] * ca) + p['be']
    dh = (p['C2'] * F) * ((s / (s + cs ** 2)) - h)
    return [dca, dh, dip] if casc else [dca, dh]

def solve(d, p, cfg):
    t = d['t']
    z_curve = d['in_raw'] / (d['itam'] ** p['al'])
    z = np.interp(t - p['tau'], d['in_t'], z_curve, left=z_curve[0], right=z_curve[-1])
    y0 = [float(d['y'][0]), float(p['h0'])] + ([float(z[0])] if cfg['casc'] else [])
    ca = np.empty_like(t); ca[0] = y0[0]
    for i in range(1, len(t)):
        sol = odeint(rhs, y0, [t[i - 1], t[i]], args=(z[i], p, cfg['casc']))
        y0 = sol[1]; ca[i] = y0[0]
    return ca

# ---------------- objective ----------------
DATA = None
G = {'cfg': None, 'pnames': []}

def _params(x):
    p = dict(DEFAULTS)
    for n, xi in zip(G['pnames'], x):
        p[n] = xi if n in LIN_B else 10.0 ** xi
    return p

def _sims(x):
    p = _params(x); out = {}
    for k, d in DATA.items():
        try: mdl = solve(d, p, G['cfg'])
        except Exception: return None
        if not np.all(np.isfinite(mdl)): return None
        out[k] = mdl
    return out

def cost(x):
    sims = _sims(x)
    if sims is None: return 1e30
    tot = 0.0; obj = G['cfg']['obj']
    for k, d in DATA.items():
        r2 = (d['y'] - sims[k]) ** 2
        if obj == 'se': tot += float(np.sum(r2 / np.maximum(d['se'], 1e-6) ** 2))
        elif obj == 'peak':
            w = np.where((d['t'] >= PEAK_WIN[0]) & (d['t'] <= PEAK_WIN[1]), PEAK_W, 1.0)
            tot += float(np.sum(w * r2))
        else: tot += float(np.sum(r2))
    return tot

def raw_ssr(x):
    sims = _sims(x)
    return float('inf') if sims is None else float(sum(np.sum((d['y'] - sims[k]) ** 2) for k, d in DATA.items()))

def cost_batch(X):
    with Pool(processes=min(len(X), os.cpu_count() or 1)) as pool:
        return np.asarray(pool.map(cost, [np.asarray(r) for r in X]), float)

# ---------------- fit / plot / run ----------------
def bounds_for(pn):
    lb = [LIN_B[n][0] if n in LIN_B else LOG_B[n][0] for n in pn]
    ub = [LIN_B[n][1] if n in LIN_B else LOG_B[n][1] for n in pn]
    return np.array(lb, float), np.array(ub, float)

def fit(mode):
    from pyswarms.single.global_best import GlobalBestPSO
    lb, ub = bounds_for(G['pnames'])
    best_x, best_f = None, np.inf
    for s in range(N_STARTS):
        np.random.seed(BASE_SEED + 1000 * s + abs(hash(mode)) % 997)
        opt = GlobalBestPSO(n_particles=PARTICLES, dimensions=len(G['pnames']),
                            options={'c1': 1.5, 'c2': 1.5, 'w': 0.6}, bounds=(lb, ub))
        f, x = opt.optimize(cost_batch, iters=ITERS); x = np.asarray(x)
        r1 = minimize(cost, x0=x, method='L-BFGS-B', bounds=list(zip(lb, ub)), options={'maxiter': 3000, 'ftol': 1e-13})
        x1 = r1.x if r1.fun < f else x
        r2 = minimize(cost, x0=x1, method='Nelder-Mead', options={'maxiter': 5000, 'xatol': 1e-9, 'fatol': 1e-9})
        cand = min([(f, x), (r1.fun, r1.x), (r2.fun, r2.x)], key=lambda q: q[0])
        xf = np.clip(cand[1], lb, ub); ff = cost(xf)
        if ff < best_f: best_x, best_f = xf, ff
        print(f'  [{mode}] start {s+1}/{N_STARTS}: obj={ff:.4e} (best {best_f:.4e})', flush=True)
    return best_x, best_f

def plot(mode, x):
    p = _params(x); ssr = raw_ssr(x)
    fig, (aL, aR) = plt.subplots(1, 2, figsize=(19.2, 7.2))
    for k in CONDS:
        d = DATA[k]
        aL.plot(d['in_t'], d['in_raw'] / (d['itam'] ** p['al']), color=COLORS[k], label=DISPLAY[k])
    aL.set_xlabel('Time (s)'); aL.set_ylabel('effective input (uM per ITAM)')
    aL.set_title(f'Input ({G["cfg"]["src"]}) / ITAM^{p["al"]:.2f}'); aL.legend()
    for k in CONDS:
        d = DATA[k]
        aR.errorbar(d['t'], d['y'], yerr=d['se'], fmt='.', ms=4, color=COLORS[k], alpha=0.6,
                    elinewidth=0.8, capsize=2, label=f'{DISPLAY[k]} exp (mean±SE)')
        aR.plot(d['t'], solve(d, p, G['cfg']), '-', lw=2.5, color=COLORS[k], label=f'{DISPLAY[k]} model')
    aR.set_xlabel('Time (s)'); aR.set_ylabel('Ca signal')
    aR.set_title(f'{mode}   total SSR={ssr:.3e}'); aR.legend(ncol=2, fontsize=8)
    ptxt = '  '.join(f'{n}={p[n]:.3g}' for n in G['pnames'])
    fig.suptitle(f'{mode} (all params shared, per-ITAM input)   {ptxt}', fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    path = os.path.join(OUT_DIR, f'plot_ca_{mode}.png'); fig.savefig(path, dpi=100); plt.close(fig)
    return path, ssr, p

def run_mode(mode):
    global DATA
    cfg = MODES[mode]; G['cfg'] = cfg; G['pnames'] = BASE + cfg['extra']
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f'=== {mode}: params={G["pnames"]} src={cfg["src"]} casc={cfg["casc"]} obj={cfg["obj"]} ===', flush=True)
    t0 = time.time()
    try:
        DATA = load(cfg)
        x, f = fit(mode)
        path, ssr, p = plot(mode, x)
        d = DATA['gamma']; mdl = solve(d, p, cfg); i = int(np.argmax(d['y']))
        gap = (d['y'][i] - mdl[i]) / max(d['se'][i], 1e-9)
        dz = DATA['zeta']; mz = solve(dz, p, cfg); iz = int(np.argmax(dz['y']))
        tshift = float(dz['t'][int(np.argmax(mz))] - dz['t'][iz])
        res = dict(ok=True, mode=mode, n_params=len(G['pnames']), raw_ssr=float(ssr),
                   gamma_peak_gap_in_SE=round(float(gap), 2), zeta_peak_time_shift_s=round(tshift, 1),
                   params={n: float(p[n]) for n in G['pnames']}, plot=path,
                   minutes=round((time.time() - t0) / 60, 1))
        print(f'  DONE SSR={ssr:.4e} gamma gap={gap:+.1f} SE zeta peak shift={tshift:+.0f}s ({res["minutes"]} min)', flush=True)
    except Exception as e:
        res = dict(ok=False, mode=mode, error=repr(e)); print(f'  FAILED: {e!r}', flush=True)
    with open(os.path.join(OUT_DIR, f'results_{mode}.json'), 'w') as f: json.dump(res, f, indent=2)

def submit():
    os.makedirs('logs', exist_ok=True)
    for obs, path in OBS_PATHS.items():
        try: export_obs(obs, path)
        except Exception as e: print(f'WARNING: {obs} export failed: {e!r}')
    script = os.path.abspath(__file__)
    for m in MODES:
        cmd = ['sbatch', '-J', m, '-N1', '-n1', '-c', '32', '--mem=16G', '-t', '8:00:00',
               '-o', f'logs/{m}.log', '--wrap', f'{ENV}; cd {os.getcwd()}; python {script} {m}']
        out = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        print(f'{m}: {out.stdout.strip()}', flush=True)
    print('\nAll submitted.  squeue -u $USER   |   later: python ' + script + ' report')

def report():
    rows = []
    for m in MODES:
        f = os.path.join(OUT_DIR, f'results_{m}.json')
        if not os.path.exists(f): rows.append((m, None, None, None, None, 'running/queued')); continue
        r = json.load(open(f))
        if r.get('ok'): rows.append((m, r['raw_ssr'], r['gamma_peak_gap_in_SE'], r['zeta_peak_time_shift_s'], r['n_params'], 'OK'))
        else: rows.append((m, None, None, None, None, 'FAILED ' + r.get('error', '')[:35]))
    rows.sort(key=lambda r: (r[1] is None, r[1] or 0))
    print(f'{"mode":<15} {"#par":>4} {"raw SSR":>11} {"gamma gap":>10} {"zeta t-shift":>13}  status')
    print('-' * 72)
    for m, s, g, ts, n, st in rows:
        if s is None: print(f'{m:<15} {"-":>4} {"-":>11} {"-":>10} {"-":>13}  {st}')
        else: print(f'{m:<15} {n:>4} {s:>11.4e} {g:>7.1f} SE {ts:>+10.0f} s  {st}')
    print('-' * 72)
    for k, v in REF.items(): print(f'{k:<22} {v:.3e}')

def main():
    a = sys.argv[1:]
    if 'submit' in a: submit(); return
    if 'report' in a: report(); return
    todo = [x for x in a if x in MODES]
    if not todo: print('usage: filesD.py submit | report | <mode>'); print(' '.join(MODES)); return
    for m in todo: run_mode(m)

if __name__ == '__main__':
    main()
