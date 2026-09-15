"""
filesB.py -- "make shared params fit" campaign: 13 structural/mechanistic
variants of the Ca fit, ALL with fully shared parameters, run as PARALLEL jobs.

    cd ~/Ca_fit_c02
    python ~/CD16_NK92_project/filesCC/filesB.py submit   # submits all 13
    python ~/CD16_NK92_project/filesCC/filesB.py report   # summary table later

Tier 1 (structural normalization, 0-1 extra shared params):
  m01_itam        input = pZAP / ITAM count (6/2/4)  -- the per-ITAM story
  m02_itam_alpha  input = pZAP / N^alpha, shared fitted alpha in [0,2]
  m12_itam_lin    per-ITAM input + linear coupling (drops Hill: max parsimony)
  m13_itam_iff    per-ITAM input + feedforward inhibition
Tier 2 (nonmonotonic coupling, 1 extra shared param; cf. Das 2010):
  m03_iff         drive = F(pZAP) * Ki/(Ki+pZAP)   (shared Ki)
  m04_bell        drive = F(pZAP) - C5*pZAP, clipped >= 0 (shared C5)
Tier 3 (different BNG inputs, 0 extra params):
  m05_psykb       input = pSYK_bound (receptor-bound pSYK only)
  m06_boundzs     input = bound_ZAP_SYK (all ITAM-bound kinase)
  m07_sum         input = pZAP + pSYK (total phospho-kinase)
  m08_integral    input = cumulative integral of pZAP (exposure, not level)
Tier 4 (simplified ODEs, FEWER params):
  m09_linear      F = k4*pZAP only (no Hill): 4 params
  m10_noh         no h gate: 4 params
  m11_h0          shared fitted initial h0: 6 params

All fits: PSO (40 particles x 120 iters, 3 starts) + L-BFGS-B + Nelder-Mead.
Raw-SSR objective throughout (apples-to-apples with c02 = 8.402e4). Each mode
writes results_{mode}.json + plot_ca_{mode}.png in out_filesB/ and reports the
gamma-peak gap in SE units. Needed BNG observable exports (pSYK_bound,
bound_ZAP_SYK, pSYK_total) are generated automatically from the v69 .gdat
files before submission.
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
OBS_PATHS = {'pSYK_total': 'optimized_model_pzap/model_output_psyk.csv',
             'pSYK_bound': 'optimized_model_pzap/model_output_psykbound.csv',
             'bound_ZAP_SYK': 'optimized_model_pzap/model_output_boundzs.csv'}
CA_PATH   = 'ca_data/Ca_NK92.csv'
V69_DIR   = os.path.expanduser('~/NK92_fit_v69/estimate_params_pzap_cleaned_up')
RUNS_MAP  = {'zeta': 'zeta_runs', 'gamma': 'gamma_runs', 'hetero': 'mixed_runs'}
K_RECENT  = 10
T0, VE, Z = 30.0, 25.0, 602.0
OUT_DIR   = 'out_filesB'
ITAMS     = {'zeta': 6.0, 'gamma': 2.0, 'hetero': 4.0}

LOG_BOUNDS = {'C1': [-2.0, 5.0], 'C2': [-3.0, 2.0], 'g': [-5.0, 0.0],
              'k3': [-4.0, 2.0], 'k4': [-4.0, 2.0],
              'Ki': [-3.0, 1.0], 'C5': [-2.0, 3.0]}
LIN_BOUNDS = {'al': [0.0, 2.0], 'h0': [0.05, 1.0]}
FIXED = {'nn': 9.0, 'be': 0.0}

N_STARTS, PARTICLES, ITERS = 3, 40, 120
BASE_SEED = 4242
C02_REF = 8.402e4

CONDS   = ['zeta', 'gamma', 'hetero']
DISPLAY = {'zeta': 'CD3z', 'gamma': 'FceRIg', 'hetero': 'Hetero'}
COLORS  = {'zeta': 'tab:blue', 'gamma': 'tab:green', 'hetero': 'tab:orange'}

ENV_SETUP = ('module load Miniconda3/4.9.2; '
             'source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh; '
             'conda activate CD16_v2')

# mode = (input source, transform, ode variant, extra params)
MODES = {
    'm01_itam':       dict(src='pzap', tr='itam',    ode='std',    extra=[]),
    'm02_itam_alpha': dict(src='pzap', tr='alpha',   ode='std',    extra=['al']),
    'm03_iff':        dict(src='pzap', tr='none',    ode='iff',    extra=['Ki']),
    'm04_bell':       dict(src='pzap', tr='none',    ode='bell',   extra=['C5']),
    'm05_psykb':      dict(src='pSYK_bound',    tr='none', ode='std', extra=[]),
    'm06_boundzs':    dict(src='bound_ZAP_SYK', tr='none', ode='std', extra=[]),
    'm07_sum':        dict(src='sum',  tr='none',    ode='std',    extra=[]),
    'm08_integral':   dict(src='pzap', tr='integral', ode='std',   extra=[]),
    'm09_linear':     dict(src='pzap', tr='none',    ode='linear', extra=[]),
    'm10_noh':        dict(src='pzap', tr='none',    ode='noh',    extra=[]),
    'm11_h0':         dict(src='pzap', tr='none',    ode='h0',     extra=['h0']),
    'm12_itam_lin':   dict(src='pzap', tr='itam',    ode='linear', extra=[]),
    'm13_itam_iff':   dict(src='pzap', tr='itam',    ode='iff',    extra=['Ki']),
}

def pnames_for(cfg):
    base = ['C1', 'C2', 'g', 'k3', 'k4']
    if cfg['ode'] == 'linear': base = ['C1', 'C2', 'g', 'k4']
    if cfg['ode'] == 'noh':    base = ['C1', 'g', 'k3', 'k4']
    return base + cfg['extra']

# ---------------- observable export from v69 gdat ----------------
def _read_gdat_obs(path, obs):
    try:
        with open(path) as f:
            names = f.readline().lstrip('#').split()
        arr = np.loadtxt(path, comments='#')
        if arr.ndim != 2 or arr.shape[1] != len(names): return None, None
        cols = {n: arr[:, i] for i, n in enumerate(names)}
        if 'time' not in cols or obs not in cols: return None, None
        return cols['time'], cols[obs]
    except Exception:
        return None, None


def export_obs(obs, outpath):
    if os.path.exists(outpath): return
    print(f'Exporting {obs} -> {outpath}', flush=True)
    series, tgrid = {}, None
    for cond, sub in RUNS_MAP.items():
        adirs = sorted(glob.glob(os.path.join(V69_DIR, sub, 'analysis*')),
                       key=os.path.getmtime)[-K_RECENT:]
        traces = []
        for dd in adirs:
            for g in glob.glob(os.path.join(dd, '**', '*.gdat'), recursive=True):
                t, v = _read_gdat_obs(g, obs)
                if t is not None: traces.append((t, v))
        if not traces:
            raise FileNotFoundError(f'no gdat with observable {obs} under {sub}')
        t0 = traces[0][0]
        series[cond] = np.mean([np.interp(t0, t, v) for t, v in traces], axis=0)
        tgrid = t0
    pd.DataFrame({'time': tgrid,
                  **{f'mean_{c}': series[c] for c in CONDS}}).to_csv(outpath, index=False)

# ---------------- data ----------------
def pick(df, cands):
    for c in cands:
        if c in df.columns: return c
    raise KeyError(f'missing {cands}; have {list(df.columns)}')


def load(cfg):
    ca_df = pd.read_csv(CA_PATH)
    pz_df = pd.read_csv(PZAP_PATH)
    src_df = None
    if cfg['src'] in OBS_PATHS:
        export_obs(cfg['src'], OBS_PATHS[cfg['src']])
        src_df = pd.read_csv(OBS_PATHS[cfg['src']])
    if cfg['src'] == 'sum':
        export_obs('pSYK_total', OBS_PATHS['pSYK_total'])
        src_df = pd.read_csv(OBS_PATHS['pSYK_total'])
    tc = pick(ca_df, ['time_seconds', 'time'])
    tp = pick(pz_df, ['time', 'time_seconds'])
    data = {}
    for k in CONDS:
        t_abs = ca_df[tc].to_numpy(float)
        sig = ca_df[pick(ca_df, [f'mean_{k}'])].to_numpy(float)
        se = ca_df[f'SE_{k}'].to_numpy(float) if f'SE_{k}' in ca_df.columns else np.ones_like(sig)
        m = t_abs >= T0
        pz_t = pz_df[tp].to_numpy(float)
        pz = pz_df[pick(pz_df, [f'mean_pZAP_{k}', f'mean_{k}'])].to_numpy(float) / (VE * Z)
        if cfg['src'] == 'pzap':
            z = pz
        elif cfg['src'] == 'sum':
            ts = pick(src_df, ['time', 'time_seconds'])
            ps = np.interp(pz_t, src_df[ts].to_numpy(float),
                           src_df[f'mean_{k}'].to_numpy(float)) / (VE * Z)
            z = pz + ps
        else:
            ts = pick(src_df, ['time', 'time_seconds'])
            z = np.interp(pz_t, src_df[ts].to_numpy(float),
                          src_df[f'mean_{k}'].to_numpy(float)) / (VE * Z)
        # transform
        if cfg['tr'] == 'itam':
            z = z / ITAMS[k]
        elif cfg['tr'] == 'integral':
            dt = np.gradient(pz_t)
            z = np.cumsum(z * dt) / 300.0
        # ('alpha' transform is parameter-dependent -> applied in solve via p)
        d = dict(t=t_abs[m], y=sig[m], se=se[m], in_t=pz_t, in_raw=z, itam=ITAMS[k])
        data[k] = d
    tmax = min(d['t'].max() for d in data.values())
    for k, d in data.items():
        keep = d['t'] <= tmax
        d['t'], d['y'], d['se'] = d['t'][keep], d['y'][keep], d['se'][keep]
        d['z_at_t'] = np.interp(d['t'], d['in_t'], d['in_raw'],
                                left=d['in_raw'][0], right=d['in_raw'][-1])
    return data

# ---------------- ODE ----------------
def rhs(y, t, z, p, ode):
    b, k1, k2 = 0.111, 0.7, 0.7
    s = k2 * k2
    ca, h = y
    nn = p['nn']
    if ode == 'linear':
        F = p['k4'] * z
    else:
        F = (z ** nn / (z ** nn + p['k3'] ** nn)) + (p['k4'] * z)
    if ode == 'iff':
        F = F * (p['Ki'] / (p['Ki'] + z))
    elif ode == 'bell':
        F = max(F - p['C5'] * z, 0.0)
    hh = 1.0 if ode == 'noh' else h
    dca = (p['C1'] * hh * F) * (((b * k1) + ca) / (k1 + ca)) - (p['g'] * ca) + p['be']
    dh = 0.0 if ode == 'noh' else (p['C2'] * F) * ((s / (s + ca ** 2)) - h)
    return [dca, dh]


def solve(d, p, ode, tr):
    t = d['t']; z = d['z_at_t'].copy()
    if tr == 'alpha':
        z = z * (ITAMS_REF ** 0)  # no-op guard
        z = d['z_at_t'] / (d['itam'] ** p['al'])
    h0 = p.get('h0', 1.0)
    y0 = [float(d['y'][0]), float(h0)]
    ca = np.empty_like(t); ca[0] = y0[0]
    for i in range(1, len(t)):
        sol = odeint(rhs, y0, [t[i - 1], t[i]], args=(z[i], p, ode))
        y0 = sol[1]; ca[i] = y0[0]
    return ca

ITAMS_REF = 1.0

# ---------------- objective ----------------
DATA = None
G = {'cfg': None, 'pnames': []}


def _params(x):
    p = dict(FIXED)
    for n, xi in zip(G['pnames'], x):
        p[n] = xi if n in LIN_BOUNDS else 10.0 ** xi
    return p


def _sims(x):
    p = _params(x)
    out = {}
    for k, d in DATA.items():
        try:
            mdl = solve(d, p, G['cfg']['ode'], G['cfg']['tr'])
        except Exception:
            return None
        if not np.all(np.isfinite(mdl)): return None
        out[k] = mdl
    return out


def cost(x):
    sims = _sims(x)
    if sims is None: return 1e30
    return float(sum(np.sum((d['y'] - sims[k]) ** 2) for k, d in DATA.items()))


def cost_batch(X):
    with Pool(processes=min(len(X), os.cpu_count() or 1)) as pool:
        return np.asarray(pool.map(cost, [np.asarray(r) for r in X]), float)

# ---------------- fit + plot ----------------
def bounds_for(pn):
    lb, ub = [], []
    for n in pn:
        if n in LIN_BOUNDS:
            lb.append(LIN_BOUNDS[n][0]); ub.append(LIN_BOUNDS[n][1])
        else:
            lb.append(LOG_BOUNDS[n][0]); ub.append(LOG_BOUNDS[n][1])
    return np.array(lb), np.array(ub)


def fit(mode):
    from pyswarms.single.global_best import GlobalBestPSO
    lb, ub = bounds_for(G['pnames'])
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
        cand = min([(f, x), (r1.fun, r1.x), (r2.fun, r2.x)], key=lambda q: q[0])
        xf = np.clip(cand[1], lb, ub); ff = cost(xf)
        if ff < best_f: best_x, best_f = xf, ff
        print(f'  [{mode}] start {s+1}/{N_STARTS}: SSR={ff:.4e} (best {best_f:.4e})', flush=True)
    return best_x, best_f


def plot(mode, x):
    p = _params(x)
    ssr = cost(x)
    fig, (aL, aR) = plt.subplots(1, 2, figsize=(19.2, 7.2))
    for k in CONDS:
        d = DATA[k]
        z = d['in_raw'] if G['cfg']['tr'] != 'alpha' else d['in_raw'] / (d['itam'] ** p.get('al', 0))
        aL.plot(d['in_t'], z, color=COLORS[k], label=DISPLAY[k])
    aL.set_xlabel('Time (s)'); aL.set_ylabel('effective input (uM)')
    aL.set_title(f'Input to Ca ODE ({G["cfg"]["src"]}, tr={G["cfg"]["tr"]})'); aL.legend()
    for k in CONDS:
        d = DATA[k]
        aR.errorbar(d['t'], d['y'], yerr=d['se'], fmt='.', ms=4, color=COLORS[k],
                    alpha=0.6, elinewidth=0.8, capsize=2, label=f'{DISPLAY[k]} exp (mean±SE)')
        mdl = solve(d, p, G['cfg']['ode'], G['cfg']['tr'])
        aR.plot(d['t'], mdl, '-', lw=2.5, color=COLORS[k], label=f'{DISPLAY[k]} model')
    aR.set_xlabel('Time (s)'); aR.set_ylabel('Ca signal')
    aR.set_title(f'{mode}   total SSR={ssr:.3e}'); aR.legend(ncol=2, fontsize=8)
    ptxt = '  '.join(f'{n}={p[n]:.3g}' for n in G['pnames'])
    fig.suptitle(f'{mode} (all params shared)   {ptxt}', fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    path = os.path.join(OUT_DIR, f'plot_ca_{mode}.png')
    fig.savefig(path, dpi=100); plt.close(fig)
    return path, ssr, p


def run_mode(mode):
    global DATA
    cfg = MODES[mode]
    G['cfg'] = cfg
    G['pnames'] = pnames_for(cfg)
    os.makedirs(OUT_DIR, exist_ok=True)
    print(f'=== {mode}: src={cfg["src"]} tr={cfg["tr"]} ode={cfg["ode"]} '
          f'params={G["pnames"]} ===', flush=True)
    t0 = time.time()
    try:
        DATA = load(cfg)
        x, f = fit(mode)
        path, ssr, p = plot(mode, x)
        d = DATA['gamma']
        mdl = solve(d, p, cfg['ode'], cfg['tr'])
        i = int(np.argmax(d['y']))
        gap = (d['y'][i] - mdl[i]) / max(d['se'][i], 1e-9)
        res = dict(ok=True, mode=mode, n_params=len(G['pnames']),
                   raw_ssr=float(ssr), gamma_peak_gap_in_SE=round(float(gap), 2),
                   params={n: float(p[n]) for n in G['pnames']}, plot=path,
                   minutes=round((time.time() - t0) / 60, 1))
        print(f'  DONE SSR={ssr:.4e} gamma gap={gap:+.1f} SE ({res["minutes"]} min)', flush=True)
    except Exception as e:
        res = dict(ok=False, mode=mode, error=repr(e))
        print(f'  FAILED: {e!r}', flush=True)
    with open(os.path.join(OUT_DIR, f'results_{mode}.json'), 'w') as fjs:
        json.dump(res, fjs, indent=2)


def submit():
    os.makedirs('logs', exist_ok=True)
    # generate all needed observable exports up-front (avoids job races)
    for obs, path in OBS_PATHS.items():
        try:
            export_obs(obs, path)
        except Exception as e:
            print(f'WARNING: could not export {obs}: {e!r} (dependent modes will fail)')
    script = os.path.abspath(__file__)
    for m in MODES:
        cmd = ['sbatch', '-J', m, '-N1', '-n1', '-c', '32', '--mem=16G',
               '-t', '4:00:00', '-o', f'logs/{m}.log',
               '--wrap', f'{ENV_SETUP}; cd {os.getcwd()}; python {script} {m}']
        out = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                             universal_newlines=True)
        print(f'{m}: {out.stdout.strip()}', flush=True)
    print('\nAll submitted. squeue -u $USER to watch.')
    print('When done:  python ' + script + ' report')


def report():
    rows = []
    for m in MODES:
        f = os.path.join(OUT_DIR, f'results_{m}.json')
        if not os.path.exists(f):
            rows.append((m, None, None, None, 'no result yet')); continue
        r = json.load(open(f))
        if r.get('ok'):
            rows.append((m, r['raw_ssr'], r['gamma_peak_gap_in_SE'], r['n_params'], 'OK'))
        else:
            rows.append((m, None, None, None, 'FAILED ' + r.get('error', '')[:40]))
    rows.sort(key=lambda r: (r[1] is None, r[1] if r[1] is not None else 0))
    print(f'{"mode":<15} {"#par":>4} {"raw SSR":>12} {"gamma gap":>10}   status')
    print('-' * 64)
    for m, ssr, gap, npar, st in rows:
        if ssr is None:
            print(f'{m:<15} {"-":>4} {"-":>12} {"-":>10}   {st}')
        else:
            print(f'{m:<15} {npar:>4} {ssr:>12.4e} {gap:>7.1f} SE   {st}')
    print('-' * 64)
    print(f'{"c02 ref":<15} {"5":>4} {C02_REF:>12.4e}      +7 SE   single-kinase shared')


def main():
    args = sys.argv[1:]
    if 'submit' in args: submit(); return
    if 'report' in args: report(); return
    todo = [a for a in args if a in MODES]
    if not todo:
        print('usage: filesB.py submit | report | <mode> [...modes]')
        print('modes:', ' '.join(MODES)); return
    for m in todo:
        run_mode(m)


if __name__ == '__main__':
    main()
