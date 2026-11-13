#!/usr/bin/env python
"""Plot the Ca model against the data for one variant. Read-only; refits nothing.

    python plot_ca_fit.py orig
    python plot_ca_fit.py orig full Sfix --dir /home/gddaslab/vxv016/Ca_fit_c02

Parameters come from out_boot_ca2/<variant>/boot_000.json, which is the reference
fit on un-jittered data. The printed SSR should match the 'ref SSR' column of
`bootstrap_ca2.py report` -- if it does not, the json and the report disagree and
one of them is stale.

The ODE, constants and variant table are copied verbatim from bootstrap_ca2.py so
this does not import that module (importing it runs its CLI dispatch).
"""
import argparse, json, os, sys
import numpy as np, pandas as pd
from scipy.integrate import odeint
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

T0, VE, Z = 30.0, 25.0, 602.0
PZAP_PATH = 'optimized_model_pzap/model_output_pzap_v77.csv'
CA_PATH   = 'ca_data/Ca_NK92.csv'
ROOT      = 'out_boot_ca2'
ITAMS     = {'zeta': 6.0, 'gamma': 2.0, 'hetero': 4.0}
CONDS     = ['zeta', 'gamma', 'hetero']
COLOR     = {'zeta': '#0000CC', 'gamma': '#C00000', 'hetero': '#2BB8C4'}

N01     = {'C1': 1267.0, 'C2': 7.752, 'g': 0.004959, 'k3': 12.81, 'k4': 0.1039, 'S': 31.6}
BASEFIX = {'nn': 9.0, 'be': 0.0, 'b': 0.111, 'k1': 0.7, 'k2': 0.7, 'h0': 1.0}

VARIANTS = {
    'full':     dict(fit=['C1','C2','g','k3','k4','S'], fix={},                             hill=True,  lin=True),
    'k3fix':    dict(fit=['C1','C2','g','k4','S'],      fix={'k3': N01['k3']},              hill=True,  lin=True),
    'k4fix':    dict(fit=['C1','C2','g','k3','S'],      fix={'k4': N01['k4']},              hill=True,  lin=True),
    'hillonly': dict(fit=['C1','C2','g','k3','S'],      fix={'k4': 0.0},                    hill=True,  lin=False),
    'linonly':  dict(fit=['C1','C2','g','k4','S'],      fix={'k3': N01['k3']},              hill=False, lin=True),
    'Sfix':     dict(fit=['C1','C2','g','k3','k4'],     fix={'S': N01['S']},                hill=True,  lin=True),
    'min4':     dict(fit=['C1','C2','g','k4'],          fix={'k3': N01['k3'], 'S': N01['S']}, hill=True, lin=True),
    'min4h':    dict(fit=['C1','C2','g','k3'],          fix={'k4': 0.0, 'S': N01['S']},     hill=True,  lin=False),
    'lin1':     dict(fit=['C1','C2','g','S'],           fix={'k4': 1.0},                    hill=False, lin=True),
    'lin1S':    dict(fit=['C1','C2','g'],               fix={'k4': 1.0, 'S': 31.5},         hill=False, lin=True),
    'orig':     dict(fit=['C1','C2','g','k3','k4'],     fix={'S': 1.0},                     hill=True,  lin=True),
    'orig1':    dict(fit=['C1','C2','g','k3'],          fix={'S': 1.0, 'k4': 1.0},          hill=True,  lin=True),
    'full1':    dict(fit=['C1','C2','g','k3','S'],      fix={'k4': 1.0},                    hill=True,  lin=True),
}


def pick(df, cands):
    for x in cands:
        if x in df.columns:
            return x
    raise KeyError(f'missing one of {cands}; have {list(df.columns)}')


def load_base():
    ca = pd.read_csv(CA_PATH); pz = pd.read_csv(PZAP_PATH)
    tc = pick(ca, ['time_seconds', 'time']); tp = pick(pz, ['time', 'time_seconds'])
    t_abs = ca[tc].to_numpy(float); m = t_abs >= T0
    base = {}
    for k in CONDS:
        base[k] = dict(t=t_abs[m],
                       mean=ca[pick(ca, [f'mean_{k}'])].to_numpy(float)[m],
                       se=ca[f'SE_{k}'].to_numpy(float)[m],
                       in_t=pz[tp].to_numpy(float),
                       in_raw=pz[pick(pz, [f'mean_pZAP_{k}', f'mean_{k}'])].to_numpy(float) / (VE * Z),
                       itam=ITAMS[k])
    tmax = min(d['t'].max() for d in base.values())
    for d in base.values():
        keep = d['t'] <= tmax
        d['t'], d['mean'], d['se'] = d['t'][keep], d['mean'][keep], d['se'][keep]
        d['y'] = d['mean'].copy()
    return base


def rhs(y, t, z, p, cfg):
    ca, h = y
    F = 0.0
    if cfg['hill']: F += z**p['nn'] / (z**p['nn'] + p['k3']**p['nn'])
    if cfg['lin']:  F += p['k4'] * z
    cs = ca / p['S']
    fb = ((p['b'] * p['k1']) + cs) / (p['k1'] + cs)
    s = p['k2'] * p['k2']
    return [(p['C1'] * h * F) * fb - p['g'] * ca + p['be'],
            (p['C2'] * F) * ((s / (s + cs**2)) - h)]


def solve(d, p, cfg):
    t = d['t']; zc = d['in_raw'] / d['itam']
    z = np.interp(t, d['in_t'], zc, left=zc[0], right=zc[-1])
    y0 = [float(d['y'][0]), float(p['h0'])]
    ca = np.empty_like(t); ca[0] = y0[0]
    for i in range(1, len(t)):
        y0 = odeint(rhs, y0, [t[i - 1], t[i]], args=(z[i], p, cfg))[1]
        ca[i] = y0[0]
    return ca


def params_for(variant):
    cfg = VARIANTS[variant]
    p = dict(BASEFIX); p.update(N01); p.update(cfg['fix'])
    src = os.path.join(ROOT, variant, 'boot_000.json')
    if os.path.exists(src):
        d = json.load(open(src))
        p.update(d.get('fixed', {})); p.update(d.get('params', {}))
        return p, cfg, src
    return p, cfg, None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('variants', nargs='*', default=['orig'])
    ap.add_argument('--dir', default='.', help='Ca working dir (holds out_boot_ca2, ca_data, ...)')
    a = ap.parse_args()
    os.chdir(a.dir)

    bad = [v for v in a.variants if v not in VARIANTS]
    if bad:
        sys.exit(f'unknown variant(s) {bad}\nknown: {" ".join(VARIANTS)}')

    base = load_base()
    for v in a.variants:
        p, cfg, src = params_for(v)
        if src is None:
            print(f'[{v}] no {ROOT}/{v}/boot_000.json -- plotting n01 defaults, NOT a fit')
        fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), dpi=130, sharex=True)
        tot = 0.0
        for ax, k in zip(axes, CONDS):
            d = base[k]
            mdl = solve(d, p, cfg)
            tot += float(np.sum((d['y'] - mdl)**2))
            ax.errorbar(d['t'], d['mean'], yerr=d['se'], fmt='o', ms=2.4, lw=0.7,
                        color='0.45', alpha=0.75, label='data (mean $\\pm$ SE)')
            ax.plot(d['t'], mdl, color=COLOR[k], lw=1.9, label='model')
            ax.set_title(f'{k}  ({int(ITAMS[k])} ITAM)', fontsize=10)
            ax.set_xlabel('Time (s)')
            for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
        axes[0].set_ylabel('Indo-1 ratio')
        axes[0].legend(fontsize=7, frameon=False)
        fit = ', '.join(f'{n}={p[n]:.4g}' for n in cfg['fit'])
        fig.suptitle(f'Ca variant "{v}"   SSR = {tot:.4e}   [{fit}]', fontsize=9)
        fig.tight_layout(rect=[0, 0, 1, 0.93])
        out = f'fit_ca_{v}.png'; fig.savefig(out); plt.close(fig)
        print(f'[{v}] SSR = {tot:.4e}   (compare with "ref SSR" in bootstrap_ca2.py report)')
        print(f'      fitted: {fit}')
        print(f'      wrote {os.path.abspath(out)}')


if __name__ == '__main__':
    main()
