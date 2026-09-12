#!/usr/bin/env python
"""
check_fit.py (V55 version, N-parameter aware)
Reads optimized log10 params (any number) from analysis_param_residue.dat and FIXED extras
from v_config.json, re-simulates with BioNetGen, and produces the 2-panel plot:
  LEFT  : gamma-normalized  -> absolute ordering hetero>gamma>zeta
  RIGHT : own-normalized    -> shape fit per curve
Prints pure SSR under both normalizations, peak ratios, and the model-vs-exp value at every
experimental time point so misses are explicit.

Usage: python check_fit.py [n_reps]     (default 3)
"""
import sys, os, json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from calculate_SSR_gamma import gamma_PZAP, calculate_residue
from calculate_SSR_zeta import zeta_PZAP
from calculate_SSR_mixed import mixed_PZAP

FIT_TMAX = 300.0
BASE4 = ['lig0', 'kd10', 'ZAP0', 'SYK0']
n = int(sys.argv[1]) if len(sys.argv) > 1 else 3

# ---------- experimental data ----------
data = pd.read_csv('data/pZAP70_Tyr493_mean.csv')
exp_time = data['time'].to_numpy(dtype=float)
mask = exp_time <= FIT_TMAX
t = exp_time[mask]
z = data['mean_zeta'].to_numpy(dtype=float)[mask]
g = data['mean_gamma'].to_numpy(dtype=float)[mask]
m = data['mean_hetero'].to_numpy(dtype=float)[mask]
eg, ez, em = max(g.max(), 1e-12), max(z.max(), 1e-12), max(m.max(), 1e-12)
exp_hg, exp_zg = em / eg, ez / eg

# ---------- config (FIXED extras) ----------
FIXED = {}
if os.path.exists('v_config.json'):
    with open('v_config.json') as f:
        FIXED = json.load(f).get('FIXED', {})

# ---------- read optimized params ----------
outdir = 'analysis'
with open(f'{outdir}_param_residue.dat') as f:
    lines = [ln.strip() for ln in f if ln.strip()]
names, estimates = None, None
for k, ln in enumerate(lines):
    if ln.startswith('lig0\tkd10\tZAP0\tSYK0'):
        names = ln.split('\t')[:-1]                       # drop residue_1_4th
        vals = lines[k + 1].split('\t')
        estimates = np.array([float(v) for v in vals[:len(names)]])
        break
if estimates is None:
    names = BASE4
    block = lines[0].split('[', 1)[1].rsplit(']', 1)[0]
    estimates = np.fromstring(block, sep=' ')[:4]

vals = {nm: float(10 ** v) for nm, v in zip(names, estimates)}
lig0, kd10, ZAP0, SYK0 = (vals[k] for k in BASE4)
print('Using optimized log10 parameters:', dict(zip(names, np.round(estimates, 4))))
print('Linear:', '  '.join(f'{k}={v:.4g}' for k, v in vals.items()), ' FIXED:', FIXED)


def extras_for(kind):
    out = {}
    for src in (FIXED, vals):
        for k, v in src.items():
            if k in BASE4:
                continue
            if ':' in k:
                mk, name = k.split(':', 1)
                if mk == kind:
                    out[name] = v
            else:
                out[k] = v
    return out


# ---------- simulate ----------
tg, pg = gamma_PZAP(n, lig0, kd10, ZAP0, SYK0, outdir, extra=extras_for('gamma'))
tz, pz = zeta_PZAP(n, lig0, kd10, ZAP0, SYK0, outdir, extra=extras_for('zeta'))
tm, pm = mixed_PZAP(n, lig0, kd10, ZAP0, SYK0, outdir, extra=extras_for('mixed'))
pg, pz, pm = (np.asarray(x, dtype=float) for x in (pg, pz, pm))
mg, mz, mm = max(pg.max(), 1e-12), max(pz.max(), 1e-12), max(pm.max(), 1e-12)

# ---------- metrics ----------
r_g_gam = calculate_residue(tg, pg / mg, t, g / eg)
r_z_gam = calculate_residue(tz, pz / mg, t, z / eg)
r_m_gam = calculate_residue(tm, pm / mg, t, m / eg)
ssr_gam = r_g_gam + r_z_gam + r_m_gam
r_g_own = calculate_residue(tg, pg / mg, t, g / eg)
r_z_own = calculate_residue(tz, pz / mz, t, z / ez)
r_m_own = calculate_residue(tm, pm / mm, t, m / em)
ssr_own = r_g_own + r_z_own + r_m_own

print(f'Model peak ratios : hetero/gamma={mm / mg:.4f}  zeta/gamma={mz / mg:.4f}')
print(f'Exp   peak ratios : hetero/gamma={exp_hg:.4f}  zeta/gamma={exp_zg:.4f}')
print(f'[gamma-norm] SSR zeta={r_z_gam:.4e} mixed={r_m_gam:.4e} gamma={r_g_gam:.4e}  '
      f'Total SSR={ssr_gam:.4e}  Total SSR^(1/4)={ssr_gam ** 0.25:.4f}')
print(f'[own-norm]   SSR zeta={r_z_own:.4e} mixed={r_m_own:.4e} gamma={r_g_own:.4e}  '
      f'Total SSR={ssr_own:.4e}  Total SSR^(1/4)={ssr_own ** 0.25:.4f}')


def at(tt, yy, tq):
    o = np.argsort(tt)
    return np.interp(tq, np.asarray(tt)[o], np.asarray(yy)[o])


print('Point-by-point (gamma-normalized): time  model vs exp  (diff)')
for name, tt, yy, ee in (('zeta ', tz, pz / mg, z / eg), ('hetero', tm, pm / mg, m / eg), ('gamma', tg, pg / mg, g / eg)):
    mv = at(tt, yy, t)
    print('  ' + name + '  ' + '  '.join(f't={tq:.0f}: {a:.3f} vs {b:.3f} ({a - b:+.3f})' for tq, a, b in zip(t, mv, ee)))
maxdev = max(abs(at(tt, yy, t) - ee).max() for tt, yy, ee in ((tz, pz / mg, z / eg), (tm, pm / mg, m / eg), (tg, pg / mg, g / eg)))
print(f'Largest single-point miss (gamma-norm): {maxdev:.3f}')

# ---------- plot ----------
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))
C = {'z': '#1f77b4', 'm': '#ff7f0e', 'g': '#2ca02c'}


def panel(ax, pz_, pm_, pg_, ez_, em_, eg_, title, ssr):
    ax.plot(tz, pz_, color=C['z'], lw=2.5, label='CD3z model')
    ax.plot(tm, pm_, color=C['m'], lw=2.5, label='Hetero model')
    ax.plot(tg, pg_, color=C['g'], lw=2.5, label='FceRIg model')
    ax.scatter(t, ez_, color=C['z'], marker='o', s=70, edgecolor='k', zorder=3, label='CD3z exp')
    ax.scatter(t, em_, color=C['m'], marker='s', s=70, edgecolor='k', zorder=3, label='Hetero exp')
    ax.scatter(t, eg_, color=C['g'], marker='^', s=75, edgecolor='k', zorder=3, label='FceRIg exp')
    ax.set_xlabel('Time (s)', fontsize=13)
    ax.set_ylabel('Normalized pZAP70', fontsize=13)
    ax.set_title(title, fontsize=14)
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=9, ncol=2)
    ax.text(0.02, 0.98, f'SSR^(1/4) = {ssr ** 0.25:.4f}', transform=ax.transAxes, fontsize=11,
            va='top', bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))


panel(ax1, pz / mg, pm / mg, pg / mg, z / eg, m / eg, g / eg, 'Gamma-normalized (absolute ordering)', ssr_gam)
panel(ax2, pz / mz, pm / mm, pg / mg, z / ez, m / em, g / eg, 'Own-normalized (shape fit)', ssr_own)
extra_txt = '  '.join(f'{k}={v:.3g}' for k, v in vals.items() if k not in BASE4)
fig.suptitle(f'lig0={lig0:.1f}  kd10={kd10:.2e}  ZAP0={ZAP0:.1f}  SYK0={SYK0:.1f}  {extra_txt}\n'
             f'model hetero/gamma peak={mm / mg:.3f} (exp {exp_hg:.3f})   max point miss={maxdev:.3f}', fontsize=11)
plt.tight_layout()
plt.savefig('plot_N.png', dpi=120)
print('Saved plot_N.png')
