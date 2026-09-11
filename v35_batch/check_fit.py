#!/usr/bin/env python
"""
check_fit.py (V35 version)
Reads optimized log10 params from analysis_param_residue.dat, re-simulates with BioNetGen,
and produces a 2-panel plot:
  LEFT  : gamma-normalized (everything / gamma max)  -> shows ABSOLUTE ordering hetero>gamma>zeta
  RIGHT : own-normalized (each curve / its own max)  -> shows SHAPE fit per curve
Also prints pure SSR under both normalizations (no penalties) and the peak ratios.

Usage: python check_fit.py [n_reps]     (default n_reps = 3)
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from calculate_SSR_gamma import gamma_PZAP, calculate_residue
from calculate_SSR_zeta import zeta_PZAP
from calculate_SSR_mixed import mixed_PZAP

FIT_TMAX = 300.0
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

# ---------- read optimized params ----------
outdir = 'analysis'
with open(f'{outdir}_param_residue.dat') as f:
    lines = [ln.strip() for ln in f if ln.strip()]
estimates = None
for k, ln in enumerate(lines):
    if ln.startswith('lig0\tkd10\tZAP0\tSYK0'):
        vals = lines[k + 1].split('\t')
        estimates = np.array([float(v) for v in vals[:4]])
        break
if estimates is None:
    block = lines[0].split('[', 1)[1].rsplit(']', 1)[0]
    estimates = np.fromstring(block, sep=' ')[:4]

lig0, kd10, ZAP0, SYK0 = 10 ** estimates
print('Using optimized log10 parameters:', estimates)
print(f'lig0={lig0:.4g}, kd10={kd10:.4g}, ZAP0={ZAP0:.4g}, SYK0={SYK0:.4g}')

# ---------- simulate ----------
tg, pg = gamma_PZAP(n, lig0, kd10, ZAP0, SYK0, outdir)
tz, pz = zeta_PZAP(n, lig0, kd10, ZAP0, SYK0, outdir)
tm, pm = mixed_PZAP(n, lig0, kd10, ZAP0, SYK0, outdir)
pg, pz, pm = (np.asarray(x, dtype=float) for x in (pg, pz, pm))
mg, mz, mm = max(pg.max(), 1e-12), max(pz.max(), 1e-12), max(pm.max(), 1e-12)

# ---------- metrics ----------
# gamma-normalized (original objective)
r_g_gam = calculate_residue(tg, pg / mg, t, g / eg)
r_z_gam = calculate_residue(tz, pz / mg, t, z / eg)
r_m_gam = calculate_residue(tm, pm / mg, t, m / eg)
ssr_gam = r_g_gam + r_z_gam + r_m_gam
# own-normalized (shape objective)
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

panel(ax1, pz / mg, pm / mg, pg / mg, z / eg, m / eg, g / eg,
      'Gamma-normalized (absolute ordering)', ssr_gam)
panel(ax2, pz / mz, pm / mm, pg / mg, z / ez, m / em, g / eg,
      'Own-normalized (shape fit)', ssr_own)
fig.suptitle(f'lig0={lig0:.1f}  kd10={kd10:.2e}  ZAP0={ZAP0:.1f}  SYK0={SYK0:.1f}   '
             f'model hetero/gamma peak={mm / mg:.3f} (exp {exp_hg:.3f})', fontsize=12)
plt.tight_layout()
plt.savefig('plot_N.png', dpi=120)
print('Saved plot_N.png')
