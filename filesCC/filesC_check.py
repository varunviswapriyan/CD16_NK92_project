"""
filesC_check.py -- run AFTER the four filesC jobs. One command, two things:

  1. pSYK sanity check: is gamma's pSYK genuinely above hetero's (signal),
     or within simulation noise? Prints peak / sustained values per condition
     and the gamma/hetero ratios.
  2. Seeded polish: takes the syk_se solution and polishes it under the RAW
     objective, giving the definitive two-kinase raw-SSR number and a plot
     (out_filesC/plot_ca_syk_seeded.png). Fixes the stuck fc_syk run without
     rerunning PSO.

Usage (login node is fine, takes ~2-5 min):
    cd ~/Ca_fit_c02
    python ~/CD16_NK92_project/filesCC/filesC_check.py
"""

import os, sys, json
import numpy as np
import pandas as pd
from scipy.optimize import minimize

sys.path.insert(0, os.path.expanduser('~/CD16_NK92_project/filesCC'))
import filesC as F

# ---------- 1. pSYK sanity check ----------
print('=' * 60)
print('pSYK sanity check (from model_output_psyk.csv)')
print('=' * 60)
d = pd.read_csv(F.PSYK_PATH)
tcol = 'time' if 'time' in d.columns else d.columns[0]
stats = {}
for c in ['zeta', 'gamma', 'hetero']:
    s = d[f'mean_pSYK_{c}']
    stats[c] = dict(peak=float(s.max()),
                    sustained=float(s[d[tcol] > 100].mean()))
    print(f'  {c:<7} peak={stats[c]["peak"]:.4f}   mean(t>100s)={stats[c]["sustained"]:.4f}')
rp = stats['gamma']['peak'] / stats['hetero']['peak']
rs = stats['gamma']['sustained'] / stats['hetero']['sustained']
print(f'  gamma/hetero: peak ratio={rp:.3f}   sustained ratio={rs:.3f}')
if rp > 1.0 and rs > 1.0:
    print('  -> gamma pSYK consistently ABOVE hetero: mechanism is grounded.')
elif rp < 1.0 and rs < 1.0:
    print('  -> gamma pSYK BELOW hetero: separation is NOT coming from pSYK')
    print('     ordering -- investigate before emailing.')
else:
    print('  -> mixed ordering (peak vs sustained differ): small effect,')
    print('     likely near simulation noise -- treat cautiously.')

# ---------- 2. seeded raw-objective polish of the two-kinase model ----------
print()
print('=' * 60)
print('Seeded two-kinase polish under RAW objective')
print('=' * 60)
res_path = os.path.join(F.OUT_DIR, 'results_syk_se.json')
if not os.path.exists(res_path):
    sys.exit(f'ERROR: {res_path} not found -- run the syk_se job first.')
seed = json.load(open(res_path))
if not seed.get('ok'):
    sys.exit('ERROR: syk_se result is marked failed; nothing to seed from.')

F.DATA = F.load(True)
F.G['two_kinase'] = True
F.G['weighted'] = False
F.G['pnames'] = ['C1', 'C1S', 'C2', 'g', 'k3', 'k4']
lb = np.array([F.LOG_BOUNDS[n][0] for n in F.G['pnames']])
ub = np.array([F.LOG_BOUNDS[n][1] for n in F.G['pnames']])
x0 = np.clip(np.log10([seed['params'][n] for n in F.G['pnames']]), lb, ub)

print(f'  seed (from syk_se): raw SSR would be {F.raw_ssr(x0):.4e}; polishing...', flush=True)
r1 = minimize(F.cost, x0=x0, method='Nelder-Mead',
              options={'maxiter': 8000, 'xatol': 1e-9, 'fatol': 1e-9})
r2 = minimize(F.cost, x0=np.clip(r1.x, lb, ub), method='L-BFGS-B',
              bounds=list(zip(lb, ub)), options={'maxiter': 3000, 'ftol': 1e-13})
x = np.clip(r2.x if r2.fun < r1.fun else r1.x, lb, ub)
ssr = F.raw_ssr(x)

p = F._params(x)
dg = F.DATA['gamma']
mdl = F.solve(dg['t'], dg['pz_at_t'], dg['ps_at_t'], p, dg['y'][0], True)
i = int(np.argmax(dg['y']))
gap = (dg['y'][i] - mdl[i]) / max(dg['se'][i], 1e-9)

path, ssr_plot, _ = F.plot('syk_seeded', F.MODES['syk'], x, F.DATA)
out = dict(label='two-kinase, raw SSR (seeded from syk_se)', ok=True,
           raw_ssr=float(ssr), gamma_peak_gap_in_SE=round(float(gap), 2),
           params={n: float(p[n]) for n in F.G['pnames']}, plot=path)
with open(os.path.join(F.OUT_DIR, 'results_syk_seeded.json'), 'w') as fjs:
    json.dump(out, fjs, indent=2)

print(f'  DONE  raw SSR = {ssr:.4e}   gamma peak gap = {gap:+.1f} SE')
print(f'  params: ' + '  '.join(f'{n}={p[n]:.3g}' for n in F.G['pnames']))
print(f'  plot:   {path}')
print()
print('Reference points:  se (single-kinase shared) = 8.822e4, gap +7.4 SE')
print('                   syk_se (two-kinase)       = 5.362e4, gap +5.2 SE')
print('                   c02 (old shared best)     = 8.402e4')
