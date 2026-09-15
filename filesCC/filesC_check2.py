"""
filesC_check2.py -- is the 2% gamma-vs-hetero pSYK difference REAL or noise?

Two tests, one command (~1 min, login node fine):
  1. Per-trace statistics: recomputes sustained pSYK (t>100s) for EACH of the
     ~30 individual NFsim traces per condition, prints mean +/- SEM, and the
     z-score of the gamma-hetero difference. |z| > 2 => the difference is a
     real model prediction; |z| < 2 => within stochastic simulation noise.
  2. Swap test: evaluates the seeded two-kinase fit with gamma/hetero pSYK
     inputs EXCHANGED. A large SSR jump means the fit genuinely rides on that
     small pSYK difference (mechanistically meaningful IF test 1 passes;
     fragile if it doesn't).

Usage:
    cd ~/Ca_fit_c02
    python ~/CD16_NK92_project/filesCC/filesC_check2.py
"""

import os, sys, glob, json
import numpy as np

sys.path.insert(0, os.path.expanduser('~/CD16_NK92_project/filesCC'))
import filesC as F

# ---------- 1. per-trace sustained pSYK statistics ----------
print('=' * 64)
print('Per-trace pSYK statistics (sustained = mean over t > 100 s)')
print('=' * 64)
vals = {}
for cond, sub in F.RUNS_MAP.items():
    adirs = sorted(glob.glob(os.path.join(F.V69_DIR, sub, 'analysis*')),
                   key=os.path.getmtime)[-F.K_RECENT:]
    per_trace = []
    for dd in adirs:
        for g in glob.glob(os.path.join(dd, '**', '*.gdat'), recursive=True):
            t, v = F._read_gdat_psyk(g)
            if t is None: continue
            per_trace.append(float(np.mean(v[t > 100])))
    a = np.array(per_trace)
    vals[cond] = a
    print(f'  {cond:<7} n={len(a):<3} mean={a.mean():9.2f}  SD={a.std(ddof=1):7.2f}'
          f'  SEM={a.std(ddof=1)/np.sqrt(len(a)):6.2f}')

gd, hd = vals['gamma'], vals['hetero']
diff = gd.mean() - hd.mean()
sed = np.sqrt(gd.var(ddof=1)/len(gd) + hd.var(ddof=1)/len(hd))
z = diff / sed if sed > 0 else float('inf')
print(f'\n  gamma - hetero (sustained): {diff:+.2f}  (z = {z:+.2f})')
if abs(z) > 2:
    print('  -> SIGNIFICANT: the pSYK difference is a real model prediction.')
else:
    print('  -> NOT significant at ~2 SEM: within stochastic simulation noise.')
    print('     More NFsim replicates would be needed to establish it.')

# ---------- 2. swap test on the seeded fit ----------
print()
print('=' * 64)
print('Swap test: seeded two-kinase fit with gamma/hetero pSYK exchanged')
print('=' * 64)
res_path = os.path.join(F.OUT_DIR, 'results_syk_seeded.json')
if not os.path.exists(res_path):
    res_path = os.path.join(F.OUT_DIR, 'results_syk_se.json')
seed = json.load(open(res_path))
F.DATA = F.load(True)
F.G['two_kinase'] = True; F.G['weighted'] = False
F.G['pnames'] = ['C1', 'C1S', 'C2', 'g', 'k3', 'k4']
x = np.log10([seed['params'][n] for n in F.G['pnames']])

ssr_orig = F.raw_ssr(x)
g_ps = F.DATA['gamma']['ps_at_t'].copy()
h_ps = F.DATA['hetero']['ps_at_t'].copy()
F.DATA['gamma']['ps_at_t'], F.DATA['hetero']['ps_at_t'] = h_ps, g_ps
ssr_swap = F.raw_ssr(x)
F.DATA['gamma']['ps_at_t'], F.DATA['hetero']['ps_at_t'] = g_ps, h_ps

print(f'  raw SSR, original inputs: {ssr_orig:.4e}')
print(f'  raw SSR, swapped inputs:  {ssr_swap:.4e}')
ratio = ssr_swap / ssr_orig
print(f'  ratio: {ratio:.2f}x')
if ratio > 1.2:
    print('  -> Fit depends strongly on the gamma/hetero pSYK difference.')
    print('     Meaningful mechanism IF test 1 was significant; fragile if not.')
else:
    print('  -> Fit barely changes when swapped: the SSR gain comes from the')
    print('     overall pSYK dynamics, not the gamma/hetero ordering. The')
    print('     two-kinase improvement is robust, but it does not by itself')
    print('     explain the gamma/hetero separation.')
