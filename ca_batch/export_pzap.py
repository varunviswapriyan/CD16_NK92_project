#!/usr/bin/env python
"""
export_pzap.py  -  simulate the optimized pZAP model of THIS fit folder and write the curves
in the format Indrani's Ca code expects:  time,mean_zeta,mean_gamma,mean_hetero  (molecule numbers).

Run from inside NK92_fit_vNN/estimate_params_pzap_cleaned_up (needs CD16_v2 env):
    python export_pzap.py [n_reps] [out.csv]         default n_reps=10, out=model_output_pzap.csv
Reads the fitted parameters exactly like check_fit.py (any number of params + FIXED extras).
"""
import sys, os, json
import numpy as np
import pandas as pd
from calculate_SSR_gamma import gamma_PZAP
from calculate_SSR_zeta import zeta_PZAP
from calculate_SSR_mixed import mixed_PZAP

BASE4 = ['lig0', 'kd10', 'ZAP0', 'SYK0']
n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
out = sys.argv[2] if len(sys.argv) > 2 else 'model_output_pzap.csv'

FIXED = {}
if os.path.exists('v_config.json'):
    FIXED = json.load(open('v_config.json')).get('FIXED', {})

with open('analysis_param_residue.dat') as f:
    lines = [ln.strip() for ln in f if ln.strip()]
names, est = None, None
for k, ln in enumerate(lines):
    if ln.startswith('lig0\tkd10\tZAP0\tSYK0'):
        names = ln.split('\t')[:-1]
        est = np.array([float(v) for v in lines[k + 1].split('\t')[:len(names)]])
        break
if est is None:
    names = BASE4
    est = np.fromstring(lines[0].split('[', 1)[1].rsplit(']', 1)[0], sep=' ')[:4]
vals = {nm: float(10 ** v) for nm, v in zip(names, est)}
lig0, kd10, ZAP0, SYK0 = (vals[k] for k in BASE4)
print('params:', '  '.join(f'{k}={v:.4g}' for k, v in vals.items()), ' FIXED:', FIXED)


def extras_for(kind):
    o = {}
    for src in (FIXED, vals):
        for k, v in src.items():
            if k in BASE4: continue
            if ':' in k:
                mk, nm = k.split(':', 1)
                if mk == kind: o[nm] = v
            else:
                o[k] = v
    return o


kw = {}
try:
    import inspect
    if 'extra' in inspect.signature(gamma_PZAP).parameters:
        kw = {'gamma': dict(extra=extras_for('gamma')), 'zeta': dict(extra=extras_for('zeta')), 'mixed': dict(extra=extras_for('mixed'))}
except Exception:
    pass
tg, pg = gamma_PZAP(n, lig0, kd10, ZAP0, SYK0, 'export', **kw.get('gamma', {}))
tz, pz = zeta_PZAP(n, lig0, kd10, ZAP0, SYK0, 'export', **kw.get('zeta', {}))
tm, pm = mixed_PZAP(n, lig0, kd10, ZAP0, SYK0, 'export', **kw.get('mixed', {}))

grid = np.linspace(0.0, 300.0, 1001)
def on_grid(t, y):
    t = np.asarray(t, float); y = np.asarray(y, float); o = np.argsort(t)
    return np.interp(grid, t[o], y[o])
df = pd.DataFrame({'time': grid, 'mean_zeta': on_grid(tz, pz), 'mean_gamma': on_grid(tg, pg), 'mean_hetero': on_grid(tm, pm)})
df.to_csv(out, index=False)
print(f'wrote {out}: {len(df)} rows, {n} replicates')
print(f'peak molecules  zeta={df.mean_zeta.max():.0f}  gamma={df.mean_gamma.max():.0f}  hetero={df.mean_hetero.max():.0f}')
print(f'peak uM (/{25*602})  zeta={df.mean_zeta.max()/15050:.4f}  gamma={df.mean_gamma.max()/15050:.4f}  hetero={df.mean_hetero.max()/15050:.4f}')
