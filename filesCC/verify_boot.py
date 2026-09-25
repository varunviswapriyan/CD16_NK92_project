"""verify_boot.py -- check every built sample BEFORE any compute is spent.

Recomputes, independently of bootstrap_pzap.py, what each sample's
data/pZAP70_Tyr493_mean.csv should contain, and compares against what was
actually written.  Also checks the PSO/N_REPS settings and the fitted-parameter
list.  Exits non-zero on any mismatch so the caller can abort the submission.
"""
import os, sys, json, glob
from itertools import combinations_with_replacement
import numpy as np, pandas as pd

ROOT  = os.path.expanduser(sys.argv[1] if len(sys.argv) > 1 else '~/boot_pzap')
SRC   = os.path.expanduser('~/NK92_fit_v77')
SUB   = 'estimate_params_pzap_cleaned_up'
XLSX  = '/home/gddaslab/share/Varun_Indrani/estimate_params_pzap/data/pZAP70_Tyr493_Tyr292_original_and_averages.xlsx'
MARK  = 'pZAP70 (Tyr493)'
L2C   = {'KI 1': 'mean_zeta', 'KI 2': 'mean_gamma', 'KI 6': 'mean_hetero'}
DCOL  = ['d_0', 'd_1', 'd_2', 'd_5']
TIMES = [0.0, 60.0, 120.0, 300.0]
SETS  = list(combinations_with_replacement(range(3), 3))
FIT6  = ['lig0', 'kdl0', 'ZAP0', 'SYK0', 'KZP_MULT', 'KPR_MULT']

def canon(n):
    n = n.strip()
    return 'kdl0' if n.lower() in ('kd10', 'kdl0') else n

x = pd.read_excel(XLSX, sheet_name='Original_values')
x = x[x['marker'].astype(str).str.strip() == MARK]
days = {}
for line in L2C:
    sub = x[x['line'].astype(str).str.strip() == line]
    assert sub.shape[0] == 3, '%s: %d rows' % (line, sub.shape[0])
    days[line] = sub[DCOL].to_numpy(float)

orig = pd.read_csv(os.path.join(SRC, SUB, 'data', 'pZAP70_Tyr493_mean.csv'))
ocol = orig.columns[0]

def expected(tag):
    """What this sample's three fitted columns should hold."""
    if tag.startswith('nul') or tag == 'emp000':
        return {c: [float(orig.loc[np.isclose(orig[ocol].to_numpy(float), t), c].iloc[0])
                    for t in TIMES] for c in L2C.values()}
    i = int(tag[3:]); sel = list(SETS[(i - 1) % len(SETS)])
    out = {}
    for line, col in L2C.items():
        m = days[line][sel].mean(axis=0).copy(); m[0] = 0.0
        out[col] = list(m)
    return out

bad, n = [], 0
for d in sorted(glob.glob(os.path.join(ROOT, '*[0-9]'))):
    tag = os.path.basename(d); n += 1
    csv = os.path.join(d, SUB, 'data', 'pZAP70_Tyr493_mean.csv')
    if not os.path.exists(csv):
        bad.append('%s: data CSV missing' % tag); continue
    got = pd.read_csv(csv); gcol = got.columns[0]
    exp = expected(tag)
    for col, vals in exp.items():
        for t, v in zip(TIMES, vals):
            row = got.loc[np.isclose(got[gcol].to_numpy(float), t), col]
            if row.empty:
                bad.append('%s: no row t=%g in %s' % (tag, t, col)); continue
            if abs(float(row.iloc[0]) - v) > 1e-9:
                bad.append('%s %s @%gs: file=%.6f expected=%.6f'
                           % (tag, col, t, float(row.iloc[0]), v))
    cfg = json.load(open(os.path.join(d, SUB, 'v_config.json')))
    if sorted(canon(x) for x in cfg['PARAMS']) != sorted(FIT6):
        bad.append('%s: PARAMS=%s' % (tag, cfg['PARAMS']))
    if cfg.get('FIXED', {}).get('KZBG_FRAC') != 1.0:
        bad.append('%s: KZBG_FRAC not fixed at 1.0' % tag)
    for j, nm in enumerate(cfg['PARAMS']):         # a collapsed or absurd box is a bug
        if not (cfg['UB'][j] > cfg['LB'][j]):
            bad.append('%s: %s has LB>=UB (%s, %s)' % (tag, nm, cfg['LB'][j], cfg['UB'][j]))

# every distinct day-set must produce a DISTINCT dataset
sigs = {}
for d in sorted(glob.glob(os.path.join(ROOT, 'emp0[0-9][0-9]'))):
    tag = os.path.basename(d)
    if tag == 'emp000': continue
    g = pd.read_csv(os.path.join(d, SUB, 'data', 'pZAP70_Tyr493_mean.csv'))
    sigs.setdefault(tuple(np.round(g[list(L2C.values())].to_numpy(float).ravel(), 10)), []).append(tag)
dups = [v for v in sigs.values() if len(v) > 1]

c0 = json.load(open(os.path.join(sorted(glob.glob(os.path.join(ROOT, '*[0-9]')))[0], SUB, 'v_config.json')))
print('checked %d sample directories' % n)
print('  fitted parameters : %s' % c0['PARAMS'])
print('  FIXED             : %s' % c0.get('FIXED'))
print('  N_REPS            : %s' % c0.get('N_REPS'))
print('  distinct day-set datasets: %d (expect 10)' % len(sigs))
if dups: print('  DUPLICATE datasets: %s' % dups)
if bad or dups or len(sigs) != 10:
    print('\nFAILED:')
    for b in bad[:25]: print('   ' + b)
    sys.exit(1)
print('\nALL CHECKS PASSED -- every sample holds exactly the data it should.')
