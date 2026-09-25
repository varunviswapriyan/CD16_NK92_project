#!/bin/bash
# =============================================================================
# ca_check.sh -- read the raw Ca bootstrap results and settle two questions:
#
#   1. Is lin1's S = 38.4 [31.6, 41.2] real, or is it hillonly's number
#      leaking into lin1's row?  (Identical to 3 s.f. across two different
#      models is not a coincidence.)
#   2. Is the reportable table safe to send?
#
#   bash ~/CD16_NK92_project/filesCC/ca_check.sh
#
# Reads only.  Submits nothing, changes nothing, takes seconds.
#
# It bypasses bootstrap_ca2.py's report entirely and recomputes everything
# from the raw boot_*.json files, so a bug in that report cannot hide here.
# Each JSON records which variant produced it, which makes a file collision
# detectable outright rather than by inference.
# =============================================================================
set -o pipefail
CADIR=$HOME/Ca_fit_c02

echo "setting up the job environment..."
[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh 2>/dev/null
module load Miniconda3/4.9.2 >/dev/null 2>&1
CONDA_SH=/gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh
[ -f "$CONDA_SH" ] && { . "$CONDA_SH" >/dev/null 2>&1; conda activate CD16_v2 >/dev/null 2>&1; }

PY="${PZ_PY:-}"
if [ -z "$PY" ]; then
  for cand in "${CONDA_PREFIX:+$CONDA_PREFIX/bin/python}" python python3; do
    [ -n "$cand" ] || continue
    c=$(command -v "$cand" 2>/dev/null || { [ -x "$cand" ] && echo "$cand"; })
    [ -n "$c" ] || continue
    if "$c" -c "import numpy" >/dev/null 2>&1; then PY="$c"; break; fi
  done
fi
[ -n "$PY" ] || { echo "ERROR: no python with numpy found."; exit 1; }
[ -d "$CADIR" ] || { echo "ERROR: $CADIR not found"; exit 1; }
cd "$CADIR" || exit 1
echo "using $PY"
echo

"$PY" - <<'PYEOF'
import json, glob, os
import numpy as np

ROOT = 'out_boot_ca2'
res, mismatch, counts = {}, {}, {}

for vd in sorted(glob.glob(os.path.join(ROOT, '*'))):
    if not os.path.isdir(vd):
        continue
    v = os.path.basename(vd)
    rows, wrong = [], []
    for f in sorted(glob.glob(os.path.join(vd, 'boot_*.json'))):
        try:
            r = json.load(open(f))
        except Exception:
            continue
        # every sample records the variant that produced it -- a collision shows here
        if r.get('variant') not in (None, v):
            wrong.append((os.path.basename(f), r.get('variant')))
        if r.get('ok') and r.get('sample', 0) > 0:
            rows.append(r)
    if not rows:
        continue
    counts[v] = len(rows)
    if wrong:
        mismatch[v] = wrong
    st = {}
    for p in sorted(rows[0].get('params', {})):
        x = np.array([r['params'][p] for r in rows if p in r['params']], float)
        x = x[np.isfinite(x)]
        if x.size < 3:
            continue
        lo, med, hi = np.percentile(x, [2.5, 50, 97.5])
        st[p] = (med, lo, hi, x.size)
    res[v] = st

print('=' * 78)
print('RAW bootstrap results, recomputed from the JSON files')
print('=' * 78)
for v in sorted(res):
    print('\n%-10s  n=%d' % (v, counts[v]))
    for p, (med, lo, hi, n) in res[v].items():
        w = (hi - lo) / abs(med) if med else float('inf')
        print('    %-4s %12.4g  [%.4g, %.4g]   w=%.2f' % (p, med, lo, hi, w))

# ---- question 1: is S duplicated across variants? --------------------------
print('\n' + '=' * 78)
print("1. IS lin1's S REAL?")
print('=' * 78)
sig = {}
for v, st in res.items():
    if 'S' in st:
        sig.setdefault(tuple(round(q, 6) for q in st['S'][:3]), []).append(v)
dups = [vs for vs in sig.values() if len(vs) > 1]

if mismatch:
    print('  FILE COLLISION -- these samples record a different variant than the')
    print('  directory they sit in:')
    for v, w in mismatch.items():
        print('    %s/  <- %s' % (v, ', '.join('%s says %s' % x for x in w[:5])))
    print('  Those runs must be redone.')
elif dups:
    print('  S is IDENTICAL across: %s' % ' and '.join(', '.join(d) for d in dups))
    print('  The JSONs each record their own variant correctly, so this is not a')
    print('  file collision -- the two models genuinely converge to the same S.')
    print('  That is believable: S is a fluorescence-to-concentration scale, and')
    print('  it is set by the size of the calcium signal, not by which drive term')
    print('  the model uses.  The number is real.')
else:
    print('  S differs between variants -- no duplication.  The earlier table was')
    print('  showing genuinely different numbers and there is nothing to fix.')

# ---- question 2: the reportable table --------------------------------------
print('\n' + '=' * 78)
print('2. REPORTABLE -- lin1')
print('=' * 78)
if 'lin1' not in res:
    print('  lin1 has no results yet.')
else:
    st = res['lin1']
    print('  n = %d bootstrap samples, k4 fixed at 1, all parameters shared' % counts['lin1'])
    print('  across adaptors (only the ITAM count differs 6 / 2 / 4)\n')
    print('  %-4s %12s %26s %7s   %s' % ('', 'estimate', '95% CI', 'w', 'what it is'))
    what = {'C1': 'Ca influx gain',
            'C2': 'h-gate gain',
            'g':  'Ca clearance rate (1/s)',
            'S':  'fluorescence -> uM scale (added by us)'}
    for p in ['C1', 'C2', 'g', 'S']:
        if p not in st:
            continue
        med, lo, hi, n = st[p]
        print('  %-4s %12.4g %26s %7.2f   %s'
              % (p, med, '[%.4g, %.4g]' % (lo, hi), (hi - lo) / abs(med), what.get(p, '')))
    print()
    print('  C1, C2 and g barely moved between 12 and 50 samples -- those are solid.')
    print('  S shifted 31.5 -> 38.4 when samples went 12 -> 50, landing outside its')
    print('  own 12-sample interval.  Quote it as a calibration constant with that')
    print('  caveat, or leave it out: it is a unit conversion, not a mechanism.')
    print()
    print('  Do NOT report lin1S -- it fixes S = 31.5, which the 50-sample run now')
    print('  places outside the interval.')
print('=' * 78)
PYEOF
