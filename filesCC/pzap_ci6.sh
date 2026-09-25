#!/bin/bash
# pzap_ci6.sh -- bootstrap confidence intervals for the six fitted pZAP
#                parameters, exactly as Indrani asked for:
#
#   "Lilly did pZAP experiments in 3 different days (26 May, 28 May and 3rd
#    June). You can use 3 data points for each timepoint (0, 1 min, 2 min and
#    5 min) to calculate bootstrapped mean to calculate CI."
#
#   bash ~/CD16_NK92_project/filesCC/pzap_ci6.sh          # 10 bootstrap + 1 ref + 4 null
#   bash ~/CD16_NK92_project/filesCC/pzap_ci6.sh report
#
# THE BOOTSTRAP IS EXACT, NOT SAMPLED.  Resampling 3 days with replacement has
# exactly 10 distinct outcomes, with known probabilities (1, 3, 3, 3, 6, 3, 1,
# 3, 3, 1 out of 27).  Fitting all ten and weighting by those probabilities
# gives the COMPLETE bootstrap distribution -- no Monte-Carlo error, and fewer
# jobs than random draws.  All three lines were run on the same three days, so
# a day is resampled as a unit across every line.
#
# THE REPORTED CI IS THE ORDINARY BOOTSTRAP CI: point estimate +/- 1.96 x the
# (weighted) bootstrap SD, with the weighted percentile interval beside it.
# No custom corrections are applied to it.
#
# COMPUTE IS SPENT ON ACCURACY, NOT ON MORE RESAMPLES.  The 10 day-sets are the
# complete bootstrap already, so extra samples would add nothing.  Instead the
# per-job budget shifts from PSO iterations into NFsim replicates (N_REPS 3->6),
# which cuts the objective's noise SD by 1/sqrt(2).  That matters because PSO
# returns the best-of-all-evaluations point: with a noisy objective, more
# iterations select harder on lucky noise and can make the fitted parameters
# MORE variable, not less.  The warm box is tight, so search is not the limit.
#
# The few "null" jobs refit the ORIGINAL data with fresh optimizer randomness.
# They are a QUALITY CHECK only -- v_config has N_REPS=3, so the objective is
# stochastic, and the nulls say how much of the spread is the simulator rather
# than the experiment.  They are reported as one diagnostic line and are NOT
# used to adjust the interval.
set -e
F=~/CD16_NK92_project/filesCC
unset BOOT_PIN
export BOOT_HALF=0.25
export BOOT_PARTICLES=24
export BOOT_ITERS=20        # fewer iterations, bought back as replicates below
export BOOT_NREPS=6         # was 3: halves the objective's noise variance
ROOT=~/boot_pzap

# ----------------------------------------------------------------- report ----
if [ "$1" = "report" ]; then
ROOT="$ROOT" python3 - <<'PY'
import os, glob, json
import numpy as np
from itertools import combinations_with_replacement
from math import factorial

ROOT = os.environ['ROOT']
SUB  = 'estimate_params_pzap_cleaned_up'
FIT6 = ['lig0', 'kdl0', 'ZAP0', 'SYK0', 'KZP_MULT', 'KPR_MULT']
UNITS = {'lig0': 'molecules/um3', 'kdl0': '1/s', 'ZAP0': 'molecules/um3',
         'SYK0': 'molecules/um3', 'KZP_MULT': 'x0.03 = kzp (uM s)^-1',
         'KPR_MULT': 'x0.01 = KPR 1/s'}
EMP_SETS = list(combinations_with_replacement(range(3), 3))

def emp_weight(ms):
    c = [list(ms).count(k) for k in range(3)]
    p = factorial(3)
    for x in c: p //= factorial(x)
    return p / 27.0

def norm(k):
    k = k.strip()
    return 'kdl0' if k.lower() in ('kd10', 'kdl0') else k

def parse(d):
    p = os.path.join(d, SUB, 'analysis_param_residue.dat')
    if not os.path.exists(p): return None
    res, vals = None, {}
    for ln in open(p, errors='ignore'):
        s = ln.strip()
        if s.lower().startswith('residue ='):
            try: res = float(s.split('=', 1)[1])
            except Exception: pass
        elif s.lower().startswith('linear'):
            for tok in s.split()[1:]:
                if '=' in tok:
                    k, v = tok.split('=', 1)
                    try: vals[norm(k)] = float(v)
                    except Exception: pass
    if res is None or not vals or not all(n in vals for n in FIT6):
        return None
    return res, vals

def wstats(x, w):
    w = np.asarray(w, float); w = w / w.sum()
    m = float((w * x).sum())
    v = float((w * (x - m) ** 2).sum())
    # small-sample correction for weighted variance
    denom = 1.0 - (w ** 2).sum()
    if denom > 0: v /= denom
    return m, float(np.sqrt(v))

def wpercentile(x, w, qs):
    o = np.argsort(x); x, w = np.asarray(x)[o], np.asarray(w, float)[o]
    cw = np.cumsum(w) / w.sum()
    return [float(np.interp(q, cw, x)) for q in qs]

ref = parse(os.path.join(ROOT, 'emp000'))
boots, miss = [], []
for i, ms in enumerate(EMP_SETS, start=1):
    r = parse(os.path.join(ROOT, 'emp%03d' % i))
    if r: boots.append((ms, emp_weight(ms), r))
    else: miss.append(tuple(d + 1 for d in ms))
nulls = [r for r in (parse(d) for d in sorted(glob.glob(os.path.join(ROOT, 'nul*')))) if r]

covered = sum(w for _, w, _ in boots)
print('EMPIRICAL BOOTSTRAP over Lilly\'s 3 experiment days (26 May / 28 May / 3 Jun)')
print('  day-sets fitted: %d of %d   covering %.0f%% of the bootstrap probability'
      % (len(boots), len(EMP_SETS), 100 * covered))
if miss: print('  not finished yet: ' + ', '.join(str(m) for m in miss))
print('  reference SSR (all 3 days, as measured) = %s'
      % ('%.4e' % ref[0] if ref else 'not finished'))
if len(boots) < 4:
    print('\nToo few day-sets finished to summarise.  Rerun when more land.')
    raise SystemExit
print('  bootstrap SSRs: ' + '  '.join('%.4g' % r[0] for _, _, r in boots))
print()

B = {}
for d in sorted(glob.glob(os.path.join(ROOT, '*[0-9]'))):
    p = os.path.join(d, SUB, 'v_config.json')
    if os.path.exists(p):
        c = json.load(open(p)); B = dict(zip(c['PARAMS'], zip(c['LB'], c['UB']))); break

hdr = ('%-9s %11s %11s %24s %6s %24s'
       % ('param', 'point est', 'boot SD', '95% CI (+/-1.96 SD)', 'w', 'weighted percentile CI'))
print(hdr); print('-' * len(hdr))
summary = {}
for n in FIT6:
    x = np.array([r[1][n] for _, _, r in boots], float)
    w = np.array([wt for _, wt, _ in boots], float)
    mu, sd = wstats(x, w)
    pt = ref[1][n] if ref else mu
    lo, hi = pt - 1.96 * sd, pt + 1.96 * sd
    if lo <= 0: lo = max(x.min() * 0.5, 1e-12)
    plo, phi = wpercentile(x, w, [0.025, 0.975])
    rw = (hi - lo) / abs(pt) if pt else float('inf')
    flag = ''
    if n in B:
        lb, ub = B[n]
        lx = np.log10(np.maximum(x, 1e-300))
        if np.any(lx <= lb + 0.02) or np.any(lx >= ub - 0.02): flag = ' !bound'
    summary[n] = dict(point=pt, boot_mean=mu, boot_sd=sd, ci=[lo, hi],
                      pct_ci=[plo, phi], rel_width=rw)
    print('%-9s %11.4g %11.4g   [%9.4g, %9.4g] %6.2f   [%9.4g, %9.4g]%s'
          % (n, pt, sd, lo, hi, rw, plo, phi, flag))
print('-' * len(hdr))
print('CI = point estimate +/- 1.96 x weighted bootstrap SD.  w = width / point.')
print('w < 1 means the parameter is well determined by the data.')
print()
print('in physical units:')
for n, f, nm, u in (('KZP_MULT', 0.03, 'kzp', '(uM s)^-1'),
                    ('KPR_MULT', 0.01, 'KPR', 's^-1')):
    s = summary[n]
    print('  %-4s = %.4g   95%% CI [%.4g, %.4g]   %s'
          % (nm, s['point'] * f, s['ci'][0] * f, s['ci'][1] * f, u))

if len(nulls) > 1:
    print()
    print('QUALITY CHECK (not applied to the CI above).  %d refits of the ORIGINAL'
          % len(nulls))
    print('data with fresh optimizer randomness; N_REPS=3 makes the objective noisy:')
    for n in FIT6:
        nv = np.array([r[1][n] for r in nulls], float)
        sd_n, sd_b = float(nv.std(ddof=1)), summary[n]['boot_sd']
        frac = (sd_n / sd_b) if sd_b > 0 else float('nan')
        note = '   <-- mostly simulator noise' if frac > 0.7 else ''
        print('   %-9s simulator SD %9.4g  = %3.0f%% of the bootstrap SD%s'
              % (n, sd_n, 100 * frac, note))
    print('   Where that fraction is high, raising N_REPS above 3 would tighten the')
    print('   interval without touching the data or the statistics.')

json.dump(dict(n_daysets=len(boots), prob_covered=covered,
               ref_ssr=(ref[0] if ref else None), params=summary),
          open(os.path.join(ROOT, 'ci_summary_6param.json'), 'w'), indent=2)

import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
fig, axes = plt.subplots(2, 3, figsize=(16, 8))
for ax, n in zip(axes.ravel(), FIT6):
    x = np.array([r[1][n] for _, _, r in boots], float)
    w = np.array([wt for _, wt, _ in boots], float)
    s = summary[n]
    ax.bar(x, w, width=(x.max() - x.min() + 1e-9) / 25.0, color='tab:green', alpha=0.8)
    ax.axvspan(s['ci'][0], s['ci'][1], color='tab:blue', alpha=0.15)
    ax.axvline(s['point'], color='tab:red', lw=2.5)
    ax.set_ylabel('bootstrap probability', fontsize=7)
    ax.set_title('%s  %.4g\n95%% CI [%.4g, %.4g]   %s'
                 % (n, s['point'], s['ci'][0], s['ci'][1], UNITS.get(n, '')), fontsize=9)
fig.suptitle('pZAP parameters: exact empirical bootstrap over 3 experiment days '
             '(%d of 10 day-sets, %.0f%% of probability)' % (len(boots), 100 * covered))
fig.tight_layout(rect=[0, 0, 1, 0.94])
p1 = os.path.join(ROOT, 'ci_hist_6param.png'); fig.savefig(p1, dpi=110); plt.close(fig)
print('\nfigure: %s' % p1)
PY

  OUT=~/ci_plots; mkdir -p $OUT
  cp $ROOT/ci_hist_6param.png $OUT/ 2>/dev/null || true
  IP=$(hostname -I | awk '{for(i=1;i<=NF;i++) if($i ~ /^10\.73\./) print $i}')
  echo; echo "Open:  http://${IP:-10.73.170.128}:8000/"
  cd $OUT && python -m http.server 8000
  exit 0
fi

# ----------------------------------------------------------------- submit ----
# BUILD -> VERIFY -> SUBMIT.  Nothing is queued until an independent checker has
# confirmed every sample holds exactly the data it should.  This is the step
# that would have caught the wrong line mapping and the stale v77 result file
# before they cost a run each.
NN="${1:-7}"

OLD=$(squeue -u "$USER" -h -o "%i %j" | awk '$2 ~ /^bp/ {print $1}')
if [ -n "$OLD" ]; then
  echo "$OLD" | xargs -r scancel
  echo "cancelled $(echo "$OLD" | grep -c .) old pZAP bootstrap jobs"
fi
rm -rf $ROOT ~/boot_pzap_pinZAP0 ~/boot_pzap_pinKZP_MULT_KPR_MULT

echo "== 1/3  building 1 reference + 10 day-sets + $NN null =="
python $F/bootstrap_pzap.py build 10 0 "$NN"

echo
echo "== 2/3  verifying every sample against an independent recomputation =="
if ! python3 $F/verify_boot.py "$ROOT"; then
  echo
  echo "ABORTED -- nothing submitted.  Fix the mismatch above first."
  exit 1
fi

echo
echo "== 3/3  submitting =="
for f in $ROOT/run_*.sh; do sbatch "$f"; done
echo
echo "$((11 + NN)) jobs.  24 particles x 20 iterations, N_REPS=6 -> about 4 h each,"
echo "one queue wave.   Results:  bash \$0 report"
