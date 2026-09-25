#!/bin/bash
# pzap_ci6.sh -- confidence intervals for the SIX fitted pZAP parameters
#                (lig0, kdl0, ZAP0, SYK0, kzp, KPR).
#
#   bash ~/CD16_NK92_project/filesCC/pzap_ci6.sh           # submit 24 + 1 reference
#   bash ~/CD16_NK92_project/filesCC/pzap_ci6.sh 16        # smaller if the queue is busy
#   bash ~/CD16_NK92_project/filesCC/pzap_ci6.sh report    # safe to run on partial results
#
# DESIGN NOTES
#  * Parametric only.  3 days admit just 10 distinct empirical resamples, so
#    that arm's percentiles are the min/max of a tiny discrete set.  Sample 0
#    (emp000) is the reference fit and leaves the shipped CSV untouched.
#  * Warm box 0.45 log10 -- each job is 16 particles x 20 iterations = 320
#    evaluations for 6 dimensions, so the box has to be small enough to converge
#    in that budget.  The report flags any replicate that lands on a bound.
#  * Replicates whose SSR is far worse than the reference did not find their own
#    optimum; they are trimmed, and the count is printed.
#  * CI = point estimate +/- 1.96 x bootstrap SD.  At n ~ 20 the 2.5th
#    percentile is essentially the sample minimum; the SD is stable.  Percentiles
#    are printed alongside so nothing is hidden.
#  * A result counts only if it reports ALL SIX parameters.  A directory that is
#    still running holds no result at all now, but an older tree may still carry
#    v77's 4-parameter file -- that is what the completeness check rejects.
set -e
F=~/CD16_NK92_project/filesCC
unset BOOT_PIN
export BOOT_HALF=0.45
ROOT=~/boot_pzap
TRIM=${TRIM:-2.0}

# ----------------------------------------------------------------- report ----
if [ "$1" = "report" ]; then
ROOT="$ROOT" TRIM="$TRIM" python3 - <<'PY'
import os, glob, json
import numpy as np

ROOT = os.environ['ROOT']; TRIM = float(os.environ['TRIM'])
SUB  = 'estimate_params_pzap_cleaned_up'
FIT6 = ['lig0', 'kdl0', 'ZAP0', 'SYK0', 'KZP_MULT', 'KPR_MULT']
V77  = {'lig0': 85.7, 'kdl0': 0.00272, 'ZAP0': 146.8, 'SYK0': 52.9,
        'KZP_MULT': 2.2, 'KPR_MULT': 0.54}
UNITS = {'lig0': 'molecules/um3', 'kdl0': '1/s', 'ZAP0': 'molecules/um3',
         'SYK0': 'molecules/um3', 'KZP_MULT': 'x0.03 = kzp (uM s)^-1',
         'KPR_MULT': 'x0.01 = KPR 1/s'}

def norm(k):
    k = k.strip()
    return 'kdl0' if k.lower() in ('kd10', 'kdl0') else k

def parse(d):
    p = os.path.join(d, SUB, 'analysis_param_residue.dat')
    if not os.path.exists(p):
        return None
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
    if res is None or not vals:
        return None
    return res, vals

def complete(v):
    return all(n in v for n in FIT6)

def box():
    for d in sorted(glob.glob(os.path.join(ROOT, '*[0-9]'))):
        p = os.path.join(d, SUB, 'v_config.json')
        if os.path.exists(p):
            c = json.load(open(p))
            return dict(zip(c['PARAMS'], zip(c['LB'], c['UB'])))
    return {}

# ---- collect, rejecting anything that is not a full 6-parameter result ----
raw_ref = parse(os.path.join(ROOT, 'emp000'))
ref = raw_ref if (raw_ref and complete(raw_ref[1])) else None
stale_ref = bool(raw_ref) and ref is None

boots, stale, pending = [], 0, 0
for d in sorted(glob.glob(os.path.join(ROOT, 'par*'))):
    r = parse(d)
    if r is None:
        pending += 1
    elif complete(r[1]):
        boots.append(r)
    else:
        stale += 1

ndirs = len(glob.glob(os.path.join(ROOT, 'par*')))
print('parametric directories: %d   real 6-param results: %d   still running: %d   '
      '4-param leftovers ignored: %d' % (ndirs, len(boots), pending, stale))
if stale_ref:
    print('reference: directory still holds a 4-parameter file (v77 template), not a result.')
elif ref is None:
    print('reference: not finished yet.')
else:
    print('reference (un-jittered, all 6 fitted) SSR = %.4e' % ref[0])
    print('   ' + '  '.join('%s=%.4g' % (n, ref[1][n]) for n in FIT6))
print()

if len(boots) < 3:
    print('Fewer than 3 real replicates so far -- nothing to summarise yet.')
    print('Rerun this once more jobs land:  bash $0 report')
    raise SystemExit

allres = np.array(sorted(b[0] for b in boots), float)
print('replicate SSRs (sorted): ' + '  '.join('%.4g' % x for x in allres))
base = ref[0] if ref else float(np.median(allres))
kept = [b for b in boots if b[0] <= TRIM * base]
drop = len(boots) - len(kept)
if len(kept) < 3:
    kept, drop = boots, 0
    print('NOTE: trim would leave <3 replicates, so nothing was trimmed.')
res = np.array([b[0] for b in kept], float)
print('kept %d of %d (SSR <= %.1f x %s);  kept SSR median %.4e  range [%.4e, %.4e]'
      % (len(kept), len(boots), TRIM, 'reference' if ref else 'median', 
         np.median(res), res.min(), res.max()))
print()

B = box()
hdr = ('%-9s %10s %10s %10s %22s %6s %22s'
       % ('param', 'point', 'boot mean', 'boot SD', '95% CI (point +/- 1.96SD)',
          'w', 'percentile CI'))
print(hdr); print('-' * len(hdr))

summary, allv = {}, {}
for n in FIT6:
    a = np.array([b[1][n] for b in kept], float)
    allv[n] = a
    mu = float(a.mean())
    sd = float(a.std(ddof=1)) if a.size > 1 else 0.0
    pt = ref[1][n] if ref else mu
    lo, hi = pt - 1.96 * sd, pt + 1.96 * sd
    if lo <= 0: lo = max(a.min() * 0.5, 1e-12)
    plo, phi = np.percentile(a, [2.5, 97.5])
    w = (hi - lo) / abs(pt) if pt else float('inf')
    flag = ''
    if n in B:
        lb, ub = B[n]
        la = np.log10(np.maximum(a, 1e-300))
        if np.any(la <= lb + 0.02) or np.any(la >= ub - 0.02):
            flag = ' !bound'
    summary[n] = dict(point=pt, mean=mu, sd=sd, ci=[lo, hi],
                      pct=[float(plo), float(phi)], rel_width=w, n=int(a.size))
    print('%-9s %10.4g %10.4g %10.4g   [%8.4g, %8.4g] %6.2f   [%8.4g, %8.4g]%s'
          % (n, pt, mu, sd, lo, hi, w, plo, phi, flag))
print('-' * len(hdr))
print('w = CI width / point estimate.   w < 1 = well determined.')
print('!bound = a replicate sat on the edge of the search box, so that interval')
print('         is set by the box rather than by the data.')
print()
print('in physical units:')
for n, f, nm, u in (('KZP_MULT', 0.03, 'kzp', '(uM s)^-1'),
                    ('KPR_MULT', 0.01, 'KPR', 's^-1')):
    if n in summary:
        s = summary[n]
        print('  %-4s = %.4g  [%.4g, %.4g]  %s'
              % (nm, s['point'] * f, s['ci'][0] * f, s['ci'][1] * f, u))
print()
print('v77 (4-parameter fit) shown for orientation only -- different model:')
print('  ' + '  '.join('%s=%g' % (k, v) for k, v in V77.items()))

json.dump(dict(n_results=len(boots), n_kept=len(kept), n_dropped=drop,
               ref_ssr=(ref[0] if ref else None), params=summary),
          open(os.path.join(ROOT, 'ci_summary_6param.json'), 'w'), indent=2)

# ---- figures ----
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 3, figsize=(16, 8))
for ax, n in zip(axes.ravel(), FIT6):
    a, s = allv[n], summary[n]
    ax.hist(a, bins=max(6, len(a) // 3), color='tab:green', alpha=0.75)
    ax.axvspan(s['ci'][0], s['ci'][1], color='tab:blue', alpha=0.12)
    ax.axvline(s['ci'][0], color='k', ls='--', lw=1)
    ax.axvline(s['ci'][1], color='k', ls='--', lw=1)
    ax.axvline(s['point'], color='tab:red', lw=2.5)
    ax.set_title('%s   %.4g [%.4g, %.4g]\n%s'
                 % (n, s['point'], s['ci'][0], s['ci'][1], UNITS.get(n, '')), fontsize=9)
fig.suptitle('pZAP 6-parameter bootstrap, n=%d   red = point estimate, '
             'band = point +/- 1.96 SD' % len(kept))
fig.tight_layout(rect=[0, 0, 1, 0.94])
p1 = os.path.join(ROOT, 'ci_hist_6param.png'); fig.savefig(p1, dpi=110); plt.close(fig)

pairs = [('lig0', 'ZAP0'), ('ZAP0', 'SYK0'), ('ZAP0', 'KZP_MULT'),
         ('lig0', 'KZP_MULT'), ('kdl0', 'KPR_MULT'), ('SYK0', 'KZP_MULT')]
fig, axes = plt.subplots(1, len(pairs), figsize=(3.1 * len(pairs), 3.3))
print('\npairwise correlations (|r| near 1 = constrained only in combination):')
for ax, (x, y) in zip(np.atleast_1d(axes), pairs):
    ax.scatter(allv[x], allv[y], s=30, alpha=0.75, color='tab:blue')
    if ref:
        ax.scatter([ref[1][x]], [ref[1][y]], marker='*', s=220, color='tab:red', zorder=5)
    sx, sy = allv[x].std(), allv[y].std()
    r = float(np.corrcoef(allv[x], allv[y])[0, 1]) if sx > 0 and sy > 0 else float('nan')
    print('  %-10s vs %-10s  r = %+.2f' % (x, y, r))
    ax.set_xlabel(x, fontsize=8); ax.set_ylabel(y, fontsize=8)
    ax.set_title('r = %+.2f' % r, fontsize=9)
fig.suptitle('pairwise: tilted cloud = only the combination is constrained; '
             'red star = point estimate', fontsize=10)
fig.tight_layout(rect=[0, 0, 1, 0.90])
p2 = os.path.join(ROOT, 'ci_pairs_6param.png'); fig.savefig(p2, dpi=110); plt.close(fig)
print('\nfigures: %s\n          %s' % (p1, p2))
PY

  OUT=~/ci_plots; mkdir -p $OUT
  cp $ROOT/ci_hist_6param.png  $OUT/ 2>/dev/null || true
  cp $ROOT/ci_pairs_6param.png $OUT/ 2>/dev/null || true
  IP=$(hostname -I | awk '{for(i=1;i<=NF;i++) if($i ~ /^10\.73\./) print $i}')
  echo; ls -1 $OUT; echo
  echo "Open:  http://${IP:-10.73.170.128}:8000/"
  cd $OUT && python -m http.server 8000
  exit 0
fi

# ----------------------------------------------------------------- submit ----
NP="${1:-24}"
OLD=$(squeue -u "$USER" -h -o "%i %j" | awk '$2 ~ /^bp/ {print $1}')
if [ -n "$OLD" ]; then
  echo "$OLD" | xargs -r scancel
  echo "cancelled $(echo "$OLD" | grep -c .) old pZAP bootstrap jobs"
fi
rm -rf $ROOT ~/boot_pzap_pinZAP0 ~/boot_pzap_pinKZP_MULT_KPR_MULT
python $F/bootstrap_pzap.py submit 0 "$NP"
echo
echo "$((NP + 1)) jobs: 1 reference + $NP parametric replicates, all 6 parameters fitted."
echo "Roughly 1.5-2.5 h per job.   Partial results:  bash \$0 report"
