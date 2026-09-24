#!/bin/bash
# pzap_ci6.sh -- confidence intervals for the SIX fitted pZAP parameters
#                (lig0, kdl0, ZAP0, SYK0, kzp, KPR), inside a 5-hour window.
#
#   bash ~/CD16_NK92_project/filesCC/pzap_ci6.sh           # submit 24 + 1 reference
#   bash ~/CD16_NK92_project/filesCC/pzap_ci6.sh 16        # smaller if the queue is busy
#   bash ~/CD16_NK92_project/filesCC/pzap_ci6.sh report
#
# FOUR DELIBERATE CHOICES, each aimed at the CIs actually meaning something:
#
# 1. PARAMETRIC ONLY.  The empirical arm resamples 3 days, and 3 days admit just
#    10 distinct resamples in the whole universe -- its percentiles are the
#    min/max of a tiny discrete set, which is where the bimodal ZAP0 histogram
#    came from.  Every job goes to the parametric arm, which can draw unlimited
#    distinct datasets from the per-timepoint SD.  Sample 0 (emp000) is the
#    reference fit on the real un-resampled means.
#
# 2. WARM BOX 0.45 log10, not 0.60.  Each job is 16 particles x 20 iterations
#    = 320 cost evaluations for 6 dimensions -- sparse.  A narrower box lets PSO
#    actually converge in that budget, so the spread across replicates reflects
#    the DATA rather than where the optimizer happened to stop.  0.45 still puts
#    every wall well outside last run's observed spread; the report checks and
#    flags any replicate that lands on a bound anyway.
#
# 3. NON-CONVERGED REPLICATES ARE TRIMMED.  A replicate whose SSR is far worse
#    than the reference fit did not find its own optimum; including it measures
#    optimizer failure, not uncertainty.  They are dropped, and the count is
#    printed so the trimming is visible rather than silent.
#
# 4. CI = point estimate +/- 1.96 x bootstrap SD, not raw percentiles.  At
#    n ~ 20 the 2.5th percentile is essentially the sample minimum (one unlucky
#    replicate sets it), while the SD is estimated to about +/-16%.  The SD form
#    is stable at this n, is centred on the point estimate by construction, and
#    is the error-bar presentation people expect.  Percentiles are printed
#    beside it so nothing is hidden.
set -e
F=~/CD16_NK92_project/filesCC
unset BOOT_PIN                 # all six fitted
export BOOT_HALF=0.45
ROOT=~/boot_pzap
TRIM=${TRIM:-2.0}              # drop replicates with SSR > TRIM x reference SSR

# ----------------------------------------------------------------- report ----
if [ "$1" = "report" ]; then
ROOT="$ROOT" TRIM="$TRIM" python3 - <<'PY'
import os, glob, json
import numpy as np

ROOT = os.environ['ROOT']; TRIM = float(os.environ['TRIM'])
SUB  = 'estimate_params_pzap_cleaned_up'
FIT6 = ['lig0', 'kdl0', 'ZAP0', 'SYK0', 'KZP_MULT', 'KPR_MULT']
UNITS = {'lig0': 'molecules/um3', 'kdl0': '1/s', 'ZAP0': 'molecules/um3',
         'SYK0': 'molecules/um3', 'KZP_MULT': 'x0.03 = kzp (uM s)^-1',
         'KPR_MULT': 'x0.01 = KPR 1/s'}
V77 = {'lig0': 85.7, 'kdl0': 0.00272, 'ZAP0': 146.8, 'SYK0': 52.9,
       'KZP_MULT': 2.2, 'KPR_MULT': 0.54}

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
    if not vals or res is None:
        return None
    return res, vals

def box():
    """log10 search box actually used (same for every sample)."""
    for d in sorted(glob.glob(os.path.join(ROOT, '*[0-9]'))):
        p = os.path.join(d, SUB, 'v_config.json')
        if os.path.exists(p):
            c = json.load(open(p))
            return dict(zip(c['PARAMS'], zip(c['LB'], c['UB'])))
    return {}

ref = parse(os.path.join(ROOT, 'emp000'))
boots = []
for d in sorted(glob.glob(os.path.join(ROOT, 'par*'))):
    r = parse(d)
    if r: boots.append(r)

if ref is None:
    print('reference fit (emp000) not finished yet -- rerun when it lands.')
if not boots:
    print('no bootstrap replicates finished yet.'); raise SystemExit

ref_res = ref[0] if ref else float(np.median([b[0] for b in boots]))
kept = [b for b in boots if b[0] <= TRIM * ref_res]
drop = len(boots) - len(kept)
if len(kept) < 5:                      # trimming too aggressive to be useful
    kept, drop = boots, 0
    print('NOTE: trim would leave <5 replicates, so nothing was trimmed.')

res = np.array([b[0] for b in kept], float)
print('reference (un-jittered) SSR = %.4e' % ref_res if ref else 'no reference fit')
print('replicates: %d finished, %d kept, %d dropped as non-converged (SSR > %.1f x reference)'
      % (len(boots), len(kept), drop, TRIM))
print('kept SSR: median %.4e   range [%.4e, %.4e]' % (np.median(res), res.min(), res.max()))
print()

B = box()
hdr = '%-9s %10s %10s %10s %22s %6s %22s' % (
    'param', 'point', 'boot mean', 'boot SD', '95% CI (point +/- 1.96SD)',
    'w', 'percentile CI')
print(hdr); print('-' * len(hdr))

# keep only replicates that reported every parameter, so all arrays align
kept = [b for b in kept if all(n in b[1] for n in FIT6)] or kept

summary, allv = {}, {}
for n in FIT6:
    a = np.array([b[1][n] for b in kept if n in b[1]], float)
    if a.size == 0:
        print('%-9s %s' % (n, '(not parsed)')); continue
    allv[n] = a
    mu, sd = float(a.mean()), float(a.std(ddof=1)) if a.size > 1 else 0.0
    pt = ref[1].get(n, mu) if ref else mu
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
    summary[n] = dict(point=pt, mean=mu, sd=sd, ci=[lo, hi], pct=[float(plo), float(phi)],
                      rel_width=w, n=int(a.size))
    print('%-9s %10.4g %10.4g %10.4g   [%8.4g, %8.4g] %6.2f   [%8.4g, %8.4g]%s'
          % (n, pt, mu, sd, lo, hi, w, plo, phi, flag))
print('-' * len(hdr))
print('w = CI width / point estimate.  w < 1 = well determined.')
print('!bound = at least one replicate sat on the edge of the search box;')
print('         that parameter\'s interval is set by the box, not by the data.')
print()
print('in physical units:')
for n in ('KZP_MULT', 'KPR_MULT'):
    if n in summary:
        s, f = summary[n], (0.03 if n == 'KZP_MULT' else 0.01)
        nm = 'kzp' if n == 'KZP_MULT' else 'KPR'
        print('  %-4s = %.4g  [%.4g, %.4g]   %s'
              % (nm, s['point'] * f, s['ci'][0] * f, s['ci'][1] * f,
                 '(uM s)^-1' if n == 'KZP_MULT' else 's^-1'))
print()
print('v77 (4-param fit) for reference only -- different model, not a target:')
print('  ' + '  '.join('%s=%g' % (k, v) for k, v in V77.items()))

json.dump(dict(n_kept=len(kept), n_dropped=drop, ref_ssr=ref_res, params=summary),
          open(os.path.join(ROOT, 'ci_summary_6param.json'), 'w'), indent=2)

# ---- figures ----
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

fig, axes = plt.subplots(2, 3, figsize=(16, 8))
for ax, n in zip(axes.ravel(), FIT6):
    a = allv.get(n)
    if a is None: ax.axis('off'); continue
    s = summary[n]
    ax.hist(a, bins=max(6, len(a) // 3), color='tab:green', alpha=0.75)
    ax.axvspan(s['ci'][0], s['ci'][1], color='tab:blue', alpha=0.12)
    ax.axvline(s['ci'][0], color='k', ls='--', lw=1)
    ax.axvline(s['ci'][1], color='k', ls='--', lw=1)
    ax.axvline(s['point'], color='tab:red', lw=2.5)
    ax.set_title('%s   %.4g [%.4g, %.4g]\n%s' % (n, s['point'], s['ci'][0], s['ci'][1],
                                                 UNITS.get(n, '')), fontsize=9)
fig.suptitle('pZAP 6-parameter bootstrap, n=%d kept   red = point estimate, '
             'band = point +/- 1.96 SD' % len(kept))
fig.tight_layout(rect=[0, 0, 1, 0.94])
p1 = os.path.join(ROOT, 'ci_hist_6param.png'); fig.savefig(p1, dpi=110); plt.close(fig)

pairs = [('lig0', 'ZAP0'), ('ZAP0', 'SYK0'), ('ZAP0', 'KZP_MULT'),
         ('lig0', 'KZP_MULT'), ('kdl0', 'KPR_MULT'), ('SYK0', 'KZP_MULT')]
fig, axes = plt.subplots(1, len(pairs), figsize=(3.1 * len(pairs), 3.3))
for ax, (x, y) in zip(np.atleast_1d(axes), pairs):
    ax_ = ax
    if x not in allv or y not in allv or len(allv[x]) != len(allv[y]):
        ax_.axis('off'); continue
    ax_.scatter(allv[x], allv[y], s=30, alpha=0.75, color='tab:blue')
    if ref:
        ax_.scatter([ref[1].get(x)], [ref[1].get(y)], marker='*', s=220,
                    color='tab:red', zorder=5)
    r = np.corrcoef(allv[x], allv[y])[0, 1] if len(allv[x]) > 2 else float('nan')
    ax_.set_xlabel(x, fontsize=8); ax_.set_ylabel(y, fontsize=8)
    ax_.set_title('r = %+.2f' % r, fontsize=9)
fig.suptitle('pairwise: a tilted cloud (|r| near 1) means only the COMBINATION '
             'is constrained; red star = point estimate', fontsize=10)
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

python $F/bootstrap_pzap.py submit 0 "$NP"     # 0 empirical -> emp000 is the reference only
echo
echo "$((NP + 1)) jobs: 1 reference + $NP parametric replicates, all 6 parameters fitted."
echo "Each is 16 particles x 20 iterations, warm-started in a 0.45 log10 box:"
echo "expect roughly 1.5-2.5 h per job, so 1-2 queue waves."
echo
echo "Check progress:   squeue -u \$USER"
echo "Partial results:  bash \$0 report      (works on however many have landed)"
