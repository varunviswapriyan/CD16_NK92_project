#!/bin/bash
# =============================================================================
# finish.sh -- everything still outstanding for Indrani's five points, in one
#              submit and one report.
#
#   bash ~/CD16_NK92_project/filesCC/finish.sh submit    # ~3.5 h of cluster time
#   bash ~/CD16_NK92_project/filesCC/finish.sh status
#   bash ~/CD16_NK92_project/filesCC/finish.sh report     # all five answers
#
# WHAT IT SUBMITS  (two independent sets, different directories, no interference)
#   A. 18 jobs  -- the last three kzp columns of the (ZAP0, kzp) map.  The first
#                  42 cells are already done; 4 of 6 rows had their minimum on
#                  the right edge, so the map needs kzp out to 2.38 before the
#                  trade-off question is closed.
#   B. 51 jobs  -- a 50-sample bootstrap of the Ca variant 'orig', i.e. HER
#                  model exactly (C1 C2 g k3 k4, no S).  The reconverge showed
#                  orig fits as well as anything, so that is what should be
#                  reported -- but orig currently has NO confidence intervals,
#                  and the Ca numbers already sent were lin1's.  Each sample is
#                  warm-started from orig's converged reference (CA_WARM_FROM),
#                  which is the fix that rescued that variant.
#
# WHAT IT DOES NOT COVER
#   Point 5 (extra Ca replicates) is blocked on Oscar's .fcs files.
#   Point 2's "look into the literature for bounds" is a reading task, not a run.
# =============================================================================
set -o pipefail
MODE=${1:-submit}
F=~/CD16_NK92_project/filesCC
CADIR=$HOME/Ca_fit_c02

echo "setting up the job environment..."
[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh 2>/dev/null
module load Miniconda3/4.9.2 >/dev/null 2>&1
CONDA_SH=/gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh
[ -f "$CONDA_SH" ] && { . "$CONDA_SH" >/dev/null 2>&1; conda activate CD16_v2 >/dev/null 2>&1; }
PY="${PZ_PY:-}"
if [ -z "$PY" ]; then
  for c in "${CONDA_PREFIX:+$CONDA_PREFIX/bin/python}" python python3; do
    [ -n "$c" ] || continue
    q=$(command -v "$c" 2>/dev/null || { [ -x "$c" ] && echo "$c"; })
    [ -n "$q" ] || continue
    "$q" -c "import numpy,pandas" >/dev/null 2>&1 && { PY="$q"; break; }
  done
fi
[ -n "$PY" ] || { echo "ERROR: no python with numpy/pandas"; exit 1; }
echo "using $PY"; echo

case "$MODE" in submit|status|report) ;; *) echo "use submit|status|report"; exit 2 ;; esac

# ------------------------------------------------------------------- status --
if [ "$MODE" = status ]; then
  n=0
  for d in "$HOME"/boot_pzap_km*_pinZAP0_KZP_MULT; do
    [ -d "$d" ] || continue
    f="$d/emp000/estimate_params_pzap_cleaned_up/analysis_param_residue.dat"
    [ -f "$f" ] && grep -qi '^linear' "$f" 2>/dev/null && n=$((n+1))
  done
  echo "  A. kzp map        : $n of 60 cells finished"
  c=$(ls "$CADIR/out_boot_ca2/orig"/boot_*.json 2>/dev/null | wc -l)
  echo "  B. Ca orig bootstrap: $c of 51 samples written"
  echo
  echo "  queue: $(squeue -u "$USER" -h 2>/dev/null | wc -l) job(s)"
  exit 0
fi

# ------------------------------------------------------------------- submit --
if [ "$MODE" = submit ]; then
  nrun=$(squeue -u "$USER" -h 2>/dev/null | wc -l)
  if [ "${nrun:-0}" -gt 0 ] && [ "${FORCE:-0}" != "1" ]; then
    echo "STOP: $nrun job(s) already queued.  Submitting now can delete their"
    echo "      working directories.  Wait, or FORCE=1 bash \$0 submit"
    exit 1
  fi

  echo "############ A. kzp map, last three columns (18 jobs) ############"
  KGRID="1.30 1.60 1.90" bash "$F/zap_kzp_map.sh" submit || {
    echo "kzp map submission failed -- read the error above"; exit 1; }

  echo
  echo "############ B. Ca bootstrap of HER model, 'orig' (51 jobs) ############"
  echo "  each sample warm-started from orig's converged reference"
  [ -d "$CADIR" ] || { echo "ERROR: $CADIR not found"; exit 1; }
  ( cd "$CADIR" && CA_WARM_FROM=orig CA_STARTS=4 \
      "$PY" "$F/bootstrap_ca2.py" submit 50 orig | tail -3 )

  echo
  echo "=================================================="
  echo "  both sets submitted.  roughly 3-4 h of cluster time."
  echo "  later:   bash \$0 status      then      bash \$0 report"
  echo "=================================================="
  exit 0
fi

# ------------------------------------------------------------------- report --
echo "#############################################################"
echo "#   ANSWERS TO INDRANI'S FIVE POINTS                        #"
echo "#############################################################"

echo
echo "============================================================="
echo " POINT 1.  Does the SYK binding rate also need estimating?"
echo "============================================================="
cat <<'EOF'
  No -- and only one rate is being estimated, not two.

  Indrani's own clarification (15 Sep): "kzp is the ITAM-bound ZAP70 and SYK
  phosphorylation rate by LCK."  So kzp is a single PHOSPHORYLATION rate that
  already covers both ZAP70 and SYK.  The ZAP and SYK BINDING rates are still
  the PNAS 1995 values used in the JI paper and were not re-estimated.
EOF
echo "  fitted value:"
BOOT_TAG=_conv "$PY" "$F/ci_report.py" 2>/dev/null \
  | grep -E '^ *(kzp|KPR) *=' | sed 's/^/    /' \
  || echo "    (run ci_report.py on boot_pzap_conv for the number)"

echo
echo "============================================================="
echo " POINT 2.  Are ZAP/SYK low, and does higher ZAP lower the rate?"
echo "============================================================="
bash "$F/zap_kzp_map.sh" report 2>/dev/null | sed -n '/POINT 2/,$p' | sed 's/^/  /'
cat <<'EOF'

  NOTE on the half of this point that is not a computation: Indrani also asked
  that we "look into the literature to determine appropriate biological bounds."
  That has not been done.  The map above stands on its own, but the literature
  values would still be useful for setting the bounds.
EOF

echo
echo "============================================================="
echo " POINT 3.  Was the Ca model changed?  (S, and the missing k3/k4)"
echo "============================================================="
if [ -d "$CADIR" ]; then
  ( cd "$CADIR" && bash "$F/indrani_q.sh" refs 2>/dev/null \
      | sed -n '/variant/,$p' | sed 's/^/  /' )
fi
echo
echo "  --- confidence intervals for HER model (orig), 50 bootstrap samples ---"
( cd "$CADIR" 2>/dev/null && "$PY" - <<'PYCA'
import json, glob, os
import numpy as np
rows = []
for f in sorted(glob.glob('out_boot_ca2/orig/boot_*.json')):
    try:
        r = json.load(open(f))
    except Exception:
        continue
    if r.get('ok') and r.get('sample', 0) > 0:
        rows.append(r)
if not rows:
    print('    no bootstrap samples yet -- run `bash finish.sh status`')
else:
    ss = np.array([r['ssr'] for r in rows], float)
    keep = [r for r, s in zip(rows, ss) if s < 5 * np.median(ss)]
    print('    n = %d samples (%d kept after outlier cut), k4 and S as in her model'
          % (len(rows), len(keep)))
    print('    %-5s %12s %26s %7s' % ('', 'estimate', '95% CI', 'w'))
    for p in sorted(keep[0]['params']):
        x = np.array([r['params'][p] for r in keep if p in r['params']], float)
        x = x[np.isfinite(x)]
        if x.size < 3:
            continue
        lo, med, hi = np.percentile(x, [2.5, 50, 97.5])
        print('    %-5s %12.4g %26s %7.2f'
              % (p, med, '[%.4g, %.4g]' % (lo, hi), (hi - lo) / abs(med) if med else float('inf')))
    print('    w = (97.5%-2.5%)/median;  w < 1 means well determined.')
PYCA
)
cat <<'EOF'

  READ THIS AS: her model as written (C1, C2, g, k3, k4; no S) fits the data as
  well as every variant we tried.  S is NOT needed and should be dropped.  The
  Ca numbers sent earlier (C1 131.5, C2 0.804, g 0.00504, S 32.1) were from the
  reduced 'lin1' variant and are SUPERSEDED by the orig values above.
EOF

echo
echo "============================================================="
echo " POINT 4.  How was the pZAP bootstrap done?"
echo "============================================================="
cat <<'EOF'
  Resampled across DAYS, exactly as asked.  Lilly ran all three cell lines on
  the same three days (26 May, 28 May, 3 Jun), so a day is resampled as a unit
  and the same day-set is applied to every line.  Drawing 3 days from 3 with
  replacement has exactly 10 distinct outcomes, so rather than drawing at
  random we fit ALL TEN and weight each by its exact multinomial probability.
  That gives the bootstrap distribution with no Monte-Carlo error.  For each
  draw the mean over the resampled days is taken and the model refit.
EOF
echo
echo "  --- pZAP confidence intervals (the reportable table) ---"
BOOT_TAG=_conv "$PY" "$F/ci_report.py" 2>/dev/null | sed -n '/REPORTABLE/,$p' | sed 's/^/  /' \
  || echo "  (could not read boot_pzap_conv)"

echo
echo "============================================================="
echo " POINT 5.  Additional Ca experiments and replicates"
echo "============================================================="
cat <<'EOF'
  Blocked, waiting on Oscar.  Needed: the raw .fcs for Experiment 1 Replicate 2
  (truncate at 280 s) and Experiment 2 Replicate 1.  Once those arrive they go
  through the same FlowKit notebook in ca_data to get mean and SE, and with a
  third replicate the Ca intervals can use the same day-resampling bootstrap as
  pZAP instead of the parametric draw from mean +/- SE.
EOF
echo "  .fcs files currently on the cluster:"
find "$HOME" -maxdepth 5 -name '*.fcs' 2>/dev/null | xargs -r -n1 basename \
  | sort -u | sed 's/^/    /' || echo "    none found"

echo
echo "#############################################################"
echo "#   END                                                      #"
echo "#############################################################"
