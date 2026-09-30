#!/bin/bash
# =============================================================================
# das_rerun.sh -- Dr. Das point 2: "the estimated SYK and ZAP concentrations are
# low compared with values reported in the T cell literature... A higher ZAP
# concentration may lead to a lower estimated ZAP binding rate."
#
#   bash ~/CD16_NK92_project/filesCC/das_rerun.sh submit
#   bash ~/CD16_NK92_project/filesCC/das_rerun.sh status
#   bash ~/CD16_NK92_project/filesCC/das_rerun.sh report
#
# WHAT THIS RUN DOES
#   Raises the ZAP0 ceiling  1585 -> 31600  (log10 3.2 -> 4.5)
#   Raises the SYK0 ceiling   316 -> 10000  (log10 2.5 -> 4.0)
#   Lowers the ZAP0 floor     100 -> 50     so the old floor stops binding
#   and searches ZAP0 and SYK0 over that FULL range rather than warm-starting
#   them near their old values.  Everything else is unchanged.
#
#   Widening bounds alone would have done nothing: the search box is the v77
#   optimum +/- BOOT_HALF clipped to the bounds, so ZAP0 would have stayed
#   within +/-0.25 log10 of 138 no matter how wide the bounds were.  BOOT_FREE
#   is what lets it move.
#
# WHAT IT ANSWERS
#   * Does the fit WANT a higher ZAP0/SYK0, or does it stay near 138 / 58.6?
#     If it stays, the data prefers low values and that is the finding to
#     report, whatever the literature says.
#   * Is Dr. Das right that higher ZAP0 pulls the binding rate down?  If ZAP0
#     rises and kzp falls, that is the predicted trade-off, and it means the
#     two are only identifiable as a product -- which sets up the argument for
#     fixing ZAP0 at a literature value rather than fitting it.
#   * Full bootstrap, so the new numbers arrive with CIs, not as point values.
# =============================================================================
set -o pipefail
MODE=${1:-submit}
F=~/CD16_NK92_project/filesCC

export BOOT_TAG=${BOOT_TAG:-_wide}
export BOOT_TMAX=300
export BOOT_NREPS=${BOOT_NREPS:-6}
export BOOT_HALF=${BOOT_HALF:-0.25}
export BOOT_PARTICLES=${BOOT_PARTICLES:-40}
export BOOT_ITERS=${BOOT_ITERS:-50}
export BOOT_WALLTIME=${BOOT_WALLTIME:-16:00:00}
export BOOT_BOUNDS=${BOOT_BOUNDS:-"ZAP0=1.7:4.5,SYK0=0.5:4.0"}
export BOOT_FREE=${BOOT_FREE:-"ZAP0,SYK0"}
NEMP=${NEMP:-10}
NNUL=${NNUL:-3}
PROOT=$HOME/boot_pzap${BOOT_TAG}

PY=$(command -v python3 || command -v python)
[ -n "$PY" ] || { echo "ERROR: no python; module load Miniconda3/4.9.2"; exit 1; }
case "$MODE" in submit|status|report) ;; *) echo "use submit|status|report"; exit 2 ;; esac

if [ "$MODE" = status ]; then
  d=0; t=0
  for x in "$PROOT"/emp[0-9][0-9][0-9] "$PROOT"/nul[0-9][0-9][0-9]; do
    [ -d "$x" ] || continue; t=$((t+1))
    f="$x/estimate_params_pzap_cleaned_up/analysis_param_residue.dat"
    [ -f "$f" ] && grep -qi '^linear' "$f" 2>/dev/null && d=$((d+1))
  done
  printf '  widened-bounds run   %s of %s fits finished\n' "$d" "$t"
  echo; echo "queue:"
  squeue -u "$USER" -h -o '%j' 2>/dev/null | sed 's/[0-9]*$//' | sort | uniq -c \
    | awk '{printf "  %5s  %s\n", $1, $2}'
  exit 0
fi

if [ "$MODE" = report ]; then
  echo "############ WIDENED ZAP0 / SYK0 BOUNDS ############"
  "$PY" "$F/ci_report.py" || echo "(not enough finished samples yet)"
  echo
  echo "############ PREVIOUS RUN, for comparison ##########"
  BOOT_TAG=_conv BOOT_BOUNDS= BOOT_FREE= "$PY" "$F/ci_report.py" || true
  cat <<'EOF'

HOW TO READ THE COMPARISON
  ZAP0 rises a lot and kzp falls   -> Dr. Das is right; the two trade off as a
                                      product, and ZAP0 should be FIXED at a
                                      literature value rather than fitted.
  ZAP0 stays near 138              -> the data genuinely prefers a low value.
                                      Report that, and the gap with the T cell
                                      literature becomes a question about the
                                      model's units or about NK92 vs T cells.
  SSR unchanged either way         -> the two are interchangeable; the fit
                                      cannot tell them apart at all.
EOF
  exit 0
fi

echo "pZAP refit with widened ZAP0 / SYK0"
echo "  root     : $PROOT"
echo "  ZAP0     : 50 - 31600      (was 100 - 1585)"
echo "  SYK0     : 3 - 10000       (was 3 - 316)"
echo "  both searched over the full range, not warm-started"
echo "  budget   : $BOOT_PARTICLES x $BOOT_ITERS   N_REPS=$BOOT_NREPS   0-300 s"
echo

rm -f "$PROOT"/run_*.sh
echo "== 1/3  building =="
"$PY" "$F/bootstrap_pzap.py" build "$NEMP" 0 "$NNUL" || exit 1

echo
echo "== 2/3  verifying =="
"$PY" "$F/verify_boot.py" "$PROOT" || { echo; echo "ABORTED -- nothing submitted."; exit 1; }

echo
echo "== 3/3  submitting =="
shopt -s nullglob; runs=("$PROOT"/run_*.sh); shopt -u nullglob
if [ ${#runs[@]} -eq 0 ]; then
  echo "nothing to submit."
else
  for r in "${runs[@]}"; do
    printf '%s: ' "$(basename "$r" .sh | sed 's/^run_//')"
    sbatch "$r" 2>&1 | tail -1
  done
  echo; echo "submitted ${#runs[@]} jobs.   later:  bash \$0 report"
fi
