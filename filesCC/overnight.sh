#!/bin/bash
# =============================================================================
# overnight.sh -- the run that fixes the point-estimate-outside-CI problem.
#
#   bash ~/CD16_NK92_project/filesCC/overnight.sh submit   # pZAP + Ca
#   bash ~/CD16_NK92_project/filesCC/overnight.sh status
#   bash ~/CD16_NK92_project/filesCC/overnight.sh report
#   bash ~/CD16_NK92_project/filesCC/overnight.sh stopa2   # free the queue first
#
# WHY THIS RUN, AND NOT MORE SAMPLES
# ----------------------------------
# The 10 day-resamples ARE the complete empirical bootstrap; more samples of
# them cannot exist.  The defect is not sample count, it is convergence:
#
#   reference SSR              5.07e-03
#   median bootstrap SSR       3.20e-02      6.3x worse
#   spread across samples      57x
#   one sample reached         4.19e-03      BETTER than the reference
#
# Resampling the same three days cannot make the achievable SSR vary 57-fold.
# The refits are stopping short of good optima at 24 particles x 20 iterations
# (480 evaluations for 6 parameters), and where they stop is systematically
# offset from the reference fit.  That offset is what pushes the bootstrap
# distribution off the point estimate.
#
# Calcium is the control: its median/reference SSR ratio is 1.31, because its
# pipeline runs a local polish after the global search.  Same statistics,
# different convergence, and only pZAP has the containment problem.
#
# So this run raises the budget ~4x and leaves everything else identical.
#
# WHAT I CAN AND CANNOT PROMISE
# -----------------------------
# I cannot guarantee containment -- that depends on whether the offset is
# convergence (fixable here) or genuine bootstrap bias (not).  What this run
# guarantees is that we will KNOW WHICH, because ci_report.py checks the
# convergence criteria separately from the intervals.  If convergence passes
# and containment still fails, the cause is bias, and the fallback is to quote
# the bootstrap median, which the percentile interval contains by construction.
# ci_report.py prints that table too, so either way you have a reportable
# answer in the morning.
# =============================================================================
set -o pipefail
MODE=${1:-submit}
F=~/CD16_NK92_project/filesCC
CADIR=$HOME/Ca_fit_c02

# ---- pZAP: same window, same N_REPS, same box -- only the budget changes ----
export BOOT_TAG=${BOOT_TAG:-_conv}
export BOOT_TMAX=300
export BOOT_NREPS=${BOOT_NREPS:-6}
export BOOT_HALF=${BOOT_HALF:-0.25}
export BOOT_PARTICLES=${BOOT_PARTICLES:-40}    # was 24
export BOOT_ITERS=${BOOT_ITERS:-50}            # was 20  -> 2000 evals, ~4.2x
export BOOT_WALLTIME=${BOOT_WALLTIME:-16:00:00}
NEMP=${NEMP:-10}
NNUL=${NNUL:-3}        # refits of the REAL data: is the reference reproducible?
NCA=${NCA:-50}
CAVAR=${CAVAR:-lin1}
PROOT=$HOME/boot_pzap${BOOT_TAG}

PY=$(command -v python3 || command -v python)
[ -n "$PY" ] || { echo "ERROR: no python on PATH; module load Miniconda3/4.9.2"; exit 1; }

case "$MODE" in submit|status|report|stopa2) ;; *)
  echo "ERROR: use submit, status, report or stopa2"; exit 2 ;; esac

# ---------------------------------------------------------------- stopa2 ----
if [ "$MODE" = stopa2 ]; then
  echo "Arm A2 (N_REPS=16) tested whether more replicates tighten SYK0/KPR."
  echo "The arm A vs B comparison already answered that: they do not."
  n=$(squeue -u "$USER" -h -o '%i %o' 2>/dev/null | grep -c 'boot_pzap_n16' || true)
  echo "jobs still queued/running for boot_pzap_n16: ${n:-0}"
  ids=$(squeue -u "$USER" -h -o '%i %o' 2>/dev/null | awk '/boot_pzap_n16/{print $1}')
  if [ -z "$ids" ]; then echo "nothing to cancel."; exit 0; fi
  echo "$ids" | xargs -r scancel && echo "cancelled."
  exit 0
fi

# ---------------------------------------------------------------- status ----
if [ "$MODE" = status ]; then
  d=0; t=0
  for x in "$PROOT"/emp[0-9][0-9][0-9] "$PROOT"/nul[0-9][0-9][0-9]; do
    [ -d "$x" ] || continue; t=$((t+1))
    f="$x/estimate_params_pzap_cleaned_up/analysis_param_residue.dat"
    [ -f "$f" ] && grep -qi '^linear' "$f" 2>/dev/null && d=$((d+1))
  done
  printf '  pZAP convergence run   %s of %s fits finished\n' "$d" "$t"
  if [ -d "$CADIR/out_boot_ca2/$CAVAR" ]; then
    printf '  Ca %-20s %s samples written\n' "$CAVAR" \
      "$(ls "$CADIR/out_boot_ca2/$CAVAR"/boot_*.json 2>/dev/null | wc -l)"
  fi
  echo; echo "queue:"
  squeue -u "$USER" -h -o '%j' 2>/dev/null | sed 's/[0-9]*$//' | sort | uniq -c \
    | awk '{printf "  %5s  %s\n", $1, $2}'
  echo "  ----- $(squeue -u "$USER" -h 2>/dev/null | wc -l) jobs total"
  exit 0
fi

# ---------------------------------------------------------------- report ----
if [ "$MODE" = report ]; then
  echo "##################  pZAP -- convergence run  ##################"
  "$PY" "$F/ci_report.py" || echo "(not enough finished samples yet)"
  echo
  echo "##################  pZAP -- previous run, for comparison  #####"
  BOOT_TAG='' "$PY" "$F/ci_report.py" || true
  echo
  echo "##################  Calcium  ##################################"
  if [ -d "$CADIR" ]; then cd "$CADIR" && "$PY" "$F/bootstrap_ca2.py" report; fi
  exit 0
fi

# ---------------------------------------------------------------- submit ----
echo "pZAP convergence run"
echo "  root       : $PROOT"
echo "  budget     : $BOOT_PARTICLES particles x $BOOT_ITERS iterations"
echo "               (was 24 x 20 = 480 evals; now $((BOOT_PARTICLES*BOOT_ITERS)))"
echo "  walltime   : $BOOT_WALLTIME   N_REPS=$BOOT_NREPS   window 0-300 s"
echo "  samples    : $NEMP day-sets + reference + $NNUL refits of the real data"
echo

rm -f "$PROOT"/run_*.sh
echo "== 1/3  building =="
"$PY" "$F/bootstrap_pzap.py" build "$NEMP" 0 "$NNUL" || exit 1

echo
echo "== 2/3  verifying =="
if ! "$PY" "$F/verify_boot.py" "$PROOT"; then
  echo; echo "ABORTED -- nothing submitted."; exit 1
fi

echo
echo "== 3/3  submitting pZAP =="
shopt -s nullglob; runs=("$PROOT"/run_*.sh); shopt -u nullglob
if [ ${#runs[@]} -eq 0 ]; then
  echo "nothing to submit."
else
  for r in "${runs[@]}"; do
    printf '%s: ' "$(basename "$r" .sh | sed 's/^run_//')"
    sbatch "$r" 2>&1 | tail -1
  done
  echo "submitted ${#runs[@]} pZAP jobs."
fi

echo
echo "== Ca: $((NCA+1)) jobs ($CAVAR) =="
if [ -d "$CADIR" ]; then
  cd "$CADIR" && "$PY" "$F/bootstrap_ca2.py" submit "$NCA" "$CAVAR"
else
  echo "$CADIR not found -- skipped"
fi

cat <<EOF

==============================================================
In the morning:

    bash \$0 status
    bash \$0 report

ci_report.py checks three things IN ORDER, and says which failed:

  1. median bootstrap SSR / reference SSR  <= 1.5   (was 6.3)
  2. SSR spread across samples             <= 10x   (was 57x)
  3. point estimate inside its interval    for all 6

If 1 and 2 pass but 3 still fails, the offset is real bootstrap bias
rather than convergence, and the report prints the bootstrap-median
table as the fallback -- contained by construction.
==============================================================
EOF
