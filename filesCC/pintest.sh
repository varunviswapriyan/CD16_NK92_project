#!/bin/bash
# =============================================================================
# pintest.sh -- "is this parameter redundant, or just under-determined?"
#
#   bash ~/CD16_NK92_project/filesCC/pintest.sh submit    # 6 pins x 4 refits = 24 jobs
#   bash ~/CD16_NK92_project/filesCC/pintest.sh submit SYK0 lig0     # only these
#   bash ~/CD16_NK92_project/filesCC/pintest.sh report
#   bash ~/CD16_NK92_project/filesCC/pintest.sh status
#
# THE TEST
# --------
# A wide confidence interval has three possible causes, and only one of them
# means the parameter is useless:
#
#   REDUNDANT      the parameter trades off exactly against another, so the
#                  model output depends only on a combination.  Holding it
#                  fixed costs NOTHING: the SSR is unchanged.  -> drop it.
#   UNDER-DETERMINED  the parameter does move the output, but 9 noisy points
#                  cannot pin it.  Holding it fixed makes the fit WORSE.
#                  -> keep it, report the wide CI, get more data.
#   BOX-LIMITED    the search hit a bound.  The CI is not a CI.  -> widen.
#
# Fixing a parameter at its fitted value and refitting the rest separates the
# first two cleanly.  This is exactly the test that produced the Ca result:
# 'lin1' fixed k4=1 and reached SSR 5.2101e4 against the 6-parameter 5.21e4 --
# identical, which proved C1 and k4 were one parameter wearing two hats.
#
# WHY THE PIN VALUE MATTERS
# -------------------------
# The earlier ZAP0 pinning test failed (SSR 0.417 -> 1.09) and I read that as
# "ZAP0 is not redundant".  That reading was wrong: it was pinned at the value
# that was optimal for the 4-parameter v77, not for the current 6-parameter
# model.  This script pins at the arm A 6-parameter point estimate, via
# BOOT_CENTRE, so the test measures redundancy and nothing else.
# =============================================================================
set -o pipefail

MODE=${1:-help}; shift 2>/dev/null
F=~/CD16_NK92_project/filesCC
PY=$(command -v python3 || command -v python)
if [ -z "$PY" ]; then
  echo "ERROR: no python on PATH.  module load Miniconda3/4.9.2 first."; exit 1
fi

# arm A 6-parameter point estimates, log10 -- the ONLY correct pin values
#   lig0 88.28  kd10 0.002471  ZAP0 173  SYK0 36.77  KZP_MULT 2.422  KPR_MULT 0.5554
CENTRE="lig0=1.9459,kd10=-2.6069,ZAP0=2.2380,SYK0=1.5654,KZP_MULT=0.3841,KPR_MULT=-0.2554"
BASE_SSR=${BASE_SSR:-5.0698e-03}        # arm A reference, unpinned, for comparison

export BOOT_CENTRE="$CENTRE"
export BOOT_TMAX=300
export BOOT_NREPS=${BOOT_NREPS:-6}
export BOOT_HALF=${BOOT_HALF:-0.25}
export BOOT_PARTICLES=${BOOT_PARTICLES:-24}
export BOOT_ITERS=${BOOT_ITERS:-20}
NREF=${NREF:-3}        # refits per pin; >1 so optimizer scatter does not fool us

ALL="lig0 kd10 ZAP0 SYK0 KZP_MULT KPR_MULT"
PARAMS="$*"; [ -z "$PARAMS" ] && PARAMS="$ALL"

root_for () { echo "$HOME/boot_pzap_pin$1"; }

submit () {
  echo "pinning at the arm A 6-parameter estimates, $NREF refits each"
  echo "baseline (nothing pinned) SSR = $BASE_SSR"
  echo
  local p root runs r n=0
  for p in $PARAMS; do
    root=$(root_for "$p")
    echo "---- $p  ->  $root"
    rm -f "$root"/run_*.sh
    # 0 day-sets, 0 parametric, NREF refits of the REAL data with $p held fixed
    if ! BOOT_PIN="$p" "$PY" "$F/bootstrap_pzap.py" build 0 0 "$NREF"; then
      echo "   build FAILED for $p -- skipping"; continue
    fi
    if ! BOOT_PIN="$p" "$PY" "$F/verify_boot.py" "$root" >/dev/null 2>&1; then
      echo "   verify FAILED for $p -- skipping (rerun verify_boot.py $root to see why)"
      continue
    fi
    shopt -s nullglob; runs=("$root"/run_*.sh); shopt -u nullglob
    for r in "${runs[@]}"; do sbatch "$r" >/dev/null 2>&1 && n=$((n+1)); done
    echo "   submitted ${#runs[@]}"
  done
  echo
  echo "$n jobs submitted.  later:  bash \$0 report"
}

status () {
  local p root d dat done tot
  for p in $PARAMS; do
    root=$(root_for "$p"); done=0; tot=0
    [ -d "$root" ] || { printf '  %-10s not started\n' "$p"; continue; }
    for d in "$root"/emp[0-9][0-9][0-9] "$root"/nul[0-9][0-9][0-9]; do
      [ -d "$d" ] || continue; tot=$((tot+1))
      dat="$d/estimate_params_pzap_cleaned_up/analysis_param_residue.dat"
      [ -f "$dat" ] && grep -qi '^linear' "$dat" 2>/dev/null && done=$((done+1))
    done
    printf '  %-10s %s of %s refits finished\n' "$p" "$done" "$tot"
  done
}

report () {
  BASE_SSR="$BASE_SSR" PARAMS="$PARAMS" "$PY" - <<'PYEOF'
import os, glob, numpy as np
HOME = os.path.expanduser('~'); SUB = 'estimate_params_pzap_cleaned_up'
base = float(os.environ['BASE_SSR'])
print('baseline SSR (nothing pinned, arm A reference) = %.4e\n' % base)
print('%-10s %4s %12s %9s   %s' % ('pinned', 'n', 'median SSR', 'x base', 'verdict'))
print('-' * 72)
any_row = False
for p in os.environ['PARAMS'].split():
    root = os.path.join(HOME, 'boot_pzap_pin' + p)
    ss = []
    for d in sorted(glob.glob(os.path.join(root, '*[0-9]'))):
        f = os.path.join(d, SUB, 'analysis_param_residue.dat')
        if not os.path.exists(f):
            continue
        for ln in open(f, errors='ignore'):
            if ln.strip().lower().startswith('residue ='):
                try: ss.append(float(ln.split('=', 1)[1]))
                except Exception: pass
                break
    if not ss:
        print('%-10s %4s %12s %9s   %s' % (p, '-', '-', '-', 'no finished refits yet'))
        continue
    any_row = True
    med = float(np.median(ss)); ratio = med / base
    if   ratio < 1.15: v = 'REDUNDANT -- fixing it costs nothing, drop it'
    elif ratio < 1.5:  v = 'mostly redundant -- worth fixing to shrink the model'
    elif ratio < 3.0:  v = 'real parameter, under-determined -- keep, report wide CI'
    else:              v = 'strongly identified -- the wide CI is NOT redundancy'
    print('%-10s %4d %12.4e %9.2f   %s' % (p, len(ss), med, ratio, v))
print('-' * 72)
if any_row:
    print('ratio = median SSR with the parameter FIXED / SSR with it free.')
    print('~1.0 means the data never needed it.  Large means it does real work')
    print('and its wide CI is a data-quantity problem, not a redundancy problem.')
    print('\nCompare with the Ca precedent: lin1 fixed k4 and the SSR did not move')
    print('(ratio 1.000), which is why dropping it there was correct.')
PYEOF
}

case "$MODE" in
  submit) submit ;;
  status) status ;;
  report) report ;;
  *)      sed -n '2,40p' "$0" ;;
esac
