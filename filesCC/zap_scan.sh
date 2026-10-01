#!/bin/bash
# =============================================================================
# zap_scan.sh -- Dr. Das point 2, done as a profile scan instead of an assertion.
#
#   bash ~/CD16_NK92_project/filesCC/zap_scan.sh submit
#   bash ~/CD16_NK92_project/filesCC/zap_scan.sh status
#   bash ~/CD16_NK92_project/filesCC/zap_scan.sh report
#
# THE QUESTION (Das, point 2)
#   "A higher ZAP concentration may lead to a lower estimated ZAP binding rate."
#   i.e. are ZAP0 and kzp a trade-off?  If so, the data fixes only their
#   combination, not each alone, and ZAP0 should be FIXED at a literature value
#   with kzp read off the curve -- rather than both being fitted.
#
# WHY A SCAN AND NOT just the widened-bounds run (das_rerun.sh)
#   das_rerun let ZAP0 float and it stayed ~144.  That answers "does the fit
#   WANT a high ZAP0 on its own" (no).  It does NOT test Das's trade-off, which
#   is conditional: "IF ZAP0 were higher, would kzp drop while the fit stayed
#   as good?"  To test that you must HOLD ZAP0 at a series of values and refit
#   everything else -- which is what this does.
#
# HOW
#   For each ZAP0 on a log10 grid, pin ZAP0 there (BOOT_PIN + BOOT_CENTRE) and
#   refit the other 5 parameters against the REAL (un-resampled) data.  Record
#   kzp and the SSR at each grid point.  This reuses Indrani's fitter unchanged;
#   nothing about the model is altered, only which parameter is held.
#
# HOW TO READ IT
#   SSR stays flat as ZAP0 rises, kzp falls ~1/ZAP0   -> Das is right: ZAP0 and
#       kzp are degenerate (only the product is fixed).  Report: pin ZAP0 at the
#       literature concentration, quote the kzp from the curve.
#   SSR RISES as ZAP0 leaves ~140                      -> ZAP0 is genuinely
#       identifiable and the data prefers the low value.  Report that; the gap
#       with T-cell literature is then a real NK92-vs-T-cell difference, not a
#       fitting artefact.
#
# NOTE ON THE "LITERATURE VALUE"
#   The PNAS 1995 paper (Bu/Shaw/Chan) Indrani cited is the source of the ZAP
#   BINDING RATE (ka = 6.0e6 M^-1 s^-1), NOT a ZAP concentration -- it is an in
#   vitro Biacore study with no cellular copy numbers.  So this scan does not
#   assume any particular literature ZAP0; it produces the whole kzp-vs-ZAP0
#   curve, and a concentration from a T-cell proteomics reference can be dropped
#   onto it later to read off the implied kzp.
# =============================================================================
set -o pipefail
MODE=${1:-submit}
F=~/CD16_NK92_project/filesCC
BP=$F/bootstrap_pzap.py

# log10(ZAP0) grid.  2.14 = the current fitted value (~138); top = 4.5 (~31600),
# the same ceiling das_rerun.sh used.  Evenly spaced so the trade-off, if there
# is one, shows as a straight kzp-vs-ZAP0 line on a log-log axis.
GRID=${GRID:-"2.14 2.45 2.75 3.05 3.35 3.70 4.10 4.50"}

# a warm-started 5-parameter refit is far easier than the free 6-param search,
# so a smaller budget than the convergence run is enough.  Override if needed.
export BOOT_PARTICLES=${BOOT_PARTICLES:-24}
export BOOT_ITERS=${BOOT_ITERS:-30}
export BOOT_NREPS=${BOOT_NREPS:-6}
export BOOT_WALLTIME=${BOOT_WALLTIME:-8:00:00}
export BOOT_HALF=${BOOT_HALF:-0.60}
export BOOT_TMAX=${BOOT_TMAX:-300}

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
[ -f "$BP" ] || { echo "ERROR: $BP not found -- sync filesCC first"; exit 1; }
echo "using $PY"
echo

case "$MODE" in submit|status|report) ;; *) echo "use submit|status|report"; exit 2 ;; esac

# index for a grid value: its position, zero-padded, so tags are stable
idx_of() { local i=0 v; for v in $GRID; do [ "$v" = "$1" ] && { printf '%02d' "$i"; return; }; i=$((i+1)); done; }

# ------------------------------------------------------------------- status --
if [ "$MODE" = status ]; then
  i=0
  for v in $GRID; do
    t=$(printf '%02d' "$i")
    d="$HOME/boot_pzap_zs${t}_pinZAP0/emp000/estimate_params_pzap_cleaned_up/analysis_param_residue.dat"
    if [ -f "$d" ] && grep -qi '^linear' "$d" 2>/dev/null; then st="done"; else st="...."; fi
    printf '  ZAP0(log10)=%-5s  ZAP0=%-8.0f  %s\n' "$v" "$("$PY" -c "print(10**$v)")" "$st"
    i=$((i+1))
  done
  echo; echo "queue: $(squeue -u "$USER" -h 2>/dev/null | grep -c bpemp) scan jobs running"
  exit 0
fi

# ------------------------------------------------------------------- report --
if [ "$MODE" = report ]; then
  GRID="$GRID" "$PY" - <<'PYREP'
import os, sys, glob
sys.path.insert(0, os.path.expanduser('~/CD16_NK92_project/filesCC'))
# clear any BOOT_* that would perturb the import-time config parsing
for k in list(os.environ):
    if k.startswith('BOOT_'): os.environ.pop(k)
import bootstrap_pzap as BP
import numpy as np, json

grid = os.environ['GRID'].split()
KZP_TO_PHYS = 0.03           # kzp = KZP_MULT * 0.03  (uM s)^-1
rows = []
for i, v in enumerate(grid):
    root = os.path.expanduser('~/boot_pzap_zs%02d_pinZAP0' % i)
    d = os.path.join(root, 'emp000')
    r = BP.parse_result(d)
    zap = 10.0**float(v)
    if not r:
        rows.append((zap, None, None, None)); continue
    ssr, m = r
    kzp = m.get('KZP_MULT')
    kzp = kzp*KZP_TO_PHYS if kzp is not None else None
    rows.append((zap, kzp, ssr, m))

done = [x for x in rows if x[1] is not None]
print('=' * 72)
print('POINT 2: kzp vs pinned ZAP0   (real data, other 5 params refitted)')
print('=' * 72)
print('%10s %14s %14s   %s' % ('ZAP0', 'kzp (uM s)^-1', 'SSR', 'ZAP0*kzp (product)'))
print('-' * 72)
for zap, kzp, ssr, m in rows:
    if kzp is None:
        print('%10.0f %14s %14s   %s' % (zap, '(pending)', '-', '-')); continue
    print('%10.0f %14.5g %14.5g   %14.5g' % (zap, kzp, ssr, zap*kzp))
print('-' * 72)

if len(done) >= 3:
    ssrs = np.array([x[2] for x in done], float)
    zaps = np.array([x[0] for x in done], float)
    kzps = np.array([x[1] for x in done], float)
    smin, smax = ssrs.min(), ssrs.max()
    spread = smax / smin if smin > 0 else np.inf
    # is kzp ~ 1/ZAP0 ?  fit log(kzp) vs log(ZAP0); slope -1 == pure product
    sl = np.polyfit(np.log(zaps), np.log(kzps), 1)[0]
    prods = zaps * kzps
    prod_cv = prods.std() / prods.mean() if prods.mean() else np.inf
    print()
    print('SSR range over the grid : %.4g -- %.4g   (max/min = %.2fx)' % (smin, smax, spread))
    print('slope d[log kzp]/d[log ZAP0] : %.2f   (-1.00 = exact 1/ZAP0 trade-off)' % sl)
    print('ZAP0*kzp spread (CV)     : %.1f%%   (small = product is what is fixed)' % (100*prod_cv))
    print()
    print('=' * 72)
    if spread < 1.3 and prod_cv < 0.15:
        print('VERDICT: SSR is essentially flat while kzp tracks 1/ZAP0 and the')
        print('product ZAP0*kzp is nearly constant.  Das is right -- ZAP0 and kzp')
        print('are degenerate.  Fix ZAP0 at a literature concentration and quote')
        print('the kzp the curve gives at that ZAP0.')
    elif spread < 1.3:
        print('VERDICT: SSR is flat but the product is not constant -- ZAP0 is')
        print('weakly constrained.  Safe to fix it at a literature value; report')
        print('the refitted kzp there.')
    else:
        print('VERDICT: SSR rises %.1fx as ZAP0 leaves ~140.  ZAP0 is genuinely' % spread)
        print('identifiable and the data prefers the low value.  Report that; the')
        print('gap with T-cell literature is a real NK92 difference, not a fitting')
        print('artefact.  Do NOT fix ZAP0 high -- it degrades the fit.')
    print('=' * 72)
else:
    print('\nonly %d of %d grid points finished -- run `status`, wait, re-report.'
          % (len(done), len(rows)))
PYREP
  exit 0
fi

# ------------------------------------------------------------------- submit --
# NOTE ON A BUG THAT COST A NIGHT: this loop used to pipe each point through
#   grep -E 'emp000|Submitted|ERROR|box:'
# A Python traceback contains none of those words -- "FileNotFoundError" does
# not match the case-sensitive 'ERROR' -- so when submission failed, the filter
# ate the error and the script still printed "submitted 8 scan points".  The
# queue was empty and nothing said why.  Output is now shown in full, and the
# script VERIFIES that sbatch actually returned a job id before claiming
# anything was submitted.
echo "=============== POINT 2: ZAP0 profile scan ==============="
echo "  grid log10(ZAP0): $GRID"
echo "  budget: ${BOOT_PARTICLES}x${BOOT_ITERS}  N_REPS=$BOOT_NREPS  window 0-${BOOT_TMAX}s"
echo "  at each point ZAP0 is PINNED and lig0,kd10,SYK0,KZP_MULT,KPR_MULT refit"
echo

echo "--- preflight ---"
SRCD=$HOME/NK92_fit_v77
XL=/home/gddaslab/share/Varun_Indrani/estimate_params_pzap/data/pZAP70_Tyr493_Tyr292_original_and_averages.xlsx
fail=0
[ -d "$SRCD" ]  && echo "  ok   source tree   $SRCD" || { echo "  FAIL source tree MISSING: $SRCD"; fail=1; }
[ -f "$XL" ]    && echo "  ok   day-level xlsx" || { echo "  FAIL xlsx MISSING: $XL"; fail=1; }
command -v sbatch >/dev/null && echo "  ok   sbatch on PATH" || { echo "  FAIL sbatch not found (are you on a login node?)"; fail=1; }
"$PY" -c "import sys; sys.path.insert(0,'$F'); import bootstrap_pzap" 2>&1 \
  && echo "  ok   bootstrap_pzap imports" || { echo "  FAIL bootstrap_pzap does not import (full error above)"; fail=1; }
avail=$(df -Pk "$HOME" 2>/dev/null | awk 'NR==2{print int($4/1048576)}')
[ -n "$avail" ] && echo "  note ${avail} GB free in \$HOME (each grid point copies the v77 tree)"
[ "$fail" = 0 ] || { echo; echo "PREFLIGHT FAILED -- nothing submitted."; exit 1; }
echo

i=0; nok=0; nbad=0
for v in $GRID; do
  t=$(printf '%02d' "$i")
  z=$("$PY" -c "print(round(10**$v))")
  echo "--- ZAP0(log10)=$v  (ZAP0=$z)  tag=_zs$t"
  log=$(mktemp)
  BOOT_TAG="_zs$t" BOOT_PIN="ZAP0" BOOT_CENTRE="ZAP0=$v" \
    "$PY" "$BP" submit 0 0 0 >"$log" 2>&1
  rc=$?
  jid=$(grep -o 'Submitted batch job [0-9]*' "$log" | head -1 | awk '{print $NF}')
  if [ "$rc" -ne 0 ] || [ -z "$jid" ]; then
    nbad=$((nbad+1))
    echo "    SUBMIT FAILED (exit $rc).  Full output:"
    sed 's/^/      /' "$log"
  else
    nok=$((nok+1))
    echo "    submitted job $jid"
  fi
  rm -f "$log"
  i=$((i+1))
done
echo
echo "=================================================="
echo "  submitted OK : $nok of $i"
[ "$nbad" -gt 0 ] && echo "  FAILED       : $nbad   (error text printed above each one)"
if [ "$nok" -gt 0 ]; then
  echo
  echo "  check:   squeue -u $USER"
  echo "  later:   bash \$0 status    then    bash \$0 report"
else
  echo
  echo "  NOTHING was submitted.  Read the error above before re-running."
fi
echo "=================================================="
