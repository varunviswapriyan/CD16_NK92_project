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
#
# ALLGRID is the CANONICAL list and fixes each point's tag (_zs00 .. _zs07).
# GRID is only which points to submit.  Tags come from a value's position in
# ALLGRID, never from its position in GRID -- otherwise resubmitting a subset
# like GRID="4.10 4.50" would reuse _zs00/_zs01 and overwrite finished points.
ALLGRID="2.14 2.45 2.75 3.05 3.35 3.70 4.10 4.50"
GRID=${GRID:-$ALLGRID}

# memory per point.  NFsim holds every molecule and complex explicitly, so the
# high-ZAP0 points need much more than the 64G default.
export BOOT_MEM=${BOOT_MEM:-64G}

# -----------------------------------------------------------------------------
# SCAN_FREE=1 -- THE RUN THAT CAN ACTUALLY TEST DAS'S TRADE-OFF
#
# The default scan warm-starts the five free parameters in a box of v77 +/-
# BOOT_HALF (0.6 in log10), so kzp may move only 4x either way.  But the ZAP0
# sweep spans 138 -> 5012, which is 36x (1.56 log10).  For the product
# ZAP0*kzp to stay constant -- exactly Das's hypothesis -- kzp would have to
# fall ~36x.  The box FORBIDS that, so the default run cannot see the
# trade-off even if it is real, and its rising SSR partly measures the box
# walls rather than the biology.
#
# SCAN_FREE lets every non-pinned parameter search its full range (BOOT_FREE)
# and widens KZP_MULT's range to 10^-2 .. 10^1.5, so kzp can drop ~100x.  It
# writes to its own _zf** tags so the constrained results stay intact and the
# two can be compared.
#
#   SCAN_FREE=1 bash $0 submit     then     SCAN_FREE=1 bash $0 report
#
# READ THE COMPARISON LIKE THIS
#   free-run SSR now FLAT and kzp ~ 1/ZAP0  -> the constrained verdict was a box
#        artefact; Das is right, only the product is identifiable.
#   free-run SSR STILL rises steeply        -> the trade-off genuinely does not
#        exist; ZAP0 is identifiable and NK92 really does prefer ~140.  That is
#        now a safe thing to report, because kzp was free to compensate and
#        chose not to.
# -----------------------------------------------------------------------------
if [ "${SCAN_FREE:-0}" = "1" ]; then
  PFX=zf
  export BOOT_FREE=${BOOT_FREE:-"lig0,kdl0,SYK0,KZP_MULT,KPR_MULT"}
  export BOOT_BOUNDS=${BOOT_BOUNDS:-"KZP_MULT=-2.0:1.5"}
  export BOOT_PARTICLES=${BOOT_PARTICLES:-40}   # full-range search needs more
  export BOOT_ITERS=${BOOT_ITERS:-50}
else
  PFX=zs
fi

# tag index for a log10 value = its position in ALLGRID (stable across subsets)
tag_of() {
  local k=0 w
  for w in $ALLGRID; do
    [ "$w" = "$1" ] && { printf '%02d' "$k"; return 0; }
    k=$((k+1))
  done
  # not in the canonical grid: derive a stable tag from the value itself
  printf 'v%s' "$(printf '%s' "$1" | tr -d '.')"
}

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

case "$MODE" in submit|status|report|logs) ;; *) echo "use submit|status|report|logs"; exit 2 ;; esac

# --------------------------------------------------------------------- logs --
# Why did a grid point die?  Slurm's accounting record gives the verdict
# (OUT_OF_MEMORY / TIMEOUT / FAILED + exit code) and the job log gives the text.
# Points were observed dying in ZAP0 order, highest first, which is the
# signature of NFsim running out of memory: it tracks every molecule and
# complex individually, so raising ZAP0 from 138 to 31623 can blow the species
# count up combinatorially.
if [ "$MODE" = logs ]; then
  echo "=============== slurm accounting (today) ==============="
  if command -v sacct >/dev/null; then
    sacct -u "$USER" --format=JobID%14,JobName%12,State%18,ExitCode%8,MaxRSS%10,Elapsed%10 \
          -S today 2>/dev/null | grep -Ev '\.(batch|extern)' | tail -25
  else
    echo "  sacct not available"
  fi

  echo
  echo "=============== per grid point ==============="
  i=0; noom=0; ndone=0; nfail=0
  for v in $ALLGRID; do
    t=$(tag_of "$v")
    d="$HOME/boot_pzap_${PFX}${t}_pinZAP0/emp000"
    lg="$d/emp000.log"
    z=$("$PY" -c "print(round(10**$v))" 2>/dev/null)
    printf -- '--- ZAP0=%-7s tag=_zs%s\n' "$z" "$t"
    if [ ! -d "$d" ]; then
      echo "    no directory -- never built"; nfail=$((nfail+1)); i=$((i+1)); continue
    fi
    if [ -f "$d/estimate_params_pzap_cleaned_up/analysis_param_residue.dat" ] \
       && grep -qi '^linear' "$d/estimate_params_pzap_cleaned_up/analysis_param_residue.dat" 2>/dev/null; then
      echo "    FINISHED (has a fitted result)"; ndone=$((ndone+1)); i=$((i+1)); continue
    fi
    if [ ! -f "$lg" ]; then
      echo "    no log at $lg"; nfail=$((nfail+1)); i=$((i+1)); continue
    fi
    # classify
    if grep -qiE 'out of memory|oom-kill|MemoryError|std::bad_alloc|Killed' "$lg" 2>/dev/null; then
      echo "    >>> OUT OF MEMORY <<<"; noom=$((noom+1))
    elif grep -qiE 'DUE TO TIME LIMIT|CANCELLED' "$lg" 2>/dev/null; then
      echo "    >>> HIT THE WALLTIME <<<"
    fi
    echo "    last lines of $lg:"
    tail -12 "$lg" 2>/dev/null | sed 's/^/      /'
    nfail=$((nfail+1))
    i=$((i+1))
  done

  echo
  echo "=================================================="
  echo "  finished $ndone   not finished $nfail   of $i"
  if [ "$noom" -gt 0 ]; then
    cat <<'EOF'

  DIAGNOSIS: at least one point ran out of memory.  NFsim holds every molecule
  and complex explicitly, so the high-ZAP0 points need far more RAM than the
  64G the job asks for.  Two ways forward, both fine:

    1. give the big points more memory and rerun just those
         BOOT_MEM=250G GRID="3.70 4.10 4.50" bash $0 submit
    2. accept the range that works.  Points from 138 to ~5000 already span 36x,
       which is enough to see whether SSR is flat and whether kzp tracks
       1/ZAP0.  The answer to Das's point 2 does not need the top of the grid.
EOF
  fi
  echo "=================================================="
  exit 0
fi

# index for a grid value: its position, zero-padded, so tags are stable
# ------------------------------------------------------------------- status --
if [ "$MODE" = status ]; then
  i=0
  for v in $ALLGRID; do
    t=$(tag_of "$v")
    d="$HOME/boot_pzap_${PFX}${t}_pinZAP0/emp000/estimate_params_pzap_cleaned_up/analysis_param_residue.dat"
    if [ -f "$d" ] && grep -qi '^linear' "$d" 2>/dev/null; then st="done"; else st="...."; fi
    printf '  ZAP0(log10)=%-5s  ZAP0=%-8.0f  %s\n' "$v" "$("$PY" -c "print(10**$v)")" "$st"
    i=$((i+1))
  done
  echo; echo "queue: $(squeue -u "$USER" -h 2>/dev/null | grep -c bpemp) scan jobs running"
  exit 0
fi

# ------------------------------------------------------------------- report --
if [ "$MODE" = report ]; then
  GRID="$ALLGRID" PFX="$PFX" "$PY" - <<'PYREP'
import os, sys, glob
sys.path.insert(0, os.path.expanduser('~/CD16_NK92_project/filesCC'))
# clear any BOOT_* that would perturb the import-time config parsing
for k in list(os.environ):
    if k.startswith('BOOT_'): os.environ.pop(k)
import bootstrap_pzap as BP
import numpy as np, json

grid = os.environ['GRID'].split()
PFX = os.environ.get('PFX', 'zs')
KZP_TO_PHYS = 0.03           # kzp = KZP_MULT * 0.03  (uM s)^-1
rows = []
for i, v in enumerate(grid):
    root = os.path.expanduser('~/boot_pzap_%s%02d_pinZAP0' % (PFX, i))
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
# GUARD: bootstrap_pzap.build_one starts each point with
#     if os.path.isdir(dst): shutil.rmtree(dst)
# so resubmitting while jobs are live DELETES the working directory out from
# under them and they die.  That already happened once (8679684-89 were wiped
# by a rerun).  Refuse unless the user really means it.
nrun=$(squeue -u "$USER" -h -o '%j' 2>/dev/null | grep -c '^bpemp' || true)
if [ "${nrun:-0}" -gt 0 ] && [ "${FORCE:-0}" != "1" ]; then
  echo "  STOP: $nrun pZAP job(s) are still in the queue."
  echo "        Submitting now would delete their working directories and kill"
  echo "        them, losing the progress they have made."
  echo
  echo "        To look at them instead:"
  echo "           bash \$0 status      bash \$0 report      bash \$0 logs"
  echo
  echo "        If you really do want to scrap them and start over:"
  echo "           scancel -u $USER --name=bpemp000     # then resubmit"
  echo "        or force it in one step:"
  echo "           FORCE=1 bash \$0 submit"
  exit 1
fi
[ "${nrun:-0}" -gt 0 ] && echo "  warn FORCE=1 -- clobbering $nrun running job(s) on purpose"

SRCD=$HOME/NK92_fit_v77
XL=/home/gddaslab/share/Varun_Indrani/estimate_params_pzap/data/pZAP70_Tyr493_Tyr292_original_and_averages.xlsx
fail=0
[ -d "$SRCD" ]  && echo "  ok   source tree   $SRCD" || { echo "  FAIL source tree MISSING: $SRCD"; fail=1; }
[ -f "$XL" ]    && echo "  ok   day-level xlsx" || { echo "  FAIL xlsx MISSING: $XL"; fail=1; }
command -v sbatch >/dev/null && echo "  ok   sbatch on PATH" || { echo "  FAIL sbatch not found (are you on a login node?)"; fail=1; }

# openpyxl: pandas needs it to read the day-level .xlsx, and read_days() is the
# FIRST thing build() calls -- so without it all 8 points die before anything is
# built.  It was present for the earlier runs and went missing later (most
# likely collateral from `pip install flowkit` resolving its own deps), so this
# repairs it rather than just complaining.
if "$PY" -c "import openpyxl" >/dev/null 2>&1; then
  echo "  ok   openpyxl   $("$PY" -c 'import openpyxl;print(openpyxl.__version__)' 2>/dev/null)"
else
  echo "  --   openpyxl MISSING -- installing into $("$PY" -c 'import sys;print(sys.prefix)')"
  "$PY" -m pip install --quiet openpyxl 2>&1 | tail -3
  if "$PY" -c "import openpyxl" >/dev/null 2>&1; then
    echo "  ok   openpyxl now $("$PY" -c 'import openpyxl;print(openpyxl.__version__)')"
  else
    echo "  FAIL pip could not install openpyxl.  Try:"
    echo "         conda install -n CD16_v2 -c conda-forge openpyxl"
    fail=1
  fi
fi

# the rest of the env, so a broken resolver shows up here and not 8 hours in
"$PY" - <<'PYCHK' || fail=1
import importlib, sys
bad = []
for m in ('numpy', 'pandas', 'scipy', 'pyswarms'):
    try:
        v = getattr(importlib.import_module(m), '__version__', '?')
        print('  ok   %-9s %s' % (m, v))
    except Exception as e:
        print('  FAIL %-9s %s' % (m, e)); bad.append(m)
sys.exit(1 if bad else 0)
PYCHK

"$PY" -c "import sys; sys.path.insert(0,'$F'); import bootstrap_pzap" 2>&1 \
  && echo "  ok   bootstrap_pzap imports" || { echo "  FAIL bootstrap_pzap does not import (full error above)"; fail=1; }
avail=$(df -Pk "$HOME" 2>/dev/null | awk 'NR==2{print int($4/1048576)}')
[ -n "$avail" ] && echo "  note ${avail} GB free in \$HOME (each grid point copies the v77 tree)"
[ "$fail" = 0 ] || { echo; echo "PREFLIGHT FAILED -- nothing submitted."; exit 1; }
echo

i=0; nok=0; nbad=0
for v in $GRID; do
  t=$(tag_of "$v")
  z=$("$PY" -c "print(round(10**$v))")
  echo "--- ZAP0(log10)=$v  (ZAP0=$z)  tag=_${PFX}$t"
  log=$(mktemp)
  BOOT_TAG="_${PFX}$t" BOOT_PIN="ZAP0" BOOT_CENTRE="ZAP0=$v" \
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
