#!/bin/bash
# =============================================================================
# zap_kzp_map.sh -- Indrani point 2, the direct test, done so it cannot fail.
#
#   bash ~/CD16_NK92_project/filesCC/zap_kzp_map.sh submit
#   bash ~/CD16_NK92_project/filesCC/zap_kzp_map.sh status
#   bash ~/CD16_NK92_project/filesCC/zap_kzp_map.sh report
#
# THE QUESTION, VERBATIM
#   "A higher ZAP concentration may lead to a lower estimated ZAP binding rate
#    compared with the value used in our JI paper."
#
# WHY THIS DESIGN AND NOT THE PREVIOUS TWO
#   Both earlier attempts asked an OPTIMIZER to find kzp at each pinned ZAP0.
#   Attempt 1 boxed kzp to +/-0.6 log10, so it could not fall far enough to
#   show the trade-off even if it existed.  Attempt 2 freed five parameters at
#   once and came back under-converged -- it fit WORSE than attempt 1 from a
#   strictly larger search space, which is impossible at a true optimum.
#
#   So this does not ask an optimizer for kzp at all.  It PINS kzp as well as
#   ZAP0 and simply measures the SSR on a 2D grid.  There is no kzp search to
#   under-converge and no kzp bound to hit.  The remaining four parameters get
#   a small warm-started refit, which is the easy, well-conditioned case.
#
# HOW TO READ THE MAP
#   For each ZAP0 row, look at which kzp gives the lowest SSR.
#     best kzp FALLS as ZAP0 rises, and the row minima stay about equally low
#         -> Indrani is right.  ZAP0 and kzp trade off; only their product is
#            fixed by the data, so ZAP0 should be set from the literature and
#            kzp read off the valley.
#     best kzp stays put, and the row minima get worse as ZAP0 rises
#         -> no trade-off.  The data really does prefer ZAP0 ~ 140, and the gap
#            with the T cell literature is a real NK92 difference.
# =============================================================================
set -o pipefail
MODE=${1:-submit}
F=~/CD16_NK92_project/filesCC
BP=$F/bootstrap_pzap.py

# log10(ZAP0).  Stops at 3.70: at 4.10 and above NFsim cannot simulate the
# model at all (every parameter set returns a non-finite cost and pyswarms
# dies), so those points carry no information and are simply not run.
ZGRID=${ZGRID:-"2.14 2.45 2.75 3.05 3.35 3.70"}

# log10(KZP_MULT).  kzp = KZP_MULT * 0.03 (uM s)^-1.  0.3424 is the v77/JI
# value (kzp = 0.066).  The grid runs ~2 decades BELOW it, which is what a
# product-preserving trade-off would need over this ZAP0 range, and a little
# above.
KGRID=${KGRID:-"-1.50 -1.00 -0.50 0.00 0.3424 0.70 1.00"}

# Only four parameters are refit per point (lig0, kd10, SYK0, KPR_MULT), all
# warm-started, so a small budget is plenty and the fit is well conditioned.
export BOOT_PARTICLES=${BOOT_PARTICLES:-16}
export BOOT_ITERS=${BOOT_ITERS:-15}
export BOOT_NREPS=${BOOT_NREPS:-6}
export BOOT_WALLTIME=${BOOT_WALLTIME:-4:00:00}
export BOOT_HALF=${BOOT_HALF:-0.60}
export BOOT_TMAX=${BOOT_TMAX:-300}
export BOOT_MEM=${BOOT_MEM:-64G}

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

root_of() {  # $1 = z index, $2 = k index
  printf '%s/boot_pzap_km%02d_%02d_pinZAP0_KZP_MULT' "$HOME" "$1" "$2"
}

# ------------------------------------------------------------------- status --
if [ "$MODE" = status ]; then
  zi=0; done=0; tot=0
  for z in $ZGRID; do
    ki=0; line=""
    for k in $KGRID; do
      f="$(root_of $zi $ki)/emp000/estimate_params_pzap_cleaned_up/analysis_param_residue.dat"
      if [ -f "$f" ] && grep -qi '^linear' "$f" 2>/dev/null; then line="$line X"; done=$((done+1))
      else line="$line ."; fi
      tot=$((tot+1)); ki=$((ki+1))
    done
    printf '  ZAP0=%-7.0f %s\n' "$("$PY" -c "print(10**$z)")" "$line"
    zi=$((zi+1))
  done
  echo
  echo "  $done of $tot grid points finished   (X = done, . = pending)"
  echo "  queue: $(squeue -u "$USER" -h 2>/dev/null | grep -c bpemp) job(s)"
  exit 0
fi

# ------------------------------------------------------------------- report --
if [ "$MODE" = report ]; then
  ZGRID="$ZGRID" KGRID="$KGRID" "$PY" - <<'PYREP'
import os, sys
sys.path.insert(0, os.path.expanduser('~/CD16_NK92_project/filesCC'))
for k in list(os.environ):
    if k.startswith('BOOT_'): os.environ.pop(k)
import bootstrap_pzap as BP
import numpy as np

Z = [float(x) for x in os.environ['ZGRID'].split()]
K = [float(x) for x in os.environ['KGRID'].split()]
KZP = 0.03

S = np.full((len(Z), len(K)), np.nan)
for i in range(len(Z)):
    for j in range(len(K)):
        d = os.path.expanduser('~/boot_pzap_km%02d_%02d_pinZAP0_KZP_MULT/emp000' % (i, j))
        r = BP.parse_result(d)
        if r and r[0] is not None:
            S[i, j] = float(r[0])

nd = int(np.isfinite(S).sum())
print('=' * 78)
print('POINT 2: SSR over the (ZAP0, kzp) grid -- both PINNED, nothing searched')
print('=' * 78)
print('%d of %d grid points finished\n' % (nd, S.size))

hdr = '%9s |' % 'ZAP0\\kzp'
for k in K:
    hdr += ' %9.4g' % (10**k * KZP)
print(hdr); print('-' * len(hdr))
for i, z in enumerate(Z):
    line = '%9.0f |' % (10**z)
    for j in range(len(K)):
        line += ' %9.4g' % S[i, j] if np.isfinite(S[i, j]) else ' %9s' % '-'
    print(line)
print('-' * len(hdr))

if nd < 4:
    print('\nnot enough points yet -- run `status`, wait, re-report.')
    raise SystemExit(0)

print('\nbest kzp at each ZAP0 (the row minimum):')
print('  %10s %12s %12s' % ('ZAP0', 'best kzp', 'its SSR'))
bz, bk, bs = [], [], []
for i, z in enumerate(Z):
    row = S[i]
    if not np.isfinite(row).any():
        continue
    j = int(np.nanargmin(row))
    # a minimum sitting on the edge of the kzp grid means the real best is
    # outside it, so the trend cannot be read from this row
    edge = ' (at grid edge -- widen KGRID)' if j in (0, len(K) - 1) else ''
    print('  %10.0f %12.4g %12.4g%s' % (10**z, 10**K[j] * KZP, row[j], edge))
    bz.append(10**z); bk.append(10**K[j] * KZP); bs.append(row[j])

print('\n' + '=' * 78)
if len(bz) >= 3:
    bz, bk, bs = np.array(bz), np.array(bk), np.array(bs)
    slope = np.polyfit(np.log(bz), np.log(bk), 1)[0]
    rng = bs.max() / bs.min() if bs.min() > 0 else np.inf
    prod = bz * bk
    cv = prod.std() / prod.mean() if prod.mean() else np.inf
    print('slope d[log best-kzp]/d[log ZAP0] : %+.2f   (-1.00 = exact trade-off)' % slope)
    print('row-minimum SSR  worst/best       : %.2fx' % rng)
    print('ZAP0 * best-kzp  spread (CV)      : %.0f%%' % (100 * cv))
    print()
    if slope < -0.5 and rng < 2.0:
        print('VERDICT: the best kzp FALLS as ZAP0 rises, and the fit stays about')
        print('as good along that line.  Indrani is right -- ZAP0 and kzp trade')
        print('off, the data fixes their product and not each one.  Set ZAP0 from')
        print('the literature and quote the kzp this map gives at that ZAP0.')
    elif rng >= 2.0 and slope > -0.5:
        print('VERDICT: the best kzp does NOT fall with ZAP0, and the fit gets')
        print('%.1fx worse as ZAP0 rises.  There is no trade-off; the data' % rng)
        print('genuinely prefers ZAP0 ~ 140.  The gap with the T cell literature')
        print('is a real NK92 difference, not a fitting artefact.')
    else:
        print('VERDICT: mixed (slope %+.2f, SSR range %.2fx).  Paste this table' % (slope, rng))
        print('and we read it together rather than trusting the rule of thumb.')
print('=' * 78)
PYREP
  exit 0
fi

# ------------------------------------------------------------------- submit --
nrun=$(squeue -u "$USER" -h -o '%j' 2>/dev/null | grep -c '^bpemp' || true)
if [ "${nrun:-0}" -gt 0 ] && [ "${FORCE:-0}" != "1" ]; then
  echo "STOP: $nrun pZAP job(s) still queued.  Submitting now would delete their"
  echo "      working directories and kill them.  Use  bash \$0 status  instead,"
  echo "      or  FORCE=1 bash \$0 submit  to scrap them deliberately."
  exit 1
fi

echo "--- preflight ---"
fail=0
SRCD=$HOME/NK92_fit_v77
XL=/home/gddaslab/share/Varun_Indrani/estimate_params_pzap/data/pZAP70_Tyr493_Tyr292_original_and_averages.xlsx
[ -d "$SRCD" ] && echo "  ok   source tree" || { echo "  FAIL missing $SRCD"; fail=1; }
[ -f "$XL" ]   && echo "  ok   day-level xlsx" || { echo "  FAIL missing $XL"; fail=1; }
command -v sbatch >/dev/null && echo "  ok   sbatch" || { echo "  FAIL no sbatch"; fail=1; }
if "$PY" -c "import openpyxl" >/dev/null 2>&1; then
  echo "  ok   openpyxl"
else
  echo "  --   installing openpyxl"
  "$PY" -m pip install --quiet openpyxl 2>&1 | tail -2
  "$PY" -c "import openpyxl" >/dev/null 2>&1 && echo "  ok   openpyxl installed" \
    || { echo "  FAIL openpyxl; try: conda install -n CD16_v2 -c conda-forge openpyxl"; fail=1; }
fi
"$PY" -c "import numpy,pandas,scipy,pyswarms" 2>&1 && echo "  ok   numpy/pandas/scipy/pyswarms" \
  || { echo "  FAIL python deps (above)"; fail=1; }
[ "$fail" = 0 ] || { echo; echo "PREFLIGHT FAILED -- nothing submitted."; exit 1; }

nz=$(echo $ZGRID | wc -w); nk=$(echo $KGRID | wc -w)
echo
echo "=============== (ZAP0, kzp) MAP ==============="
echo "  ZAP0 points : $nz      kzp points : $nk      total jobs : $((nz*nk))"
echo "  both ZAP0 and KZP_MULT are PINNED at each point; lig0, kd10, SYK0 and"
echo "  KPR_MULT get a warm-started ${BOOT_PARTICLES}x${BOOT_ITERS} refit"
echo

zi=0; nok=0; nbad=0
for z in $ZGRID; do
  ki=0
  for k in $KGRID; do
    log=$(mktemp)
    BOOT_TAG="_km$(printf '%02d_%02d' $zi $ki)" \
    BOOT_PIN="ZAP0,KZP_MULT" \
    BOOT_CENTRE="ZAP0=$z,KZP_MULT=$k" \
      "$PY" "$BP" submit 0 0 0 >"$log" 2>&1
    rc=$?
    jid=$(grep -o 'Submitted batch job [0-9]*' "$log" | head -1 | awk '{print $NF}')
    if [ "$rc" -ne 0 ] || [ -z "$jid" ]; then
      nbad=$((nbad+1))
      echo "  ZAP0=$z kzp_mult=$k  FAILED (exit $rc):"
      sed 's/^/      /' "$log"
    else
      nok=$((nok+1))
      printf '  ZAP0=%-7s kzp_mult=%-8s job %s\n' "$z" "$k" "$jid"
    fi
    rm -f "$log"
    ki=$((ki+1))
  done
  zi=$((zi+1))
done

echo
echo "=================================================="
echo "  submitted OK : $nok of $((nok+nbad))"
[ "$nbad" -gt 0 ] && echo "  FAILED       : $nbad  (errors printed above)"
if [ "$nok" -gt 0 ]; then
  echo "  later:  bash \$0 status      then      bash \$0 report"
else
  echo "  NOTHING submitted -- read the error above."
fi
echo "=================================================="
