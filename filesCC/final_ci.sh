#!/bin/bash
# =============================================================================
# final_ci.sh -- the two tables that go to Indrani, side by side.
#
#   bash ~/CD16_NK92_project/filesCC/final_ci.sh          # pZAP + Ca tables
#   bash ~/CD16_NK92_project/filesCC/final_ci.sh ca
#   bash ~/CD16_NK92_project/filesCC/final_ci.sh pzap
#   bash ~/CD16_NK92_project/filesCC/final_ci.sh band        # Ca fit + CI band, lin1
#   bash ~/CD16_NK92_project/filesCC/final_ci.sh band hillonly
#
# Reads finished results only.  Submits nothing, changes nothing.
#
# WHICH COLUMN TO QUOTE
#   pZAP prints two interval forms.  Quote the WEIGHTED PERCENTILE one.
#   Because all 10 possible day-resamples were enumerated and weighted by
#   their exact multinomial probabilities, the bootstrap distribution is known
#   exactly -- a discrete distribution on 10 atoms -- so its quantiles are
#   exact rather than estimated.  The "+/- 1.96 SD" column imposes a normal
#   shape on those 10 atoms and comes out conservatively wide; it is a bound,
#   not the answer.
# =============================================================================
set -o pipefail
WHICH=${1:-both}
VARIANT=${2:-lin1}
F=~/CD16_NK92_project/filesCC
CADIR=$HOME/Ca_fit_c02

case "$WHICH" in
  both|ca|pzap|band) ;;
  *) echo "ERROR: use 'both', 'ca', 'pzap' or 'band [variant]'"; exit 2 ;;
esac

# ---- interpreter (same probe as pairs.sh; $CONDA_PREFIX first) --------------
echo "setting up the job environment (module + conda)..."
[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh 2>/dev/null
module load Miniconda3/4.9.2 >/dev/null 2>&1
CONDA_SH=/gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh
if [ -f "$CONDA_SH" ]; then
  # shellcheck disable=SC1090
  . "$CONDA_SH" >/dev/null 2>&1 && conda activate CD16_v2 >/dev/null 2>&1
fi
PY="${PZ_PY:-}"
if [ -z "$PY" ]; then
  TMO=""; command -v timeout >/dev/null 2>&1 && TMO="timeout 180"
  for cand in "${CONDA_PREFIX:+$CONDA_PREFIX/bin/python}" python python3; do
    [ -n "$cand" ] || continue
    c=$(command -v "$cand" 2>/dev/null || { [ -x "$cand" ] && echo "$cand"; })
    [ -n "$c" ] || continue
    if $TMO "$c" -c "import pandas,numpy,matplotlib" >/dev/null 2>&1; then PY="$c"; break; fi
  done
fi
if [ -z "$PY" ]; then
  echo "ERROR: no python with pandas + matplotlib."
  echo "  module load Miniconda3/4.9.2 && source $CONDA_SH && conda activate CD16_v2"
  echo "then rerun, or:  PZ_PY=\$CONDA_PREFIX/bin/python bash \$0"
  exit 1
fi
echo "using $PY"

if [ "$WHICH" = band ]; then
  echo
  echo "##############################################################"
  echo "#  Calcium fit + bootstrap prediction band -- variant '$VARIANT'"
  echo "##############################################################"
  [ -d "$CADIR" ] || { echo "$CADIR not found"; exit 1; }
  cd "$CADIR" || exit 1

  # 'band' is the first command that reads the two input CSVs -- 'report' only
  # reads the JSON results.  Check them here so a directory or a missing file
  # gives a usable message instead of an IsADirectoryError traceback.
  bad=0
  for f in ca_data/Ca_NK92.csv optimized_model_pzap/model_output_pzap_v77.csv; do
    if [ -d "$CADIR/$f" ]; then
      echo "ERROR: $CADIR/$f is a DIRECTORY, not a file."
      echo "       Something created it with mkdir -p on the full path."
      echo "       Remove it and regenerate:  rmdir '$CADIR/$f'"
      bad=1
    elif [ ! -f "$CADIR/$f" ]; then
      echo "ERROR: missing $CADIR/$f"
      bad=1
    elif [ ! -s "$CADIR/$f" ]; then
      echo "ERROR: $CADIR/$f is empty"
      bad=1
    fi
  done
  if [ "$bad" -ne 0 ]; then
    echo
    echo "ca_data/Ca_NK92.csv               = the calcium measurements"
    echo "optimized_model_pzap/..._v77.csv  = the pZAP input curve the Ca model reads"
    echo "The second one is produced by the v77 pZAP run:"
    echo "    bash $F/run_ca_v77.sh"
    exit 1
  fi

  "$PY" "$F/bootstrap_ca2.py" band "$VARIANT" || echo "(band failed for '$VARIANT')"

  IP=$(hostname -I 2>/dev/null \
       | awk '{for(i=1;i<=NF;i++) if($i ~ /^10\.73\./) print $i}' | head -1)
  [ -n "$IP" ] || IP=$(hostname -I 2>/dev/null | awk '{print $1}')
  echo
  echo "=================== CALCIUM FIGURES ==================="
  # newest first, so the band just written is at the top
  find "$CADIR" -name '*.png' -newermt '-60 days' 2>/dev/null \
    | sed "s|^$HOME/|  http://$IP:8000/|" | sort | head -40
  echo
  echo "the band figure for '$VARIANT' is the one to show:"
  echo "  data points, the fit, and the envelope the bootstrap samples span"
  echo
  [ -n "$BOOT_NOSERVE" ] && exit 0
  echo "starting web server -- ctrl-C when done"
  cd "$HOME" && exec "$PY" -m http.server 8000
fi

if [ "$WHICH" = both ] || [ "$WHICH" = ca ]; then
  echo
  echo "##############################################################"
  echo "#  Calcium   --  'lin1' is the reportable variant"
  echo "#  (C1, C2, g, S fitted; k4 fixed at 1; all parameters shared"
  echo "#   across adaptors, only the ITAM count differs 6 / 2 / 4)"
  echo "##############################################################"
  if [ -d "$CADIR" ]; then
    cd "$CADIR" && "$PY" "$F/bootstrap_ca2.py" report
  else
    echo "$CADIR not found"
  fi
fi

if [ "$WHICH" = both ] || [ "$WHICH" = pzap ]; then
  echo
  echo "##############################################################"
  echo "#  pZAP70 (Tyr493)   --  quote the WEIGHTED PERCENTILE column"
  echo "##############################################################"
  if [ -x "$F/pzap_max.sh" ] || [ -f "$F/pzap_max.sh" ]; then
    BOOT_NOSERVE=1 bash "$F/pzap_max.sh" report
  else
    echo "pzap_max.sh not found in $F"
  fi
fi

echo
echo "=============================================================="
echo "Ca is printed FIRST on purpose: the pZAP report can end in a web"
echo "server, and anything after it would never be reached."
echo
echo "Reminder before this goes in an email:"
echo "  * on a pZAP row means the reference fit falls OUTSIDE its own"
echo "  interval -- bootstrap bias, needs a BCa correction first."
echo "=============================================================="
