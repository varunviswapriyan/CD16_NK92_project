#!/bin/bash
# =============================================================================
# pairs.sh -- the FAST look: redundancy scatter plots from samples already on
#             disk.  Submits nothing, changes nothing, takes seconds.
#
#   bash ~/CD16_NK92_project/filesCC/pairs.sh          # every arm with data, then serve
#   bash ~/CD16_NK92_project/filesCC/pairs.sh a        # arm A only
#   bash ~/CD16_NK92_project/filesCC/pairs.sh b        # arm B only (currently the
#                                                      #   only complete arm)
#   BOOT_NOSERVE=1 bash ... pairs.sh                   # make the plots, no web server
#
# WHAT YOU ARE LOOKING AT
#   ci_pairs_*.png -- bootstrap samples plotted two parameters at a time.
#     points on a tilted LINE  -> only the COMBINATION is constrained.  That is
#                                 redundancy: fix one, fit the other.
#     a diffuse BLOB           -> both do real work, the data is just thin.
#                                 Keep both; the wide CI is a data problem.
#     red star                 -> the reference (un-resampled) fit.
#
#   The CI numbers this prints are UNWEIGHTED percentiles and will not match
#   pzap_max.sh, which weights each day-set by its multinomial probability.
#   Use this for the shapes; use pzap_max.sh report for numbers you quote.
# =============================================================================
set -o pipefail
WHICH=${1:-all}
F=~/CD16_NK92_project/filesCC
SUB=estimate_params_pzap_cleaned_up

# ---- find a python that actually has pandas + matplotlib --------------------
# The login shell is non-interactive here, so `module` may be undefined and the
# system python3 usually has no matplotlib.  Try the job environment first,
# then fall back, and say plainly what to do if neither works.
[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh 2>/dev/null
module load Miniconda3/4.9.2 >/dev/null 2>&1
CONDA_SH=/gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh
if [ -f "$CONDA_SH" ]; then
  # shellcheck disable=SC1090
  . "$CONDA_SH" >/dev/null 2>&1 && conda activate CD16_v2 >/dev/null 2>&1
fi

PY=""
for cand in python python3 /gpfs0/scratch/miniforge3/24.11.2/envs/CD16_v2/bin/python; do
  c=$(command -v "$cand" 2>/dev/null || { [ -x "$cand" ] && echo "$cand"; })
  [ -n "$c" ] || continue
  if "$c" -c "import pandas,numpy,matplotlib" >/dev/null 2>&1; then PY="$c"; break; fi
done
if [ -z "$PY" ]; then
  echo "ERROR: could not find a python with pandas + matplotlib."
  echo "Run these three lines first, then rerun this script:"
  echo "  module load Miniconda3/4.9.2"
  echo "  source $CONDA_SH"
  echo "  conda activate CD16_v2"
  exit 1
fi
echo "using $PY"
echo

# ---- arms ------------------------------------------------------------------
# name | directory | environment that makes bootstrap_pzap.py point at it
arms_for () {
  case "$1" in
    a)   echo "A|$HOME/boot_pzap|" ;;
    b)   echo "B|$HOME/boot_pzap_t600|BOOT_TMAX=600" ;;
    a2)  echo "A2|$HOME/boot_pzap_n16|BOOT_TAG=_n16" ;;
    p50) echo "P50|$HOME/boot_pzap_p50|BOOT_TAG=_p50" ;;
    all) arms_for a; arms_for b; arms_for a2; arms_for p50 ;;
  esac
}

n_done () {
  local root=$1 d dat n=0
  for d in "$root"/emp[0-9][0-9][0-9] "$root"/par[0-9][0-9][0-9]; do
    [ -d "$d" ] || continue
    dat="$d/$SUB/analysis_param_residue.dat"
    [ -f "$dat" ] && grep -qi '^linear' "$dat" 2>/dev/null && n=$((n+1))
  done
  echo "$n"
}

case "$WHICH" in
  a|b|a2|p50|all) ;;
  *) echo "ERROR: unknown arm '$WHICH' -- use: a, b, a2, p50, all"; exit 2 ;;
esac

HOST=$(hostname -I 2>/dev/null | awk '{print $1}')
[ -n "$HOST" ] || HOST=$(hostname)
LINKS=""
ANY=0

while IFS='|' read -r name root envset; do
  [ -n "$name" ] || continue
  if [ ! -d "$root" ]; then
    printf '%-4s %-28s not started\n' "$name" "$(basename "$root")"; continue
  fi
  n=$(n_done "$root")
  if [ "$n" -lt 3 ]; then
    printf '%-4s %-28s only %s finished fits -- too few to plot\n' \
           "$name" "$(basename "$root")" "$n"
    continue
  fi
  echo "================ arm $name  ($n finished fits) ================"
  # shellcheck disable=SC2086
  if env $envset "$PY" "$F/bootstrap_pzap.py" report; then
    ANY=1
    b=$(basename "$root")
    for img in ci_pairs_pzap.png ci_hist_pzap.png; do
      [ -f "$root/$img" ] && LINKS="$LINKS
  http://$HOST:8000/$b/$img"
    done
  else
    echo "  (report failed for arm $name)"
  fi
  echo
done <<EOF
$(arms_for "$WHICH")
EOF

if [ "$ANY" -eq 0 ]; then
  echo "No arm had enough finished fits to plot yet."
  exit 0
fi

echo "=========================== IMAGES ==========================="
echo "$LINKS"
echo
echo "The one to look at is ci_pairs_pzap.png:"
echo "  tilted line = redundant pair (fix one, fit the other)"
echo "  round blob  = both real, data just thin"
echo

[ -n "$BOOT_NOSERVE" ] && exit 0
echo "starting web server -- ctrl-C when done"
cd "$HOME" && "$PY" -m http.server 8000
