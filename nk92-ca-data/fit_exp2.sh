#!/usr/bin/env bash
# fit_exp2.sh -- fit the Ca model to Experiment 2 (18 May) in its OWN directory,
#                leaving every Experiment 1 result untouched.
#
#   bash fit_exp2.sh            set up ~/Ca_fit_exp2 and submit 51 jobs
#   bash fit_exp2.sh setup      set it up, print the submit line, submit nothing
#   bash fit_exp2.sh status     squeue + how many samples have landed
#   bash fit_exp2.sh report     the CI table, once they finish
#
# Variant fitted: 'orig' -- Indrani's five parameters (C1 C2 g k3 k4), S fixed
# at 1, no S term. Same variant reported for Experiment 1, so the parameters are
# comparable between the two experiments.
#
# NOT comparable: the SSR. Exp 2 is trimmed at 221.6 s and Exp 1 runs to 331 s,
# so Exp 2 has ~145 time points per condition against Exp 1's ~228. A smaller
# SSR here means fewer residuals, not a better fit. Compare parameters and CIs.
set -uo pipefail

CMD="${1:-submit}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DST="${DST:-$HOME/Ca_fit_exp2}"
NSAMP="${NSAMP:-50}"
NEWCSV="${NEWCSV:-$HOME/ca_check/Ca_NK92_exp2_20260518.csv}"
F="${F:-$HOME/CD16_NK92_project/filesCC}"

say()  { printf '\n\033[1m== %s\033[0m\n' "$*"; }
ok()   { printf '   \033[32mOK\033[0m   %s\n' "$*"; }
warn() { printf '   \033[33mNOTE\033[0m %s\n' "$*"; }
bad()  { printf '   \033[31mFAIL\033[0m %s\n' "$*"; }

# ------------------------------------------------------------------ python env
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
[ -n "$PY" ] || { bad "no python with numpy/pandas"; exit 1; }
echo "using $PY"

BOOT="$F/bootstrap_ca2.py"
[ -f "$BOOT" ] || { bad "bootstrap_ca2.py not at $BOOT  (set F=<dir>)"; exit 1; }

# ---------------------------------------------------------------- status/report
if [ "$CMD" = status ] || [ "$CMD" = report ]; then
  [ -d "$DST" ] || { bad "$DST does not exist -- run 'bash fit_exp2.sh' first"; exit 1; }
  say "queue"
  squeue -u "$USER" 2>/dev/null | head -20
  n=$(ls "$DST"/out_boot_ca2/orig/boot_*.json 2>/dev/null | wc -l)
  echo "   samples written so far: $n  (expect $((NSAMP + 1)))"
  if [ "$CMD" = report ]; then
    say "report"
    ( cd "$DST" && "$PY" "$BOOT" report )
    echo
    warn "SSR is NOT comparable to Experiment 1 -- fewer time points. Compare parameters."
  fi
  exit 0
fi

# ------------------------------------------------------------------ band plots
# NOTE the subcommand is "band orig", with a space. "band_orig" is not a command
# and "org" is not a variant -- both just print the usage text.
if [ "$CMD" = band ]; then
  [ -d "$DST" ] || { bad "$DST does not exist"; exit 1; }
  for v in "${@:2}"; do
    ( cd "$DST" && "$PY" "$BOOT" band "$v" ) && \
      ok "http://\$(hostname -I | awk '{print \$1}'):8000/Ca_fit_exp2/out_boot_ca2/band_$v.png"
  done
  [ $# -gt 1 ] || ( cd "$DST" && "$PY" "$BOOT" band orig )
  exit 0
fi

# -------------------------------------------------------------------- variants
# 'orig' leaves BOTH the Hill threshold k3 and the linear coefficient k4 free.
# The script's own header says that combination makes a flat valley: k3, k4 and
# C1 trade off, so C1 can roam orders of magnitude at nearly constant SSR. That
# is what the Exp 2 CIs show. These three close the valley different ways:
#   orig1  her model, S=1, k4 pinned to 1  (C1 carries the linear gain)
#   lin1   linear drive only, k4=1, S free
#   min4   both terms, k3 and S pinned at the n01 values
if [ "$CMD" = variants ]; then
  [ -d "$DST" ] || { bad "$DST does not exist -- run 'bash fit_exp2.sh' first"; exit 1; }
  V="${VARIANTS:-orig1 lin1 min4}"
  N="${NSAMP2:-20}"
  say "submitting $(echo $V | wc -w) variants x $((N + 1)) samples on Exp 2"
  echo "   variants: $V"
  ( cd "$DST" && CA_WARM_FROM=exp1seed CA_STARTS=4 "$PY" "$BOOT" submit "$N" $V )
  exit $?
fi

# ----------------------------------------------------------------------- setup
say "setting up $DST"

[ -f "$NEWCSV" ] || { bad "no $NEWCSV -- run 'bash ca_go.sh' first"; exit 1; }
ok "exp 2 data: $NEWCSV"

# find the Exp 1 fit directory: the one holding a converged orig reference
SRC=""
while IFS= read -r f; do SRC="$(dirname "$(dirname "$(dirname "$f")")")"; break; done < <(
  find "$HOME" -maxdepth 6 -path '*/out_boot_ca2/orig/boot_000.json' 2>/dev/null | sort)
if [ -z "$SRC" ]; then
  bad "no converged Exp 1 reference (out_boot_ca2/orig/boot_000.json) found."
  bad "Without it the fit cold-starts, which is what got 'orig' stuck before."
  exit 1
fi
ok "exp 1 reference: $SRC/out_boot_ca2/orig/boot_000.json"

PZ="$(find "$SRC" -maxdepth 2 -name 'model_output_pzap_v77.csv' 2>/dev/null | head -1)"
[ -n "$PZ" ] || { bad "model_output_pzap_v77.csv not found under $SRC"; exit 1; }
ok "pZAP input: $PZ"

mkdir -p "$DST/ca_data" "$DST/out_boot_ca2/exp1seed" "$DST/logs_boot_ca2"
cp "$NEWCSV" "$DST/ca_data/Ca_NK92.csv"
mkdir -p "$DST/optimized_model_pzap"
cp "$PZ" "$DST/optimized_model_pzap/model_output_pzap_v77.csv"
# seed the warm start from Exp 1's converged optimum, under its own name so the
# 51 jobs cannot race the sample-0 job rewriting it
cp "$SRC/out_boot_ca2/orig/boot_000.json" "$DST/out_boot_ca2/exp1seed/boot_000.json"
ok "staged ca_data/, optimized_model_pzap/, out_boot_ca2/exp1seed/"

# sanity: can the data actually be loaded and costed?
( cd "$DST" && "$PY" - "$BOOT" <<'PYCHK'
import sys, runpy, numpy as np
mod = runpy.run_path(sys.argv[1], run_name='_probe')
base = mod['load_base']()
for k, d in base.items():
    print(f'   {k:8s} {len(d["t"]):4d} points, t = {d["t"].min():.1f} - {d["t"].max():.1f} s, '
          f'mean {d["mean"].min():.1f} - {d["mean"].max():.1f}')
PYCHK
) || { bad "could not load the Exp 2 data through bootstrap_ca2.py"; exit 1; }
ok "data loads; above is what the fit will see (t >= T0 = 30 s)"

SUBMIT="cd $DST && CA_WARM_FROM=exp1seed CA_STARTS=4 $PY $BOOT submit $NSAMP orig"

if [ "$CMD" = setup ]; then
  say "ready -- nothing submitted"
  echo "   $SUBMIT"
  exit 0
fi

# ---------------------------------------------------------------------- submit
nrun=$(squeue -u "$USER" -h 2>/dev/null | wc -l)
if [ "${nrun:-0}" -gt 0 ] && [ "${FORCE:-0}" != "1" ]; then
  warn "$nrun job(s) already queued."
  warn "This writes to $DST only, so it cannot disturb them."
  warn "Submit anyway with:  FORCE=1 bash $0"
  exit 1
fi

say "submitting $((NSAMP + 1)) jobs (sample 0 = reference, 1..$NSAMP = bootstrap)"
( cd "$DST" && CA_WARM_FROM=exp1seed CA_STARTS=4 "$PY" "$BOOT" submit "$NSAMP" orig )
rc=$?
echo
if [ "$rc" -ne 0 ]; then bad "submission failed (exit $rc) -- output above"; exit 1; fi
cat <<EOF

   Experiment 1 results are untouched; everything new is under $DST.

   bash $0 status
   bash $0 report
EOF
