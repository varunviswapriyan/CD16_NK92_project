#!/usr/bin/env bash
# ============================================================================
# run_all.sh — build, run every simulation needed for the manuscript figures,
# and draw the figures.
#
#   ./run_all.sh                  # defaults below
#   SEEDS="8377 8378 8379" JOBS=3 ./run_all.sh
#   FEMALE_INIT=init/d14_F_4_SPPARKS_init ./run_all.sh   # real female slide
#
# Each 7-day run takes ~2–3 h of CPU on a laptop core (~10 s per 10-min unit).
# JOBS runs in parallel; set it to your number of physical cores.
# ============================================================================
set -euo pipefail
cd "$(dirname "$0")"

SEEDS="${SEEDS:-8377 8378 8379 8380 8381}"
JOBS="${JOBS:-4}"
DAYS="${DAYS:-7}"                      # day 14 -> day 21
SNAP_DAYS="${SNAP_DAYS:-0,3.5,7}"      # day 14, 17.5, 21 (first seed only)
MALE_INIT="${MALE_INIT:-init/d14_M_4_SPPARKS_init}"
FEMALE_INIT="${FEMALE_INIT:-init/d14_M_4_madeF_SPPARKS_init}"
RES="${RES:-results}"

CXX="${CXX:-g++}"
echo "== building"
$CXX -O3 -std=c++17 sex_bias_sim_v7.cpp -o sex_bias_sim_v7
$CXX -O3 -std=c++17 ifng_effective_diffusion.cpp -o ifng_effective_diffusion

mkdir -p "$RES/runs" "$RES/diffusion" "$RES/logs"

echo "== effective-diffusion model (~15 s)"
( cd "$RES/diffusion" && ../../ifng_effective_diffusion ifng_diffusion 7 0.01 )

# job list: condition init sex extra-flags
FIRST_SEED="${SEEDS%% *}"
jobs_file="$RES/jobs.txt"; : > "$jobs_file"
for seed in $SEEDS; do
  snap=""; [ "$seed" = "$FIRST_SEED" ] && snap="--snapshot-days $SNAP_DAYS"
  echo "male       $MALE_INIT   male   $seed $snap"                     >> "$jobs_file"
  echo "female     $FEMALE_INIT female $seed $snap"                     >> "$jobs_file"
  echo "male_noseq $MALE_INIT   male   $seed --no-sequestration"         >> "$jobs_file"
done

echo "== $(wc -l < "$jobs_file") simulation runs, $JOBS at a time"
run_one() {
  local cond=$1 init=$2 sex=$3 seed=$4; shift 4
  local out="$RES/runs/${cond}_s${seed}"
  if [ -f "${out}_timecourse.csv" ] && awk -F, -v d="$DAYS" 'END{exit !($1+0 >= d*144-0.01)}' "${out}_timecourse.csv"; then
    echo "   skip $cond seed $seed (done)"; return 0
  fi
  echo "   start $cond seed $seed"
  ./sex_bias_sim_v7 --init "$init" --sex "$sex" --days "$DAYS" --seed "$seed" \
      --report-hours 6 --out "$out" "$@" > "$RES/logs/${cond}_s${seed}.log" 2>&1
  echo "   done  $cond seed $seed"
}
export -f run_one; export RES DAYS
xargs -P "$JOBS" -L 1 bash -c 'run_one "$@"' _ < "$jobs_file"

echo "== figures"
python3 make_figures.py "$RES"
echo "All done: figures in $RES/figures/"
