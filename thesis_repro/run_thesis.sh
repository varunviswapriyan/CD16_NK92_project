#!/usr/bin/env bash
# run_thesis.sh -- reproduce Giuseppe's day-14 -> day-21 tumour-growth runs
# with the with-CD4s thesis model (in.full), for male / female / male-no-seq.
#
#   ./thesis_repro/run_thesis.sh                 # full: 5 seeds, run 1008
#   SMOKE=1 ./thesis_repro/run_thesis.sh         # quick: 1 seed, run 20 (build check)
#   SEEDS="8377 8378 8379" CONDS="male female" ./thesis_repro/run_thesis.sh
#
# Each full run is one MPI rank (KMC, sector no) and takes hours. Set JOBS to
# how many you want in parallel (one core each).
set -euo pipefail
cd "$(dirname "$0")/.."                      # repo root
KIT=thesis_repro
SRC=rsc_out_2026-07-03                        # where in.full + slides live
SPK=./spk_sexbias

[ -x "$SPK" ] || { echo "!! build first:  bash $KIT/build.sh"; exit 1; }

SEEDS="${SEEDS:-8377 8378 8379 8380 8381}"
CONDS="${CONDS:-male female male_noseq}"
JOBS="${JOBS:-4}"
STEPS="${STEPS:-1008}"                         # 1008 = day14->21 (1 week)
RES="${RES:-thesis_runs}"

if [ "${SMOKE:-0}" = "1" ]; then
  SEEDS="8377"; CONDS="male"; STEPS="20"; JOBS="1"
  echo "== SMOKE TEST: male, seed 8377, run $STEPS"
fi

# condition -> (slide, sex, extra make_input flags)
slide_for() { case "$1" in
    male)        echo "d14_M_4_SPPARKS_init male";;
    female)      echo "d14_M_4_madeF_SPPARKS_init female";;   # made-female (isolates AR)
    female_real) echo "d14_F_4_SPPARKS_init female";;         # real female slide
    male_noseq)  echo "d14_M_4_SPPARKS_init male";;           # + --no-seq added below
  esac; }

mkdir -p "$RES/logs"
jobs_file="$RES/jobs.txt"; : > "$jobs_file"
for cond in $CONDS; do
  read -r slide sex <<<"$(slide_for "$cond")"
  extra=""; [ "$cond" = "male_noseq" ] && extra="--no-seq"
  for seed in $SEEDS; do
    printf '%s\n' "$cond|$slide|$sex|$seed|$extra" >> "$jobs_file"
  done
done
echo "== $(wc -l < "$jobs_file") runs, $JOBS at a time, STEPS=$STEPS"

run_one() {
  IFS='|' read -r cond slide sex seed extra <<<"$1"
  local dir="$RES/${cond}_s${seed}"
  mkdir -p "$dir/sex_bias_initialization_files"
  cp -f "$SRC/sex_bias_initialization_files/$slide" "$dir/sex_bias_initialization_files/"
  python3 "$KIT/make_input.py" --infull "$SRC/in.full" --out "$dir/in.simulation" \
      --seed "$seed" --sex "$sex" --slide "$slide" $extra >/dev/null
  # temporarily force the requested number of steps
  sed -i.bak -E "s/^run[[:space:]].*/run             $STEPS/" "$dir/in.simulation"
  echo "   start $cond seed $seed"
  ( cd "$dir" && mpirun -np 1 ../../$SPK < in.simulation > log.spparks 2>&1 ) \
      && echo "   done  $cond seed $seed" \
      || echo "   FAILED $cond seed $seed -- see $dir/log.spparks"
}
export -f run_one; export RES KIT SRC SPK

# shellcheck disable=SC2016
cat "$jobs_file" | xargs -P "$JOBS" -I{} bash -c 'run_one "$@"' _ {}

echo
echo "== logs written. Next:  python3 $KIT/parse_and_plot.py $RES"
echo "   (for the smoke test, just look at $RES/male_s8377/log.spparks )"
