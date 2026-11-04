#!/usr/bin/env bash
# run_block_patched.sh -- task 4: true IFNGR blockade on ALL T cells.
#
# Receptors on AR- and AR+ T cells = 0 (cancer keeps 1000). Stock SPPARKS then
# computes 0/0 (bound/receptors, app_lattice.cpp ~line 895 and 959) -> NaN rates,
# which is why the first try "finished" instantly. This script builds a SEPARATE
# binary (spk_sexbias_block) from a copy of the source with those divisions
# guarded (x/0 -> 0, i.e. a blocked receptor gives no signal). The production
# binary and source are not touched.
#
# Run: cd ~/CD16_NK92_project && git pull --no-edit && bash run_block_patched.sh
set -uo pipefail
REPO=~/CD16_NK92_project/SPPARKS_Sex_Bias
cd "$REPO" || exit 1
SRC=rsc_out_2026-07-03
DECK="$PWD/$SRC/in.full"
B=~/spk_block_build
module load foss/2023a 2>/dev/null

echo "== 1. copy source -> $B"
rm -rf "$B"; mkdir -p "$B"; cp *.cpp *.h "$B"/; cd "$B"

echo "== 2. guard divisions by receptor counts"
perl -0pi -e '
  s{cell_bound\s*/\s*cell_receptors}{((cell_receptors>0)?(cell_bound/cell_receptors):0.0)}g;
  s{final_bound\s*/\s*total_receptors}{((total_receptors>0)?(final_bound/total_receptors):0.0)}g;
  s{(\b[A-Za-z_]\w*(?:\[[^\[\]]+\])+)\[0\]\s*/\s*\1\[1\]}{(($1\[1\]>0)?($1\[0\]/$1\[1\]):0.0)}g;
' app_lattice.cpp app_rxn_diff_custom.cpp app_rxn_diff.cpp
echo "   remaining unguarded divisions by a receptor-like quantity (review):"
grep -nE "/ *(cell_receptors|total_receptors|[A-Za-z_]+(\[[^]]+\])+\[1\])" app_lattice.cpp app_rxn_diff_custom.cpp app_rxn_diff.cpp \
  | grep -v "?" | head -20 | sed 's/^/     /'

echo "== 3. build"
mpicxx -O3 -std=c++14 -w -fpermissive *.cpp -o "$REPO/spk_sexbias_block" || { echo "BUILD FAILED"; exit 1; }
cd "$REPO"; ls -l spk_sexbias_block

BLOCK='/app_style/ s/ 2\.0 1000 1000 1000 / 2.0 0 0 1000 /'
echo "== 4. smoke tests (20 steps): unpatched-equivalent check + blocked"
for mode in std block; do
  d=smoke_block/$mode; rm -rf $d; mkdir -p $d/sex_bias_initialization_files
  cp $SRC/sex_bias_initialization_files/d14_M_4_SPPARKS_init $d/sex_bias_initialization_files/
  sed -e "s|read_sites.*|read_sites sex_bias_initialization_files/d14_M_4_SPPARKS_init|" -e "s|^run .*|run 20|" "$DECK" > $d/in.simulation
  [ $mode = block ] && sed -i "$BLOCK" $d/in.simulation
  ( cd $d && mpirun -np 1 "$REPO/spk_sexbias_block" < in.simulation > log.spparks 2>&1 )
  echo "   -- $mode (patched binary):"; grep -A3 "sum(i1)" $d/log.spparks | awk '{print "      ",$1,$6,$7,$8,$15}'
done
if grep -qi nan smoke_block/block/log.spparks; then
  echo "!! still NaN with the patch -- not submitting. Send Claude the lines above."; exit 1
fi
cmp -s <(grep -A3 "sum(i1)" smoke_block/std/log.spparks | awk '{print $7}') \
       <(grep -A3 "sum(i1)" runs_big/prod/male/s9001/log.spparks | awk '{print $7}') \
  && echo "   patched binary reproduces standard run start: yes" \
  || echo "   (std smoke start differs from prod s9001 -- expected if seed differs; check the numbers look sane)"

echo "== 5. submit 40 male + 40 female"
scancel -u "$USER" -n blk_m,blk_f 2>/dev/null; rm -rf runs_big/blockall
JOB=$PWD/runs_big/_job_generic.sh
[ -f "$JOB" ] || { echo "missing $JOB (run do_tasks.sh once first)"; exit 1; }
for sex in male female; do
  base=runs_big/blockall/$sex; mkdir -p $base
  [ $sex = male ] && { slide=d14_M_4_SPPARKS_init; fem=0; } || { slide=d14_F_4_SPPARKS_init; fem=1; }
  printf '%s\n' "$BLOCK" > $base/extra.sed
  cat > $base/cfg.sh <<EOF
BASE=$PWD/$base
SEEDBASE=9000
SLIDEPATH=$PWD/$SRC/sex_bias_initialization_files/$slide
SLIDENAME=$slide
FEM=$fem
DECK=$DECK
BIN=$PWD/spk_sexbias_block
EXTRASED=$PWD/$base/extra.sed
EOF
  sbatch --job-name=blk_${sex:0:1} --array=1-40 --time=24:00:00 --ntasks=1 \
         --output=$PWD/$base/slurm_%A_%a.out --export=ALL,CFG=$PWD/$base/cfg.sh "$JOB" \
    && echo "   submitted blockall $sex x40"
done
echo "When the queue is empty:  cd ~/CD16_NK92_project && bash do_tasks.sh plot"
