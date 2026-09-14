#!/bin/bash
# run_filesCC.sh -- submit the full c10-c20 + c16 sweep. No smoke test, no options.
#
#   cd ~/CD16_NK92_project && git pull
#   bash ~/CD16_NK92_project/filesCC/run_filesCC.sh
#
# That's it. Two jobs get submitted:
#   ccA  shared family    (c13 c14 c15 c17 c18)      -> logs/ccA.log
#   ccB  per-cond + MAIN  (c10 c11 c19 c20 c12 c16)  -> logs/ccB.log
# ccB only starts after ccA finishes (-w), so results.json never collides:
# ccA's summary is saved aside as results_sharedfam.json before ccB begins.
#
# Next morning:
#   tail -40 ~/Ca_fit_c02/logs/ccA.log ~/Ca_fit_c02/logs/ccB.log
#   ls ~/Ca_fit_c02/out_filesCC/
#   cd ~/Ca_fit_c02/out_filesCC && python -m http.server 8000   # view plots

set -e

RUN_DIR="$HOME/Ca_fit_c02"                       # has ca_data/ + optimized_model_pzap/
SCRIPT="$HOME/CD16_NK92_project/filesCC/filesCC.py"
ENV_SETUP='module load Miniconda3/4.9.2
source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh
conda activate CD16_v2'

# ---- sanity checks so a bad submit fails NOW, not at 2am ----
[ -f "$SCRIPT" ] || { echo "ERROR: $SCRIPT not found (git pull first?)"; exit 1; }
[ -d "$RUN_DIR/ca_data" ] || { echo "ERROR: $RUN_DIR/ca_data missing"; exit 1; }
[ -d "$RUN_DIR/optimized_model_pzap" ] || { echo "ERROR: $RUN_DIR/optimized_model_pzap missing"; exit 1; }
command -v bsub >/dev/null || { echo "ERROR: bsub not found -- run this on the cluster login node"; exit 1; }

mkdir -p "$RUN_DIR/logs"

# ---- job A: shared-parameter family ----
bsub -J ccA -n 32 -R "span[hosts=1]" -W 24:00 \
     -oo "$RUN_DIR/logs/ccA.log" <<EOF
$ENV_SETUP
cd $RUN_DIR
python $SCRIPT c13 c14 c15 c17 c18
mv out_filesCC/results.json out_filesCC/results_sharedfam.json
EOF

# ---- job B: per-condition family + auto-MAIN c16, waits for A ----
bsub -J ccB -w "done(ccA)" -n 32 -R "span[hosts=1]" -W 48:00 \
     -oo "$RUN_DIR/logs/ccB.log" <<EOF
$ENV_SETUP
cd $RUN_DIR
python $SCRIPT c10 c11 c19 c20 c12 c16
EOF

echo ""
echo "Submitted. Check queue with:  bjobs"
echo "Watch progress with:          tail -f $RUN_DIR/logs/ccA.log"
echo "Results land in:              $RUN_DIR/out_filesCC/"
echo "  shared family summary  -> results_sharedfam.json + plots"
echo "  per-cond + c16 summary -> results.json + plots"
