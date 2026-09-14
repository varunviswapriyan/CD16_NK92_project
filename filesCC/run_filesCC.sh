#!/bin/bash
# run_filesCC.sh -- submit the full c10-c20 + c16 sweep on Slurm. One command:
#
#   bash ~/CD16_NK92_project/filesCC/run_filesCC.sh
#
# Submits two jobs:
#   ccA  shared family    (c13 c14 c15 c17 c18)      -> logs/ccA.log
#   ccB  per-cond + MAIN  (c10 c11 c19 c20 c12 c16)  -> logs/ccB.log
# ccB starts only after ccA finishes OK (dependency), so results.json never
# collides: ccA's summary is renamed to results_sharedfam.json before ccB runs.
#
# Next morning:
#   tail -60 ~/Ca_fit_c02/logs/ccA.log ~/Ca_fit_c02/logs/ccB.log
#   cd ~/Ca_fit_c02/out_filesCC && python -m http.server 8000

set -e

RUN_DIR="$HOME/Ca_fit_c02"                       # has ca_data/ + optimized_model_pzap/
SCRIPT="$HOME/CD16_NK92_project/filesCC/filesCC.py"

# ---- sanity checks so a bad submit fails NOW, not at 2am ----
[ -f "$SCRIPT" ] || { echo "ERROR: $SCRIPT not found (git pull first?)"; exit 1; }
[ -d "$RUN_DIR/ca_data" ] || { echo "ERROR: $RUN_DIR/ca_data missing"; exit 1; }
[ -d "$RUN_DIR/optimized_model_pzap" ] || { echo "ERROR: $RUN_DIR/optimized_model_pzap missing"; exit 1; }
command -v sbatch >/dev/null || { echo "ERROR: sbatch not found"; exit 1; }

mkdir -p "$RUN_DIR/logs"

# ---- write job scripts ----
cat > "$RUN_DIR/logs/ccA.sh" <<'EOF'
#!/bin/bash
#SBATCH --job-name=ccA
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=24:00:00
module load Miniconda3/4.9.2
source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh
conda activate CD16_v2
cd ~/Ca_fit_c02
python ~/CD16_NK92_project/filesCC/filesCC.py c13 c14 c15 c17 c18
mv out_filesCC/results.json out_filesCC/results_sharedfam.json
EOF

cat > "$RUN_DIR/logs/ccB.sh" <<'EOF'
#!/bin/bash
#SBATCH --job-name=ccB
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=16G
#SBATCH --time=48:00:00
module load Miniconda3/4.9.2
source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh
conda activate CD16_v2
cd ~/Ca_fit_c02
python ~/CD16_NK92_project/filesCC/filesCC.py c10 c11 c19 c20 c12 c16
EOF

# ---- submit, B depends on A ----
JID_A=$(sbatch --parsable -o "$RUN_DIR/logs/ccA.log" "$RUN_DIR/logs/ccA.sh")
echo "Submitted ccA: job $JID_A"
JID_B=$(sbatch --parsable --dependency=afterok:$JID_A -o "$RUN_DIR/logs/ccB.log" "$RUN_DIR/logs/ccB.sh")
echo "Submitted ccB: job $JID_B (starts after ccA finishes)"

echo ""
echo "Check queue:     squeue -u \$USER"
echo "Watch progress:  tail -f $RUN_DIR/logs/ccA.log"
echo "Results:         $RUN_DIR/out_filesCC/  (plots + results_sharedfam.json + results.json)"
