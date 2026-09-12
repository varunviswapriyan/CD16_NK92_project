#!/bin/bash
# run_ca_from_fit.sh <V> [n_reps]
#   1. export pZAP curves from ~/NK92_fit_v<V> (n_reps replicates, default 10)
#   2. make a working copy ~/Ca_fit_v<V> of Indrani's Ca folder with that pZAP as input
#   3. sbatch the Ca fit there
# Needs CD16_v2 env active (the export runs BioNetGen on the login node, ~5-10 min for 10 reps).
set -e
V=$1; N=${2:-10}
[ -n "$V" ] || { echo "usage: bash run_ca_from_fit.sh <version> [n_reps]"; exit 1; }
HOME_DIR=/home/gddaslab/vxv016
SRC=$(cd "$(dirname "$0")" && pwd)
FIT=$HOME_DIR/NK92_fit_v$V/estimate_params_pzap_cleaned_up
CA=$HOME_DIR/Ca_fit_v$V

cd $FIT
cp "$SRC/export_pzap.py" .
python export_pzap.py $N model_output_pzap_v$V.csv

rm -rf $CA
cp -r $HOME_DIR/Ca_fitting_NK92 $CA
rm -f $CA/Ca_model.csv $CA/slurm*.out
cp $FIT/model_output_pzap_v$V.csv $CA/optimized_model_pzap/model_output_pzap.csv
cp $HOME_DIR/Ca_fitting_NK92/optimized_model_pzap/model_output_pzap.csv $CA/optimized_model_pzap/model_output_pzap_INDRANI_ORIGINAL.csv 2>/dev/null || true

cat > $CA/run_ca.sh <<WRAP
#!/bin/bash
#SBATCH --mail-user=Varun.Viswapriyan@nationwidechildrens.org
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --job-name="Ca_v$V"
#SBATCH --output=$CA/slurm_ca_v${V}_%j.out
#SBATCH --time=4:00:00
#SBATCH --chdir=$CA

ml load Miniconda3/4.9.2
source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh
conda activate CD16_v2
python estimate_ca_params.py
WRAP
cd ~ && sbatch $CA/run_ca.sh
echo "Result will be in: $CA/slurm_ca_v${V}_<jobid>.out  and  $CA/Ca_model.csv"
