#!/bin/bash
# Baseline: Indrani's Ca fit with HER pZAP input, unchanged. Submit from ~/Ca_fitting_NK92:
#   cd ~/Ca_fitting_NK92 && sbatch ~/CD16_NK92_project/ca_batch/run_ca.sh
#SBATCH --mail-user=Varun.Viswapriyan@nationwidechildrens.org
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --job-name="Ca_base"
#SBATCH --output=/home/gddaslab/vxv016/Ca_fitting_NK92/slurm_ca_base_%j.out
#SBATCH --time=4:00:00
#SBATCH --chdir=/home/gddaslab/vxv016/Ca_fitting_NK92

ml load Miniconda3/4.9.2
source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh
conda activate CD16_v2
python estimate_ca_params.py
