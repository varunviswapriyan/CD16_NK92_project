#!/bin/bash
#SBATCH --mail-user=Varun.Viswapriyan@nationwidechildrens.org
#SBATCH --mail-type=END,FAIL
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --job-name=sexbias
#SBATCH --output=/home/gddaslab/vxv016/CD16_NK92_project/sexbias_paper_code/slurm_sexbias_%j.out
#SBATCH --time=1-00:00:00
#SBATCH --chdir=/home/gddaslab/vxv016/CD16_NK92_project/sexbias_paper_code

module load GCC/13.2.0
source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh
conda activate CD16_v2
JOBS=15 ./run_all.sh
