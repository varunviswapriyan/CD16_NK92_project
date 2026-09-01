#!/bin/bash
#SBATCH --mail-user=indrani.nayak@nationwidechildrens.org
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=32
#SBATCH --job-name="Fit_Calcium_Params"
#SBATCH --output="slurm%j.out"
#SBATCH --time=2-00:00:00
# ml purge
# ml load miniforge3
# conda activate CD16_sig



srun python ca_param_estimation_BWZ_less_param.py 32 100

