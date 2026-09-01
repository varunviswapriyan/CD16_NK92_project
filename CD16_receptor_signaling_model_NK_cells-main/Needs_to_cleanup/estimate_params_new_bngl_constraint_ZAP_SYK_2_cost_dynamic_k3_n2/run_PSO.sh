#!/bin/bash
#SBATCH --mail-user=indrani.nayak@nationwidechildrens.org
#SBATCH --mail-type=FAIL
#SBATCH --cpus-per-task=32
#SBATCH --job-name="n_2_ZAP_SYK_Fit_Calcium_Params"
#SBATCH --output="slurm%j.out"
#SBATCH --time=4-00:00:00
#ml purge
# ml load miniforge3
# conda activate CD16_sig

srun BWZ_zeta_only.py 32 22
#srun python param_estimation_BWZ_zeta_hetero_k3.py 32 90
#srun python ca_param_estimation_BWZ_less_param_constraint.py 32 90
#srun python ca_param_estimation_BWZ_gamma.py 5 3
