#!/usr/bin/env bash
# run_diffusion.sh -- effective-diffusion figure (Fig4 #5) in the NEW model.
# Compiles Giuseppe's standalone single-source diffusion program, runs it
# (bound + bare calibration), and plots front^2 vs time with the fit -> D_eff.
# FAST (minutes) -- NOT a Slurm job. Run on a login node.
#   cd ~/CD16_NK92_project && git pull && bash run_diffusion.sh
set -uo pipefail
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || cd ~/CD16_NK92_project || { echo "repo not found"; exit 1; }
CPP=ifng_effective_diffusion.cpp
[ -f "$CPP" ] || CPP=~/CD16_NK92_project/ifng_effective_diffusion.cpp
[ -f "$CPP" ] || { echo "ifng_effective_diffusion.cpp not found (add it to the repo)"; exit 1; }

module load foss/2023a 2>/dev/null || true
echo ">> compiling..."
g++ -O3 -std=c++17 "$CPP" -o ifng_eff_diff || { echo "compile failed"; exit 1; }

echo ">> running (bound + bare, 5 days)..."
./ifng_eff_diff diffusion 5 0.01    # prefix=diffusion, 5 days, front fraction 0.01

source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh 2>/dev/null && conda activate CD16_v2 2>/dev/null
python3 ~/CD16_NK92_project/thesis_repro/plot_diffusion.py diffusion_front.csv
echo ">> done. figure: Fig_effective_diffusion.png"
