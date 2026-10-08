#!/bin/bash
#SBATCH --mail-user=Varun.Viswapriyan@nationwidechildrens.org
#SBATCH --mail-type=END,FAIL
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --job-name=sexbias_rf
#SBATCH --output=/home/gddaslab/vxv016/CD16_NK92_project/sexbias_paper_code/slurm_sexbias_rf_%j.out
#SBATCH --time=1-00:00:00
#SBATCH --chdir=/home/gddaslab/vxv016/CD16_NK92_project/sexbias_paper_code

# Male vs REAL female slide (d14_F_4), as a side comparison.
# Writes to results_realfemale/ and NEVER touches results/.
# The male and male_noseq runs don't depend on the female slide, so they are
# copied from results/ and skipped — only the 5 real-female runs execute (~1 h).

module load GCC/13.2.0
source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh
conda activate CD16_v2

RES=results_realfemale
mkdir -p "$RES/runs" "$RES/logs" "$RES/diffusion"
if [ -d results/runs ]; then
  cp -n results/runs/male_s*_timecourse.csv  results/runs/male_s*_snap_*.txt \
        results/runs/male_noseq_s*_timecourse.csv "$RES/runs/" 2>/dev/null || true
fi

FEMALE_INIT=init/d14_F_4_SPPARKS_init RES="$RES" JOBS=15 ./run_all.sh
