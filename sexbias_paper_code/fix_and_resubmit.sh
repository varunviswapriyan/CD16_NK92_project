#!/bin/bash
# Run ON THE CLUSTER after `git pull`. Re-runs the 12 runs (seeds 8378-8381)
# that never started because of the old run_all.sh bug. Seed 8377 is left
# alone (it is valid and may still be running in the first job).
set -e
cd "$(dirname "$0")"
for s in 8378 8379 8380 8381; do
  rm -f results/logs/male_s$s.log results/logs/female_s$s.log results/logs/male_noseq_s$s.log
done
SEEDS="8378 8379 8380 8381" sbatch run_sexbias.sh
echo "Submitted. Check with: squeue -u \$USER"
