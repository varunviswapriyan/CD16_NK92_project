#!/usr/bin/env bash
# rerun_timeouts.sh -- after the queue is EMPTY, resubmit any seed that did not
# reach Time>=1000 (e.g. killed by the 24 h limit) with a 48 h limit.
# Uses the same cfg.sh/job script as the original runs, so settings are identical.
# Run: cd ~/CD16_NK92_project && git pull --no-edit && bash rerun_timeouts.sh
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || exit 1
if squeue -u "$USER" -h -o "%j" | grep -qE "^(blk_|cd4v2_)"; then
  echo "blk_/cd4v2_ jobs are still in the queue -- wait until they finish, then rerun this."
  squeue -u "$USER" -h -o "%j %T" | sort | uniq -c; exit 0
fi
JOB=$PWD/runs_big/_job_generic.sh
for base in runs_big/blockall/male runs_big/blockall/female runs_big/cd4off_v2/male runs_big/cd4off_v2/female; do
  [ -f "$base/cfg.sh" ] || continue
  idx=""
  for i in $(seq 1 40); do
    f=$base/s$((9000+i))/log.spparks
    t=$(awk '/sum\(i1\)/{h=1;next} h&&NF>20&&$1+0==$1{t=$1} END{print t+0}' "$f" 2>/dev/null)
    awk -v t="${t:-0}" 'BEGIN{exit !(t>=1000)}' || idx="$idx${idx:+,}$i"
  done
  if [ -z "$idx" ]; then echo "  $base: all 40 finished"; continue; fi
  echo "  $base: resubmitting seeds (array index = seed-9000): $idx"
  for i in ${idx//,/ }; do mv "$base/s$((9000+i))" "$base/s$((9000+i))_timedout" 2>/dev/null; done
  name=re_$(echo $base | tr '/' '_' | sed 's/runs_big_//')
  sbatch --job-name="$name" --array="$idx" --time=48:00:00 --ntasks=1 \
         --output="$PWD/$base/slurm_re_%A_%a.out" --export=ALL,CFG="$PWD/$base/cfg.sh" "$JOB"
done
echo "When these finish:  bash do_tasks.sh plot && bash plot_growth_gstyle.sh"
