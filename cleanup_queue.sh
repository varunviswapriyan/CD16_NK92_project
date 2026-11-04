#!/usr/bin/env bash
# cleanup_queue.sh -- cancel the stray kin_* jobs and any blk_* jobs that are not
# the patched blockade arrays (8800621 male, 8800622 female). Keeps cd4v2_*.
# Run: cd ~/CD16_NK92_project && git pull --no-edit && bash cleanup_queue.sh
GOOD="8800621 8800622"
ids=$(squeue -u "$USER" -h -o "%A %F %j" | awk -v good="$GOOD" '
  BEGIN{n=split(good,g," "); for(i=1;i<=n;i++) keep[g[i]]=1}
  $3 ~ /^kin_/ || ($3 ~ /^blk_/ && !($2 in keep)) {print $1}' | sort -u)
if [ -n "$ids" ]; then echo "cancelling: $ids"; scancel $ids; else echo "nothing to cancel"; fi
sleep 3
echo "=== remaining queue (job : count) ==="
squeue -u "$USER" -h -o "%j" | sort | uniq -c
