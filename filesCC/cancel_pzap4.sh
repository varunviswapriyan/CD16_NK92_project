#!/bin/bash
# cancel_pzap4.sh -- cancel the 4-parameter pZAP bootstrap and clear its dirs,
# so the 6-parameter version can be submitted cleanly.
#   bash ~/CD16_NK92_project/filesCC/cancel_pzap4.sh          # show what would go
#   bash ~/CD16_NK92_project/filesCC/cancel_pzap4.sh yes      # actually cancel
set -e
MODE="${1:-dry}"

IDS=$(squeue -u "$USER" -h -o "%i %j" | awk '$2 ~ /^bp(emp|par)/ {print $1}')
N=$(echo "$IDS" | grep -c . || true)

echo "pZAP bootstrap jobs found (name starts with bpemp/bppar): $N"
squeue -u "$USER" -h -o "%i %j %T %M" | awk '$2 ~ /^bp(emp|par)/ {printf "  %-10s %-10s %-9s %s\n", $1, $2, $3, $4}'

if [ "$MODE" != "yes" ]; then
  echo
  echo "Ca jobs (bca*, and the ca2 variants) are NOT touched."
  echo "To go ahead:  bash \$0 yes"
  exit 0
fi

if [ -n "$IDS" ]; then
  echo "$IDS" | xargs -r scancel
  echo "cancelled $N jobs"
else
  echo "nothing running to cancel"
fi

rm -rf ~/boot_pzap
echo "removed ~/boot_pzap"
echo
echo "Now submit the 6-parameter version:"
echo "  python ~/CD16_NK92_project/filesCC/bootstrap_pzap.py submit 12 12"
