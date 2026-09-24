#!/bin/bash
# go_lin1.sh -- verify the updated script is here, clear the duplicate jobs,
# and submit only the lin1 / lin1S variants.
#   bash ~/CD16_NK92_project/filesCC/go_lin1.sh
set -e
F=~/CD16_NK92_project/filesCC

if ! grep -q "'lin1'" $F/bootstrap_ca2.py; then
  echo "STOP: this copy of bootstrap_ca2.py has no lin1 variant."
  echo
  echo "On the Mac, the updated file is still in Downloads. Run:"
  echo "  mv ~/Downloads/bootstrap_ca2.py ~/Downloads/indrani_code/filesCC/"
  echo "  cd ~/Downloads/indrani_code && git add filesCC/ && git commit -m lin1 && git push"
  echo "Then here:"
  echo "  cd ~/CD16_NK92_project && git pull && bash \$0"
  exit 1
fi
echo "OK: lin1 variant present."

OLD=$(squeue -u "$USER" -h -o "%i %j" | awk '$2 ~ /^(full|k3fix|k4fix|hillo|linon|Sfix|min4)/ {print $1}')
if [ -n "$OLD" ]; then
  echo "$OLD" | xargs -r scancel
  echo "cancelled $(echo "$OLD" | grep -c .) duplicate jobs from the earlier variants"
else
  echo "no duplicate jobs to cancel"
fi

cd ~/Ca_fit_c02
python $F/bootstrap_ca2.py submit 12 lin1 lin1S
echo
echo "~30-45 min.  Then:  bash ~/CD16_NK92_project/filesCC/run_lin1.sh check"
