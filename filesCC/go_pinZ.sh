#!/bin/bash
# go_pinZ.sh -- pZAP bootstrap with ZAP0 pinned. 25 jobs, one arm.
#
# The pZAP cost function divides every curve by g_max (or its own max), so it
# is scale-free: parameters that move amplitude together (lig0, ZAP0,
# KZP_MULT) slide freely and get wide CIs no matter how many replicates run.
# Pinning ZAP0 at its v77 value removes that sliding direction -- same model,
# same rules, same optimizer, one fewer dimension searched.
#
#   bash ~/CD16_NK92_project/filesCC/go_pinZ.sh          # submit 25 jobs
#   bash ~/CD16_NK92_project/filesCC/go_pinZ.sh report   # table + plots + server
set -e
F=~/CD16_NK92_project/filesCC
export BOOT_PIN=ZAP0

if [ "$1" = "report" ]; then
  python $F/bootstrap_pzap.py report
  OUT=~/ci_plots; mkdir -p $OUT
  cp ~/boot_pzap_pinZAP0/ci_hist_pzap.png  $OUT/pzap_pinZ_hist.png  2>/dev/null || true
  cp ~/boot_pzap_pinZAP0/ci_pairs_pzap.png $OUT/pzap_pinZ_pairs.png 2>/dev/null || true
  IP=$(hostname -I | awk '{for(i=1;i<=NF;i++) if($i ~ /^10\.73\./) print $i}')
  echo; ls -1 $OUT
  echo; echo "Open:  http://${IP:-10.73.170.129}:8000/"
  cd $OUT && python -m http.server 8000
  exit 0
fi

OLD=$(squeue -u "$USER" -h -o "%i %j" | awk '$2 ~ /^bp/ {print $1}')
[ -n "$OLD" ] && { echo "$OLD" | xargs -r scancel; echo "cancelled $(echo "$OLD"|grep -c .) old pZAP jobs"; }
rm -rf ~/boot_pzap_pinZAP0

python $F/bootstrap_pzap.py submit 12 12
echo
echo "25 jobs (12 emp + 12 par + 1 reference), 4-6 h each in parallel."
echo "When done:  bash \$0 report"
echo
echo "Read the reference SSR first: if it matches the free-arm run, pinning"
echo "cost nothing in fit quality and the tighter CIs are real."
