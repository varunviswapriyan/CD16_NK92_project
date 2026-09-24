#!/bin/bash
# run_lin1.sh -- fit the identifiable linear-drive Ca variants and check them.
#
#   bash ~/CD16_NK92_project/filesCC/run_lin1.sh          # submit lin1 + lin1S
#   bash ~/CD16_NK92_project/filesCC/run_lin1.sh check    # ref SSR + CIs + bands + server
#
# lin1  : F = z (linear drive, k4 fixed at 1 so C1 carries the gain)
#         fits C1, C2, g, S  -- 4 shared params, same model as n01/k3fix
# lin1S : same with S fixed at 31.5 -- 3 shared params
#
# The decisive number is "ref SSR" (sample 0, un-jittered): it should come out
# at ~5.21e4, i.e. identical to n01.  If it does, the fit is unchanged and the
# CIs are the payoff.  Ca inputs are the FIXED v77 pZAP curves throughout.
set -e
F=~/CD16_NK92_project/filesCC
MODE="${1:-submit}"

cd ~/Ca_fit_c02
if [ "$MODE" = "submit" ]; then
  python $F/bootstrap_ca2.py submit 12 lin1 lin1S
  echo
  echo "~30-45 min per job, all parallel.  Then:  bash \$0 check"
  exit 0
fi

echo "=== how many samples have landed ==="
for V in lin1 lin1S; do
  N=$(ls out_boot_ca2/$V/boot_*.json 2>/dev/null | wc -l)
  echo "  $V: $N"
done

echo
echo "=== reference fit (sample 0, un-jittered) -- compare with n01 = 5.210e+04 ==="
for V in lin1 lin1S; do
  python3 -c "
import json, os
p='out_boot_ca2/$V/boot_000.json'
if os.path.exists(p):
    r=json.load(open(p))
    if r.get('ok'):
        print('  $V  SSR=%.4e  ' % r['ssr'] + '  '.join('%s=%.4g'%(k,v) for k,v in r['params'].items()))
    else:
        print('  $V  FAILED:', r.get('error','')[:60])
else:
    print('  $V  (not finished)')
"
done

echo
python $F/bootstrap_ca2.py report

echo
for V in lin1 lin1S; do python $F/bootstrap_ca2.py band $V || true; done

OUT=~/ci_plots; mkdir -p $OUT
cp out_boot_ca2/band_lin1*.png $OUT/ 2>/dev/null || true
cp out_filesD/plot_ca_n01_S.png $OUT/ca_n01_fit.png 2>/dev/null || true
IP=$(hostname -I | awk '{for(i=1;i<=NF;i++) if($i ~ /^10\.73\./) print $i}')
echo
ls -1 $OUT
echo
echo "Open:  http://${IP:-10.73.170.129}:8000/"
cd $OUT && python -m http.server 8000
