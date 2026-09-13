#!/bin/bash
# results_ca.sh            -> SSR + params for every Ca_fit_c* folder
# results_ca.sh plot c02 c07  -> plots -> ~/plot_ca_cNN.png
HOME_DIR=/home/gddaslab/vxv016
if [ "$1" == "plot" ]; then
  shift; for ID in "$@"; do
    cd $HOME_DIR/Ca_fit_$ID || continue
    python check_ca.py $HOME_DIR/plot_ca_$ID.png 2>&1 | grep -vE "^Saved" ; echo "    -> ~/plot_ca_$ID.png"
  done; exit 0
fi
echo "Indrani baseline (her pZAP, her bounds, L-BFGS): total SSR 8.9066e+04   zeta 1.727e4  gamma 4.369e4  hetero 2.811e4   C1=10 (on bound)"
echo "v58 pZAP, her bounds, L-BFGS:                    total SSR 8.4195e+04   zeta 1.462e4  gamma 4.592e4  hetero 2.366e4   C1=10 (on bound)"
echo
for D in $HOME_DIR/Ca_fit_c*; do
  ID=$(basename $D | sed 's/Ca_fit_//')
  echo "=== $ID  $(python -c "import json;print(json.load(open('$D/ca_config.json'))['note'].split(':',1)[1].strip())" 2>/dev/null)"
  if [ -f $D/ca_result.json ]; then
    python - <<PY
import json; r=json.load(open('$D/ca_result.json'))
print(f"  total SSR (unweighted) = {r['ssr_unweighted']:.4e}   objective = {r['ssr_objective']:.4e}")
print('  ' + '  '.join(f"{k}={v:.4g}" for k,v in r['linear'].items()) + '   ' + '  '.join(f"{k}={v:.4g}(fixed)" for k,v in r['fixed'].items() if k not in r['params']))
print('  per-curve: ' + '  '.join(f"{k}={v:.3e}" for k,v in r['per_curve'].items()))
PY
    grep -h "ON BOUND" $D/slurm_*.out 2>/dev/null | sed 's/^/  /'
  else
    echo "  (no result yet)"; L=$(ls -t $D/slurm_*.out 2>/dev/null | head -1); [ -n "$L" ] && tail -1 "$L" | cut -c1-120
  fi
done
