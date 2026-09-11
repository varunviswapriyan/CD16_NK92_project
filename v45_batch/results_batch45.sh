#!/bin/bash
# results_batch45.sh
#   bash results_batch45.sh                 -> residues + config + model patch for v45..v54
#   bash results_batch45.sh plot 48 50 52   -> run check_fit.py on chosen versions, copy plots to ~/plot_vNN.png
HOME_DIR=/home/gddaslab/vxv016
SUB=estimate_params_pzap_cleaned_up

if [ "$1" == "plot" ]; then
  shift
  for V in "$@"; do
    echo "=== plotting v$V ==="
    cd $HOME_DIR/NK92_fit_v$V/$SUB || continue
    python check_fit.py 2>&1 | grep -iE "lig0=|peak|SSR"
    cp plot_N.png $HOME_DIR/plot_v$V.png && echo "    -> ~/plot_v$V.png"
  done
  echo; echo "View: cd ~ && python -m http.server 8000  (browser: http://<login-node-ip>:8000/)"
  exit 0
fi

for V in 45 46 47 48 49 50 51 52 53 54; do
  D=$HOME_DIR/NK92_fit_v$V/$SUB
  F=$D/analysis_param_residue.dat
  echo "=== v$V ==="
  [ -f $D/v_config.json ] && grep -oE '"NOTE": "[^"]+"' $D/v_config.json | cut -d'"' -f4
  if [ -f "$F" ]; then
    head -1 "$F"
    grep "fourth root" "$F"
    grep "^config" "$F" | cut -c1-160
  else
    echo "  (no result yet)"
    L=$(ls -t $HOME_DIR/slurm_v${V}_*.out 2>/dev/null | head -1)
    [ -n "$L" ] && grep -oE "[0-9]+/[0-9]+, best_cost=[0-9.]+" "$L" | tail -1
  fi
  echo
done
echo "Next: bash $(dirname $0)/results_batch45.sh plot 45 48 50 52   (the 'hetero/gamma peak' line is the diagnostic)"
