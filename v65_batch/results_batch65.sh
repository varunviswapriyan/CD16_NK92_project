#!/bin/bash
# results_batch65.sh
#   bash results_batch65.sh                    -> residues + params for v65..v74
#   bash results_batch65.sh plot 66 68 72      -> run check_fit.py on chosen versions, copy plots to ~/plot_vNN.png
HOME_DIR=/home/gddaslab/vxv016
SUB=estimate_params_pzap_cleaned_up

if [ "$1" == "plot" ]; then
  shift
  for V in "$@"; do
    echo "=== plotting v$V ==="
    cd $HOME_DIR/NK92_fit_v$V/$SUB || continue
    python check_fit.py 2>&1 | grep -E "^Linear|peak ratios|Total SSR|t=0|Largest"
    cp plot_N.png $HOME_DIR/plot_v$V.png && echo "    -> ~/plot_v$V.png"
  done
  echo; echo "View: cd ~ && python -m http.server 8000  (browser: http://<login-node-ip>:8000/)"
  exit 0
fi

for V in 65 66 67 68 69 70 71 72 73 74; do
  D=$HOME_DIR/NK92_fit_v$V/$SUB
  F=$D/analysis_param_residue.dat
  echo "=== v$V ==="
  [ -f $D/v_config.json ] && grep -oE '"NOTE": "[^"]+"' $D/v_config.json | cut -d'"' -f4
  if [ -f "$F" ]; then
    grep "fourth root" "$F"
    grep "^linear" "$F"
  else
    echo "  (no result yet)"
    L=$(ls -t $HOME_DIR/slurm_v${V}_*.out 2>/dev/null | head -1)
    [ -n "$L" ] && grep -oE "[0-9]+/[0-9]+, best_cost=[0-9.]+" "$L" | tail -1
  fi
  echo
done
echo "Next: bash $(dirname $0)/results_batch65.sh plot 66 68 69 72   ('Largest single-point miss' + point-by-point lines are the diagnostic)"
