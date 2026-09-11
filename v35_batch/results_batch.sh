#!/bin/bash
# results_batch.sh
#   bash results_batch.sh              -> print residues + config for v35..v44
#   bash results_batch.sh plot 35 40   -> run check_fit.py on v35 and v40, copy plots to ~/plot_vNN.png
HOME_DIR=/home/gddaslab/vxv016
SUB=estimate_params_pzap_cleaned_up

if [ "$1" == "plot" ]; then
  shift
  for V in "$@"; do
    echo "=== plotting v$V ==="
    cd $HOME_DIR/NK92_fit_v$V/$SUB || continue
    python check_fit.py 2>&1 | grep -E "lig0=|peak ratios|Total SSR"
    cp plot_N.png $HOME_DIR/plot_v$V.png && echo "    -> ~/plot_v$V.png"
  done
  echo; echo "View: cd ~ && python -m http.server 8000  (browser: http://<login-node-ip>:8000/)"
  exit 0
fi

for V in 35 36 37 38 39 40 41 42 43 44; do
  F=$HOME_DIR/NK92_fit_v$V/$SUB/analysis_param_residue.dat
  echo "=== v$V ==="
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
