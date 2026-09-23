#!/bin/bash
# run_bootstrap.sh -- launch both bootstrap campaigns for parameter CIs.
#
#   bash ~/CD16_NK92_project/filesCC/run_bootstrap.sh check    # dry run: build only, no jobs
#   bash ~/CD16_NK92_project/filesCC/run_bootstrap.sh ca       # 25 Ca samples (~1 h)
#   bash ~/CD16_NK92_project/filesCC/run_bootstrap.sh pzap     # 12 emp + 12 par (overnight)
#   bash ~/CD16_NK92_project/filesCC/run_bootstrap.sh all      # both
#   bash ~/CD16_NK92_project/filesCC/run_bootstrap.sh report   # CI tables for both
set -e
F=~/CD16_NK92_project/filesCC
MODE="${1:-all}"

case "$MODE" in
  check)
    echo "== Ca: single sample, local, short (sanity check) =="
    cd ~/Ca_fit_c02 && timeout 900 python $F/bootstrap_ca.py 1 || echo "(timed out or failed -- see above)"
    echo
    echo "== pZAP: build 1 emp + 1 par, no submit =="
    python $F/bootstrap_pzap.py build 1 1
    echo
    echo "Looks right? then: bash \$0 all"
    ;;
  ca)
    cd ~/Ca_fit_c02 && python $F/bootstrap_ca.py submit 25
    ;;
  pzap)
    python $F/bootstrap_pzap.py submit 12 12
    ;;
  all)
    cd ~/Ca_fit_c02 && python $F/bootstrap_ca.py submit 25
    echo
    python $F/bootstrap_pzap.py submit 12 12
    echo
    echo "All submitted. squeue -u \$USER"
    ;;
  report)
    echo "################ Ca (n01) ################"
    cd ~/Ca_fit_c02 && python $F/bootstrap_ca.py report || true
    echo
    echo "############### pZAP (v77) ###############"
    python $F/bootstrap_pzap.py report || true
    ;;
  *)
    echo "usage: $0 {check|ca|pzap|all|report}"; exit 1 ;;
esac
