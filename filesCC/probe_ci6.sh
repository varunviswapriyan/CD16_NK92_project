#!/bin/bash
# probe_ci6.sh -- does the fitting code use the mean_NK92 (KI 6) column?
#   bash ~/CD16_NK92_project/filesCC/probe_ci6.sh
D=~/NK92_fit_v77/estimate_params_pzap_cleaned_up
grep -n "mean_NK92\|NK92_\|KI 6\|KI6" $D/*.py
echo "--- (empty above = KI 6 unused; bootstrap ignores it)"
