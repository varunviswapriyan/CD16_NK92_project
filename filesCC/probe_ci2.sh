#!/bin/bash
# probe_ci2.sh -- inspect how the pZAP script loads the experimental targets,
# so the bootstrap knows what to resample.
#   bash ~/CD16_NK92_project/filesCC/probe_ci2.sh
SRC=~/NK92_fit_v77/estimate_params_pzap_cleaned_up/pzap_param_estimation_NK92.py
CSSR=~/NK92_fit_v77/estimate_params_pzap_cleaned_up/calculate_SSR_gamma.py

echo "== how pzap_param_estimation_NK92.py builds the EXP dict =="
grep -n "EXP\[\|EXP =\|exp\['\|exp =\|m_own\|z_own\|g_own\|m_gam\|z_gam\|g_gam\|hg\|zg\|extras\|read_csv\|read_excel\|np\.array" "$SRC" | head -60

echo
echo "== calculate_SSR_gamma.py: what does the residue function need =="
head -80 "$CSSR"
