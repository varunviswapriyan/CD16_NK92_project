#!/bin/bash
# probe_ci3.sh -- final probe: how are experimental peak-ratio targets loaded?
#   bash ~/CD16_NK92_project/filesCC/probe_ci3.sh
SRC=~/NK92_fit_v77/estimate_params_pzap_cleaned_up
echo "== calculate_SSR_gamma.py first 40 lines =="
sed -n '1,40p' "$SRC/calculate_SSR_gamma.py"
echo
echo "== calculate_residue definition + EXP usage =="
grep -n "def calculate_residue\|EXP\[\|exp\[" "$SRC/calculate_SSR_gamma.py"
echo
echo "== EXP / own / gam references in pzap_param_estimation_NK92.py =="
grep -n "own\|gam\|EXP\|read_csv\|read_excel" "$SRC/pzap_param_estimation_NK92.py" | head -25
echo
echo "== also SSR_zeta and SSR_mixed just in case =="
grep -n "EXP\|exp\[" "$SRC/calculate_SSR_zeta.py" "$SRC/calculate_SSR_mixed.py" | head -20
