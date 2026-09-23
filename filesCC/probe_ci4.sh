#!/bin/bash
# probe_ci4.sh -- final data-format probe before bootstrap runners.
#   bash ~/CD16_NK92_project/filesCC/probe_ci4.sh
set -e
CSV=~/NK92_fit_v77/estimate_params_pzap_cleaned_up/data/pZAP70_Tyr493_mean.csv
XLSX=/home/gddaslab/share/Varun_Indrani/estimate_params_pzap/data/pZAP70_Tyr493_Tyr292_original_and_averages.xlsx

echo "== pZAP70_Tyr493_mean.csv (first 5 rows) =="
if [ -f "$CSV" ]; then
  head -5 "$CSV"
else
  echo "NOT FOUND: $CSV"
  find ~/NK92_fit_v77 -name "*mean*.csv" 2>/dev/null
fi

echo
echo "== Original_values sheet, all 24 rows =="
python3 - <<PY
import pandas as pd
d = pd.read_excel(r"$XLSX", sheet_name='Original_values')
print('cols:', d.columns.tolist())
print('shape:', d.shape)
print(d.to_string())
PY

echo
echo "== also check ca_data csv format one more time =="
head -3 ~/Ca_fit_c02/ca_data/Ca_NK92.csv
