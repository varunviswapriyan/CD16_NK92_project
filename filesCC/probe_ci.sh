#!/bin/bash
# probe_ci.sh -- inspect the pZAP replicate xlsx and the SSR machinery, so the
# bootstrap runners can be built to hit the right columns.
#   bash ~/CD16_NK92_project/filesCC/probe_ci.sh
set -e
XLSX=/home/gddaslab/share/Varun_Indrani/estimate_params_pzap/data/pZAP70_Tyr493_Tyr292_original_and_averages.xlsx
SRC=~/NK92_fit_v77/estimate_params_pzap_cleaned_up/pzap_param_estimation_NK92.py

echo "== pZAP replicate xlsx: sheet names, columns, first 3 rows =="
python3 - <<PY
import pandas as pd
try:
    d = pd.read_excel(r"$XLSX", sheet_name=None)
except Exception as e:
    print("ERROR reading xlsx:", e); raise SystemExit
for name, df in d.items():
    print(f'-- sheet: {name}  shape={df.shape} --')
    print('cols:', df.columns.tolist())
    print(df.head(3).to_string(index=False))
    print()
PY

echo
echo "== SSR / cost function machinery in $SRC =="
grep -n "^ratio\|ratio_\|SSR\|residue\|cost\|hetero_weight\|order_weight\|ratio_weight\|peak\|def [a-z]" "$SRC" | head -60
echo
echo "== how many rows the code processes: =="
grep -n "np\.array\|for.*time\|t_data\|obs.*=\|pzap_data" "$SRC" | head -25
