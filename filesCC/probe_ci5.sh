#!/bin/bash
# probe_ci5.sh -- dump all sheets of the pZAP xlsx so we can see how the raw
# 3-day values get transformed into the mean_zeta/gamma/hetero/NK92 columns
# that the fitting code actually reads.
#   bash ~/CD16_NK92_project/filesCC/probe_ci5.sh
XLSX=/home/gddaslab/share/Varun_Indrani/estimate_params_pzap/data/pZAP70_Tyr493_Tyr292_original_and_averages.xlsx

python3 - <<PY
import pandas as pd
x = pd.read_excel(r"$XLSX", sheet_name=None)
print('sheets:', list(x.keys()))
for name, df in x.items():
    print(f'\n===== sheet: {name}  shape={df.shape} =====')
    print('cols:', df.columns.tolist())
    print(df.to_string())
PY
