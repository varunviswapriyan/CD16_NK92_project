#!/bin/bash
# final_plots.sh -- collect the two final figures (pZAP v77, Ca n03) and serve them.
#   bash ~/CD16_NK92_project/filesCC/final_plots.sh
set -e
OUT=~/final_plots
mkdir -p $OUT
cp ~/NK92_fit_v77/estimate_params_pzap_cleaned_up/plot_N.png $OUT/pzap_v77.png
cp ~/Ca_fit_c02/out_filesD/plot_ca_n03_S_nn.png              $OUT/ca_n03.png
cp ~/Ca_fit_c02/out_filesD/plot_ca_n01_S.png                 $OUT/ca_n01_minimal.png 2>/dev/null || true

# right-panel-only crop of the Ca figure (fit vs data, no input panel)
python3 - <<'PY' 2>/dev/null || echo "(crop skipped: PIL not available)"
from PIL import Image
im = Image.open('/dev/stdin'.replace('/dev/stdin', __import__('os').path.expanduser('~/final_plots/ca_n03.png')))
w, h = im.size
im.crop((int(w*0.505), 0, w, h)).save(__import__('os').path.expanduser('~/final_plots/ca_n03_fit_only.png'))
print("wrote ca_n03_fit_only.png (right panel crop)")
PY

IP=$(hostname -I | awk '{for(i=1;i<=NF;i++) if($i ~ /^10\.73\./) print $i}')
echo
ls -1 $OUT
echo
echo "Open:  http://${IP:-10.73.170.128}:8000/"
cd $OUT && python -m http.server 8000
