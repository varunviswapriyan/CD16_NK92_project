#!/bin/bash
# show_bands.sh -- generate the prediction-band figures for the Ca variants and
# serve them, alongside the pZAP CI histograms.
#   bash ~/CD16_NK92_project/filesCC/show_bands.sh                  # default set
#   bash ~/CD16_NK92_project/filesCC/show_bands.sh hillonly Sfix    # pick variants
set -e
F=~/CD16_NK92_project/filesCC
VARIANTS="${@:-hillonly min4h k3fix linonly}"

cd ~/Ca_fit_c02
for V in $VARIANTS; do
  echo "== band: $V =="
  python $F/bootstrap_ca2.py band "$V" || echo "  (skipped $V)"
done

# collect everything worth looking at into one folder
OUT=~/ci_plots
mkdir -p $OUT
cp ~/Ca_fit_c02/out_boot_ca2/band_*.png            $OUT/ 2>/dev/null || true
cp ~/Ca_fit_c02/out_boot_ca/ci_hist_ca.png         $OUT/ca_n01_hist.png 2>/dev/null || true
cp ~/boot_pzap/ci_hist_pzap.png                    $OUT/pzap_hist.png 2>/dev/null || true
cp ~/Ca_fit_c02/out_filesD/plot_ca_n01_S.png       $OUT/ca_n01_fit.png 2>/dev/null || true
cp ~/NK92_fit_v77/estimate_params_pzap_cleaned_up/plot_N.png $OUT/pzap_v77_fit.png 2>/dev/null || true

echo
ls -1 $OUT
IP=$(hostname -I | awk '{for(i=1;i<=NF;i++) if($i ~ /^10\.73\./) print $i}')
echo
echo "Open:  http://${IP:-10.73.170.129}:8000/"
cd $OUT && python -m http.server 8000
