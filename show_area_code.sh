#!/usr/bin/env bash
# show_area_code.sh -- show how Giuseppe's analysisCalcs.py computes tumour area / contours.
# Run on log01: cd ~/CD16_NK92_project && git pull --no-edit && bash show_area_code.sh
# Saves everything to ~/CD16_NK92_project/area_code_report.txt
MNT=~/rsgxg003_mount; mkdir -p $MNT
mountpoint -q $MNT || squashfuse /gpfs0/home1/gddaslab/lab/dead_users/rsgxg003/rsgxg003.sqsh $MNT || { echo "mount failed"; exit 1; }
F=$MNT/melanoIMCProject/analysisCalcs.py
{
echo "=== function names and lines about area / flood fill / contour / grid ==="
grep -n -E "^def |flood|[Aa]rea|tumorContour *=|pixel|binSize|gridSize|30" "$F" | head -90
echo
echo "=== code around each 'area' calculation (5 lines before, 8 after) ==="
grep -n -E "[Aa]rea *=|[Aa]rea\(|\.area" "$F" | head -8 | cut -d: -f1 | while read -r n; do
  echo "--- line $n"; sed -n "$((n>5?n-5:1)),$((n+8))p" "$F"
done
} 2>&1 | tee ~/CD16_NK92_project/area_code_report.txt
echo; echo "Saved: ~/CD16_NK92_project/area_code_report.txt"
