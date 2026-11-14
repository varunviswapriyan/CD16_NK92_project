#!/usr/bin/env bash
# find_share_mif.sh -- look for the mouse mIF slide data in the lab's shared folder
# (/home/gddaslab/share), which our earlier searches of Giuseppe's archive did not cover.
# Run: cd ~/CD16_NK92_project && git pull --no-edit && bash find_share_mif.sh
# Report: ~/CD16_NK92_project/share_mif_report.txt
S=/home/gddaslab/share
OUT=~/CD16_NK92_project/share_mif_report.txt
{
echo "=== 2DCancerModelData folders ==="
for d in "$S/Rajdeep_Giuseppe/2DCancerModelData" "$S/Darren_Giuseppe/2DCancerModelData"; do
  echo "-- $d"; ls "$d" 2>/dev/null | head -20
done
echo
echo "=== names that look like mouse slide / mIF data anywhere in share ==="
find "$S" -maxdepth 6 \( -iname "*mb49*" -o -iname "*mif*" -o -iname "*d14*" -o -iname "*d7_*" -o -iname "*d21*" \
  -o -iname "*day*" -o -iname "*slide*" -o -iname "*mouse*" -o -iname "*sex*" -o -iname "*zihai*" \
  -o -iname "*phenotyp*" -o -iname "*cell_seg*" -o -iname "*SPPARKS_init*" \) 2>/dev/null | head -80
echo
echo "=== data files in share whose first lines have imaging column headers ==="
find "$S" -maxdepth 6 -type f -size -2G \( -iname "*.csv" -o -iname "*.tsv" -o -iname "*.txt" \) 2>/dev/null \
  | grep -v -E "cancer_only|dump\.|log\.spparks" | while read -r f; do
    head -c 4000 "$f" 2>/dev/null | grep -q -E "CD8|CD3|TIM3|TCF1|Phenotype|Centroid|Cell X|XMin|DAPI|Opal" && echo "MATCH $f"
  done | head -60
} 2>&1 | tee "$OUT"
echo
echo "Report: $OUT"
