#!/usr/bin/env bash
# find_mif_deep.sh -- thorough search for the raw mIF per-cell data (x, y, phenotype).
# Looks INSIDE files (first lines) for imaging-style column headers, not just at names.
# Searches: all of Giuseppe's archive (no depth limit) + the lab's shared folder.
# Run on log01: cd ~/CD16_NK92_project && git pull --no-edit && bash find_mif_deep.sh
# Takes a few minutes. Report: ~/CD16_NK92_project/mif_deep_report.txt
SQSH=/gpfs0/home1/gddaslab/lab/dead_users/rsgxg003/rsgxg003.sqsh
MNT=~/rsgxg003_mount
LAB=/gpfs0/home1/gddaslab/lab
OUT=~/CD16_NK92_project/mif_deep_report.txt
mkdir -p "$MNT"
[ -z "$(ls -A "$MNT" 2>/dev/null)" ] && { squashfuse "$SQSH" "$MNT" || { echo "mount failed"; exit 1; }; }

SKIP='dump\.|log\.spparks|/sample[0-9]+/|site-packages|anaconda|miniconda|/\.cache/|/\.local/|/\.git/|CMake|\.o$|\.so$|\.png$|\.pdf$|\.jpg$|\.svg$|\.cpp$|\.h$|\.py$|\.sh$|spk_'
HDR='CD8|CD3|CD4|TIM3|TIM-3|TCF1|TCF-1|Phenotype|phenotype|Centroid|Cell X|Cell Y|XMin|YMin|Cell ID|AR\+|Androgen|Opal|DAPI|Tissue Category'

scan(){ # root label
  local root=$1 label=$2
  echo "############ $label: $root"
  echo "--- file types (count by extension, simulation outputs excluded) ---"
  find "$root" -type f 2>/dev/null | grep -v -E "$SKIP" | sed -E 's/.*\.([A-Za-z0-9]{1,6})$/\1/; t; s/.*/(no ext)/' | sort | uniq -c | sort -rn | head -25
  echo "--- files whose first lines look like imaging tables ---"
  find "$root" -type f -size -2G 2>/dev/null | grep -v -E "$SKIP" | while read -r f; do
    case "$f" in *.xlsx|*.xls|*.zip|*.gz|*.qptiff|*.tif|*.tiff|*.mat|*.pkl|*.npy|*.h5) echo "  [binary/table, check by name] $f"; continue;; esac
    head -c 4000 "$f" 2>/dev/null | grep -q -E "$HDR" && echo "  MATCH $f"
  done | head -150
  echo "--- files > 1 MB (not simulation output) ---"
  find "$root" -type f -size +1M 2>/dev/null | grep -v -E "$SKIP" | xargs -r -d '\n' ls -la 2>/dev/null | sort -k5 -n -r | awk '{print "  "$5"  "$9}' | head -40
  echo
}
{
scan "$MNT" "Giuseppe archive"
echo "############ lab shared folder: top level"; ls -la "$LAB"
for d in "$LAB"/*; do
  [ "$d" = "$LAB/dead_users" ] && continue
  [ -d "$d" ] && scan "$d" "lab folder"
done
} 2>&1 | tee "$OUT"
echo
echo "Report: $OUT  (send Claude the MATCH lines and the >1 MB lists)"
echo "Unmount when done: fusermount -u $MNT"
