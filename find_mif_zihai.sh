#!/usr/bin/env bash
# find_mif_zihai.sh -- look inside Giuseppe's sex_bias_Zihai folder (and a few other likely
# folders) for the raw mIF slide data: per-cell tables with x, y and phenotype.
# Run on log01: cd ~/CD16_NK92_project && git pull --no-edit && bash find_mif_zihai.sh
SQSH=/gpfs0/home1/gddaslab/lab/dead_users/rsgxg003/rsgxg003.sqsh
MNT=~/rsgxg003_mount
OUT=~/CD16_NK92_project/mif_zihai_report.txt
mkdir -p "$MNT"
[ -z "$(ls -A "$MNT" 2>/dev/null)" ] && { squashfuse "$SQSH" "$MNT" || { echo "mount failed"; exit 1; }; }
Z="$MNT/sex_bias_Zihai"
{
echo "=== sex_bias_Zihai top level ==="
ls -la "$Z"
echo
echo "=== what multiple_slide_reader.py reads ==="
grep -n -E "open\(|read_csv|loadtxt|genfromtxt|json|glob|listdir|_path|\.csv|\.txt|\.json" "$Z/multiple_slide_reader.py" | head -40
echo
echo "=== folders in sex_bias_Zihai ==="
(cd "$Z" && find . -maxdepth 3 -type d | grep -v -E "src|cmake|\.git|_out|run" | head -60)
echo
echo "=== data files in sex_bias_Zihai ==="
(cd "$Z" && find . -maxdepth 4 -type f \( -iname "*.csv" -o -iname "*.json" -o -iname "*.xlsx" -o -iname "*.txt" \) | grep -v -E "log\.|dump\.|CMake" | head -60)
echo
echo "=== other likely folders ==="
for d in old_and_reference_stuff porting_out_to_new_computer plots_going_out; do echo "-- $d"; ls "$MNT/$d"; done
echo
echo "=== earlier report lines mentioning zihai / sex bias / mb49 / days ==="
grep -i -E "zihai|sex_bias|mb49|d7|day" ~/CD16_NK92_project/mif_search_report.txt 2>/dev/null | head -40
} 2>&1 | tee "$OUT"
echo
echo "Report saved to $OUT. Archive still mounted; unmount when done: fusermount -u $MNT"
