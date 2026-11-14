#!/usr/bin/env bash
# find_init_builder.sh -- find the script Giuseppe used to build the SPPARKS_init slide files
# from the raw mIF data, and the raw data it reads. Searches the FULL archive, skipping
# Aaron Yang's (human) project and simulation sample folders so nothing gets crowded out.
# Run on log01: cd ~/CD16_NK92_project && git pull --no-edit && bash find_init_builder.sh
# Report: ~/CD16_NK92_project/init_builder_report.txt
SQSH=/gpfs0/home1/gddaslab/lab/dead_users/rsgxg003/rsgxg003.sqsh
MNT=~/rsgxg003_mount
OUT=~/CD16_NK92_project/init_builder_report.txt
mkdir -p "$MNT"
[ -z "$(ls -A "$MNT" 2>/dev/null)" ] && { squashfuse "$SQSH" "$MNT" || { echo "mount failed"; exit 1; }; }
PRUNE=( \( -name Aaron_Yang_Data_Analysis -o -name "sample[0-9]*" -o -name .cache -o -name .local -o -name .git -o -name site-packages -o -name anaconda3 -o -name miniconda3 \) -prune -o )
{
echo "=== HOW_TO_MOUNT.txt (Naresh) ==="
cat /gpfs0/home1/gddaslab/lab/dead_users/rsgxg003/HOW_TO_MOUNT.txt 2>/dev/null || cat /home/gddaslab/lab/dead_users/rsgxg003/HOW_TO_MOUNT.txt 2>/dev/null
echo
echo "=== 1. scripts that mention SPPARKS_init (the slide-file builders) ==="
find "$MNT" "${PRUNE[@]}" -type f \( -name "*.py" -o -name "*.ipynb" -o -name "*.cpp" -o -name "*.m" -o -name "*.R" -o -name "*.sh" -o -name "*.jl" \) -print 2>/dev/null \
  | xargs -r -d '\n' grep -l -E "SPPARKS_init|_init_file|write.*Values|Site file written" 2>/dev/null | tee /tmp/builders_$$.txt
echo
echo "=== 2. what those scripts read (paths / read calls) ==="
while read -r f; do
  echo "-- $f"
  grep -n -E "read_csv|read_excel|read_table|loadtxt|genfromtxt|open\(.*['\"]r|glob\(|listdir|_path *=|_dir *=|\.csv|\.tsv|\.xlsx|\.json|mIF|slide" "$f" 2>/dev/null | grep -v "^\s*#" | head -25
done < /tmp/builders_$$.txt
rm -f /tmp/builders_$$.txt
echo
echo "=== 3. files / folders named like mouse mIF data ==="
find "$MNT" "${PRUNE[@]}" \( -iname "*mb49*" -o -iname "*mif*" -o -iname "*multiplex*" -o -iname "*inform*" -o -iname "*halo*" -o -iname "*cell_seg*" -o -iname "*phenotyp*" -o -iname "*zihai*" -o -iname "*d14m*" -o -iname "*d14f*" -o -iname "*d7m*" -o -iname "*d7f*" -o -iname "*d21m*" -o -iname "*d21f*" -o -iname "*d7_*" -o -iname "*d21_*" -o -iname "*day7*" -o -iname "*day21*" -o -iname "*mouse*" \) -print 2>/dev/null \
  | grep -v -E "dump\.|log\.spparks|_100_samples|/sample[0-9]" | head -100
echo
echo "=== 4. large non-simulation data files (>1 MB), outside Aaron Yang ==="
find "$MNT" "${PRUNE[@]}" -type f -size +1M \( -iname "*.csv" -o -iname "*.tsv" -o -iname "*.txt" -o -iname "*.xlsx" -o -iname "*.xls" -o -iname "*.json" -o -iname "*.zip" -o -iname "*.gz" -o -iname "*.tar" -o -iname "*.pkl" -o -iname "*.npy" -o -iname "*.h5" -o -iname "*.mat" \) -print 2>/dev/null \
  | grep -v -E "dump\.|log\.spparks" | xargs -r -d '\n' ls -la 2>/dev/null | awk '{print "  "$5"  "$9}' | sort -n -r | head -50
} 2>&1 | tee "$OUT"
echo
echo "Report: $OUT   Unmount when done: fusermount -u $MNT"
