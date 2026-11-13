#!/usr/bin/env bash
# find_mif_data.sh -- look for the raw mIF slide data (x-y coordinates + phenotype per cell)
# in Giuseppe's archived home folder (SquashFS), which must be mounted first.
# The mount is only visible on the node you run this on (e.g. log01).
# Run: cd ~/CD16_NK92_project && git pull --no-edit && bash find_mif_data.sh
# Leaves the archive mounted at ~/rsgxg003_mount so you can look around.
# Unmount when done:  fusermount -u ~/rsgxg003_mount
SQSH=/gpfs0/home1/gddaslab/lab/dead_users/rsgxg003/rsgxg003.sqsh
MNT=~/rsgxg003_mount
OUT=~/CD16_NK92_project/mif_search_report.txt
mkdir -p "$MNT"
if [ -z "$(ls -A "$MNT" 2>/dev/null)" ]; then
  squashfuse "$SQSH" "$MNT" || { echo "mount failed"; exit 1; }
  echo ">> mounted $SQSH at $MNT"
else
  echo ">> already mounted at $MNT"
fi
{
echo "=== top level of Giuseppe's archive ==="
ls -la "$MNT"
echo
echo "=== what multiple_slide_reader.py / main.py read (file paths, open/read calls) ==="
for f in $(find "$MNT" -maxdepth 4 -name "multiple_slide_reader.py" -o -maxdepth 4 -name "main.py" 2>/dev/null | head -4); do
  echo "-- $f"
  grep -n -E "open\(|read_csv|read_excel|loadtxt|genfromtxt|json\.load|glob|listdir|path|\.csv|\.txt|\.xlsx|\.json" "$f" | head -40
done
echo
echo "=== folders whose names look like slide / mIF / imaging data ==="
find "$MNT" -maxdepth 6 -type d \( -iname "*slide*" -o -iname "*mif*" -o -iname "*imag*" -o -iname "*tissue*" -o -iname "*halo*" -o -iname "*inform*" -o -iname "*phenotyp*" -o -iname "*coord*" -o -iname "*data*" \) 2>/dev/null | grep -v -E "/\.|site-packages|anaconda|miniconda|/lib/" | head -60
echo
echo "=== data files mentioning day 7 / 14 / 21 (csv, txt, tsv, xlsx, json, npy, pkl) ==="
find "$MNT" -maxdepth 8 -type f \( -iname "*.csv" -o -iname "*.tsv" -o -iname "*.txt" -o -iname "*.xlsx" -o -iname "*.xls" -o -iname "*.json" -o -iname "*.npy" -o -iname "*.pkl" -o -iname "*.pickle" \) 2>/dev/null \
  | grep -v -E "/\.|site-packages|anaconda|miniconda|/lib/|CMake|log\.spparks|dump\." \
  | grep -i -E "d7|d14|d21|day|slide|mif|cell|phenotyp|coord|mouse|male|female|sex" | head -120
echo
echo "=== biggest data files anywhere (raw slide exports are usually large) ==="
find "$MNT" -maxdepth 8 -type f \( -iname "*.csv" -o -iname "*.tsv" -o -iname "*.txt" -o -iname "*.xlsx" \) -size +5M 2>/dev/null \
  | grep -v -E "site-packages|anaconda|miniconda|log\.spparks|dump\." | xargs -r ls -la 2>/dev/null | sort -k5 -n -r | head -30
} 2>&1 | tee "$OUT"
echo
echo "Report saved to $OUT. Send Claude a screenshot of it."
echo "Archive is still mounted at $MNT. Unmount when done: fusermount -u $MNT"
