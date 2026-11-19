#!/usr/bin/env bash
# find_area_method.sh -- find how Giuseppe computed tumour area from the mIF slides (thesis Fig 4.3).
# Searches his archived home folder for code that computes an area (convex hull, alpha shape,
# polygon, "tumor area", pixel counting...) and prints the matching lines.
# Run on log01: cd ~/CD16_NK92_project && git pull --no-edit && bash find_area_method.sh
# Report: ~/CD16_NK92_project/area_method_report.txt
MNT=~/rsgxg003_mount; mkdir -p $MNT
mountpoint -q $MNT || squashfuse /gpfs0/home1/gddaslab/lab/dead_users/rsgxg003/rsgxg003.sqsh $MNT || { echo "mount failed"; exit 1; }
PAT='tumou?r_?area|tumou?rArea|area_of_tum|ConvexHull|convex_hull|alphashape|alpha_shape|concave|Polygon\(|shapely|contourArea|find_contours|binary_fill_holes|tumou?r[ _]size|slide_area|calc[A-Za-z_]*[Aa]rea|[Aa]rea[A-Za-z_]*\(.*(x|X)'
{
echo "=== code files that compute an area (Python, notebooks, R, MATLAB) ==="
find $MNT \( -name "sample[0-9]*" -o -name .git -o -name site-packages -o -name anaconda3 -o -name miniconda3 -o -name .cache -o -name .local \) -prune -o \
     -type f \( -name "*.py" -o -name "*.ipynb" -o -name "*.R" -o -name "*.m" \) -print 2>/dev/null \
  | xargs -r -d '\n' grep -l -E -i "$PAT" 2>/dev/null | grep -v -i aaron_yang | tee /tmp/area_files_$$.txt
echo
echo "=== the matching lines in each file ==="
while read -r f; do
  echo "-- $f"
  grep -n -E -i "$PAT" "$f" 2>/dev/null | cut -c1-220 | head -15
done < /tmp/area_files_$$.txt
rm -f /tmp/area_files_$$.txt
echo
echo "=== output files whose names mention area (results he may have saved) ==="
find $MNT \( -name "sample[0-9]*" -o -name .git -o -name site-packages \) -prune -o -type f -iname "*area*" -print 2>/dev/null | grep -v -i -E "aaron_yang|/lib/|\.so$" | head -40
} 2>&1 | tee ~/CD16_NK92_project/area_method_report.txt
echo; echo "Report: ~/CD16_NK92_project/area_method_report.txt   Unmount when done: fusermount -u $MNT"
