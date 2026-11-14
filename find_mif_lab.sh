#!/usr/bin/env bash
# find_mif_lab.sh -- last sweep for the raw mouse mIF slide data: every lab member's home folder
# under /home/gddaslab and the rest of /home/gddaslab/lab (Giuseppe's .sqsh and the share folder
# were already searched). Unreadable folders are skipped and listed at the end.
# Run: cd ~/CD16_NK92_project && git pull --no-edit && bash find_mif_lab.sh
# Report: ~/CD16_NK92_project/mif_lab_report.txt
OUT=~/CD16_NK92_project/mif_lab_report.txt
ROOTS=$(ls -d /home/gddaslab/*/ /home/gddaslab/lab/*/ 2>/dev/null | grep -v -E "/share/$|/dead_users/$|/vxv016/$")
NAMES='-iname *mb49* -o -iname *mif* -o -iname *multiplex* -o -iname *inform* -o -iname *halo* -o -iname *cell_seg* -o -iname *phenotyp* -o -iname *zihai* -o -iname *sex_bias* -o -iname *sexbias* -o -iname *d14m* -o -iname *d14f* -o -iname *d7m* -o -iname *d7f* -o -iname *d21m* -o -iname *d21f* -o -iname *d7_* -o -iname *d21_* -o -iname *SPPARKS_init*'
{
echo "=== folders searched ==="; echo "$ROOTS"
echo
echo "=== unreadable (could not search) ==="
for r in $ROOTS; do [ -r "$r" ] && [ -x "$r" ] || echo "  $r"; done
echo
echo "=== names that look like mouse slide / mIF data ==="
for r in $ROOTS; do
  [ -r "$r" ] || continue
  set -f; find "$r" -maxdepth 8 \( -name .git -o -name "sample[0-9]*" -o -name site-packages -o -name .cache \) -prune -o \( $NAMES \) -print 2>/dev/null; set +f
done | grep -v -E "dump\.|log\.spparks|IMC_init|responder" | head -100
echo
echo "=== data files whose first lines have imaging column headers ==="
for r in $ROOTS; do
  [ -r "$r" ] || continue
  find "$r" -maxdepth 8 \( -name .git -o -name "sample[0-9]*" -o -name site-packages -o -name .cache \) -prune -o -type f -size -2G \( -iname "*.csv" -o -iname "*.tsv" -o -iname "*.txt" \) -print 2>/dev/null
done | grep -v -E "dump\.|log\.spparks|slurm|CyGMM" | while read -r f; do
  head -c 4000 "$f" 2>/dev/null | grep -q -E "TIM3|TCF1|Phenotype|Centroid|Cell X|XMin|DAPI|Opal|CD3\+|CD8\+" && echo "MATCH $f"
done | head -80
} 2>&1 | tee "$OUT"
echo
echo "Report: $OUT"
