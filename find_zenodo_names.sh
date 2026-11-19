#!/usr/bin/env bash
# find_zenodo_names.sh -- final check: are the Zenodo slide files (d7_M_1.csv ... d21_F_5.csv, folder
# Processed_Zihai_Data) anywhere on HPC under these exact names? Searches Giuseppe's archive, the lab
# share, the lab folder and every lab member's home folder we can read. Our own copy is skipped.
# Run on log01: cd ~/CD16_NK92_project && git pull --no-edit && bash find_zenodo_names.sh
MNT=~/rsgxg003_mount
mkdir -p $MNT
mountpoint -q $MNT || squashfuse /gpfs0/home1/gddaslab/lab/dead_users/rsgxg003/rsgxg003.sqsh $MNT || echo "!! could not mount Giuseppe's archive"
ROOTS="$MNT /home/gddaslab/share /home/gddaslab/lab $(ls -d /home/gddaslab/*/ 2>/dev/null | grep -v -E '/(lab|share)/$')"
echo "searching: $ROOTS" | tr ' ' '\n' | sed 's/^/  /'
echo
echo "=== files/folders named like the Zenodo data ==="
find $ROOTS \( -name "sample[0-9]*" -o -name .git -o -name site-packages \) -prune -o \
     \( -iname "*Zihai_Data*" -o -regex '.*/d\(7\|14\|21\)_[MF]_[0-9]\.\(csv\|txt\|xlsx\|tsv\)' \) -print 2>/dev/null \
  | grep -v -E "vxv016/CD16_NK92_project|SPPARKS_init" | head -60
echo "=== (end of list; nothing above = not on HPC under these names) ==="
echo
echo "=== unreadable home folders (could not be searched) ==="
for d in /home/gddaslab/*/; do [ -r "$d" ] && [ -x "$d" ] || echo "  $d"; done
echo
echo "Unmount when done: fusermount -u $MNT"
