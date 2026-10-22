#!/usr/bin/env bash
# copy_giuseppe_ref.sh -- mount Giuseppe's archived SquashFS, copy the input
# decks + major code + slides + one src/binary into ~/giuseppe_ref, confirm our
# prod in.full matches his, list his slides, then unmount.
# NOTE: a squashfuse mount is visible only on the machine you run this on
#       (e.g. log01). Run it on one login node and stay there.
# Run:  cd ~/CD16_NK92_project && git pull && bash copy_giuseppe_ref.sh
set -uo pipefail

SQSH=/gpfs0/home1/gddaslab/lab/dead_users/rsgxg003/rsgxg003.sqsh
MNT=~/rsgxg003_mount
SRC="$MNT/sex_bias_Zihai"
DST=~/giuseppe_ref
OURS=~/CD16_NK92_project/SPPARKS_Sex_Bias/rsc_out_2026-07-03/in.full

mkdir -p "$MNT" "$DST"

# mount only if not already mounted
if ! mountpoint -q "$MNT" 2>/dev/null && [ -z "$(ls -A "$MNT" 2>/dev/null)" ]; then
  echo ">> mounting squashfs ..."
  squashfuse "$SQSH" "$MNT" || { echo "mount failed"; exit 1; }
else
  echo ">> already mounted"
fi

echo ">> copying input decks + code + slides + src/binary ..."
cp -v "$SRC"/in.* "$DST"/ 2>/dev/null
cp -v "$SRC"/main.py "$SRC"/multiple_slide_reader.py "$SRC"/run_samples.py \
      "$SRC"/DumpFileReader.py "$SRC"/auto_run.sh "$SRC"/run_sim_0.sh "$DST"/ 2>/dev/null
cp -rv "$SRC"/sex_bias_initialization_files "$DST"/ 2>/dev/null
cp -v "$SRC"/spk_redsky "$SRC"/spk_redsky_O_2024-12-3_CD4s_ARSink "$DST"/ 2>/dev/null
cp -rv "$SRC"/src_2024-11-18_post-Zihai_adding_CD4s "$DST"/ 2>/dev/null

echo
echo "=== does his in.full match our prod deck? ==="
if diff "$DST/in.full" "$OURS" >/dev/null 2>&1; then
  echo "SAME DECK -- prod runs use exactly Giuseppe's in.full"
else
  echo "THEY DIFFER -- diff below:"
  diff "$DST/in.full" "$OURS"
fi

echo
echo "=== his day-14 initialization slides ==="
ls -la "$DST/sex_bias_initialization_files" 2>/dev/null

echo
echo ">> unmounting ..."
fusermount -u "$MNT" && echo ">> done. reference copy is in $DST (persists on GPFS)"
