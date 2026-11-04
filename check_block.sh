#!/usr/bin/env bash
# check_block.sh -- the IFNGR-block (receptors=0 on T cells) runs finish instantly
# with cancer=0 -> almost certainly a divide-by-zero (bound/receptors) -> NaN rates.
# This cancels them, shows the evidence, and prints the source lines that divide
# by receptor number so the code can be patched properly.
# Run: cd ~/CD16_NK92_project && git pull --no-edit && bash check_block.sh
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || exit 1
scancel -u "$USER" -n blk_m,blk_f 2>/dev/null && echo "cancelled blk_m / blk_f"
L=$(ls runs_big/blockall/male/s*/log.spparks 2>/dev/null | head -1)
echo "=== $L : first 5 and last 5 data lines ==="
grep -A5 "sum(i1)" "$L" | head -7; echo ...; tail -5 "$L"
echo "=== nan/inf anywhere in blockall logs? ==="
grep -lis "nan\|inf" runs_big/blockall/*/s*/log.spparks | head -3
echo "=== source lines dividing by receptor counts ==="
for f in app_lattice.cpp app_rxn_diff_custom.cpp app_rxn_diff.cpp; do
  p=$(find . -path ./runs_big -prune -o -name "$f" -print | head -1); [ -n "$p" ] || continue
  echo "-- $p"
  grep -nE "/ *\(?[a-z_\[\]0-9]*(recept|\[1\])" "$p" | head -30
  grep -nE "receptors_per_cell|cell_receptors|total_receptors" "$p" | head -30
done
