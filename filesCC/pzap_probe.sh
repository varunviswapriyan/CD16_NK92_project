#!/bin/bash
# =============================================================================
# pzap_probe.sh -- read-only. Shows what a finished pZAP bootstrap sample holds,
# so the prediction-band plot can be written against the real file layout
# instead of a guessed one.
#
#   bash ~/CD16_NK92_project/filesCC/pzap_probe.sh
#
# Submits nothing, changes nothing, writes nothing.  Takes seconds.
# =============================================================================
ROOT=${ROOT:-$HOME/boot_pzap_conv}
SUB=estimate_params_pzap_cleaned_up
SRC=$HOME/NK92_fit_v77

echo "root: $ROOT"
[ -d "$ROOT" ] || { echo "not found -- set ROOT=~/boot_pzap instead"; exit 1; }

# a sample that actually finished
S=""
for d in "$ROOT"/emp[0-9][0-9][0-9]; do
  [ -f "$d/$SUB/analysis_param_residue.dat" ] || continue
  grep -qi '^linear' "$d/$SUB/analysis_param_residue.dat" 2>/dev/null && { S="$d"; break; }
done
[ -n "$S" ] || { echo "no finished sample found under $ROOT"; exit 1; }
echo "using finished sample: $(basename "$S")"

echo
echo "================ 1. files in that sample (depth 3) ================"
find "$S" -maxdepth 3 \( -type f -o -type l \) 2>/dev/null \
  | sed "s|^$S/||" | sort | head -60
echo "  ... $(find "$S" -type f 2>/dev/null | wc -l) files total"

echo
echo "================ 2. every CSV, with its header ================"
find "$S" -name '*.csv' 2>/dev/null | sort | head -20 | while read -r f; do
  echo "--- ${f#$S/}   ($(wc -l < "$f") lines)"
  head -2 "$f"
done

echo
echo "================ 3. any .gdat / .cdat / .scan model output ================"
find "$S" \( -name '*.gdat' -o -name '*.cdat' -o -name '*.scan' \) 2>/dev/null \
  | sed "s|^$S/||" | sort | head -20
f=$(find "$S" \( -name '*.gdat' -o -name '*.cdat' \) 2>/dev/null | head -1)
if [ -n "$f" ]; then echo "--- head of ${f#$S/}"; head -3 "$f"; fi

echo
echo "================ 4. analysis_param_residue.dat ================"
cat "$S/$SUB/analysis_param_residue.dat" 2>/dev/null | head -20

echo
echo "================ 5. check_fit.py (this already knows how to plot) ======"
CF=$(find "$SRC" "$S" -name 'check_fit.py' 2>/dev/null | head -1)
if [ -n "$CF" ]; then
  echo "--- $CF  ($(wc -l < "$CF") lines)"
  head -80 "$CF"
else
  echo "check_fit.py not found under $SRC or the sample"
  echo "other python in the fit folder:"
  ls -1 "$SRC/$SUB"/*.py 2>/dev/null | sed 's|.*/|  |'
fi

echo
echo "================ 6. the experimental data file ================"
D="$SRC/$SUB/data/pZAP70_Tyr493_mean.csv"
[ -f "$D" ] && { echo "--- $D"; cat "$D"; } || echo "not found: $D"

echo
echo "================ 7. is there an SE/SD file to draw error bars? ========"
ls -1 "$SRC/$SUB/data/" 2>/dev/null | sed 's|^|  |'
