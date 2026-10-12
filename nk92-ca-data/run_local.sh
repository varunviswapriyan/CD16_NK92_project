#!/usr/bin/env bash
# run_local.sh -- Oscar's FCS zips -> fit-ready Ca CSVs, on your Mac.
#
#   bash run_local.sh ~/Downloads/Exported_FCS_Files_20260407NFH.zip \
#                     ~/Downloads/exported_NK_FCS_files_20260518NFH_2.zip
#
# Produces, under ./ca_out/ :
#   Ca_NK92_exp1_20260407.csv        drop-in for ca_data/Ca_NK92.csv (Indrani's schema)
#   Ca_NK92_exp2_20260518.csv        same, Experiment 2
#   *_full.csv                       adds n_* and the parental tube, for QC
#   ca_data.png                      mean +/- SE, both experiments
#
# Grid and columns match Ca_flow_analysis_mean_se.ipynb: 250 bins over
# 0..max(time) across all four tubes; columns time_seconds, mean_{zeta,gamma,
# hetero}, SE_{...}. --trim drops bins where a tube has already stopped.
set -euo pipefail

PY="${PY:-python3}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT="${OUT:-$PWD/ca_out}"
RAW="$OUT/raw"
NBINS="${NBINS:-250}"

if [ $# -lt 1 ]; then sed -n '2,18p' "$HERE/run_local.sh"; exit 1; fi

$PY -c 'import numpy, pandas, matplotlib' 2>/dev/null || {
  echo "installing numpy/pandas/matplotlib ..."
  $PY -m pip install --quiet --user numpy pandas matplotlib
}

mkdir -p "$RAW"
csvs=()
for z in "$@"; do
  base="$(basename "$z" .zip)"
  d="$RAW/$base"; rm -rf "$d"; mkdir -p "$d"
  unzip -q -o "$z" -d "$d"

  case "$base" in
    *20260407*|*0407*) name=Ca_NK92_exp1_20260407 ;;
    *20260518*|*0518*) name=Ca_NK92_exp2_20260518 ;;
    *)                 name="Ca_NK92_$base" ;;
  esac

  echo; echo "=== $base -> $name.csv"
  # drop-in file: her grid, her columns
  $PY "$HERE/parse_ca_fcs.py" "$d" -o "$OUT/$name.csv" \
      --grid indrani --nbins "$NBINS" --schema indrani --trim
  # QC file: cell counts + parental tube
  $PY "$HERE/parse_ca_fcs.py" "$d" -o "$OUT/${name}_full.csv" \
      --grid indrani --nbins "$NBINS" --schema full --trim >/dev/null
  csvs+=("$OUT/$name.csv")
done

echo
$PY "$HERE/plot_ca_data.py" "${csvs[@]}" -o "$OUT/ca_data.png"
echo; echo "done. files in $OUT:"; ls -la "$OUT"
