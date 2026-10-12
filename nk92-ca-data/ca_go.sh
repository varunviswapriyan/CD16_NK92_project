#!/usr/bin/env bash
# ca_go.sh -- one command for the whole Ca cross-check, on the cluster.
#
#   bash ca_go.sh              find everything, parse both experiments,
#                              run Indrani's notebook, diff everything, print a verdict
#   bash ca_go.sh swap         back up the old Ca_NK92.csv and put the new one in
#   bash ca_go.sh serve        start http.server so you can open the PNGs
#
# Nothing here overwrites anything unless you type `swap`.
set -uo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
WORK="${WORK:-$HOME/ca_check}"
PY="${PY:-python3}"
NBINS="${NBINS:-250}"
CMD="${1:-check}"

say()  { printf '\n\033[1m== %s\033[0m\n' "$*"; }
ok()   { printf '   \033[32mOK\033[0m   %s\n' "$*"; }
warn() { printf '   \033[33mNOTE\033[0m %s\n' "$*"; }
bad()  { printf '   \033[31mFAIL\033[0m %s\n' "$*"; }

# ---------------------------------------------------------------- serve
if [ "$CMD" = serve ]; then
  ip="$(hostname -I 2>/dev/null | awk '{print $1}')"
  cd "$HOME" || exit 1
  pkill -f 'http.server 8000' 2>/dev/null
  nohup "$PY" -m http.server 8000 >/dev/null 2>&1 &
  sleep 1
  echo "serving $HOME on port 8000"
  echo "  http://${ip:-10.73.170.129}:8000/ca_check/"
  exit 0
fi

# ---------------------------------------------------------------- discovery
say "finding things"

# the repo must have actually arrived (see the gitlink problem in the README)
if [ ! -f "$HERE/parse_ca_fcs.py" ]; then
  bad "parse_ca_fcs.py is not next to this script -- the repo did not transfer."
  bad "On your Mac: rm -rf ~/Downloads/indrani_code/nk92-ca-data/.git, then re-add and push."
  exit 1
fi
ok "repo: $HERE"

# the FCS zips
ZIPS=()
while IFS= read -r line; do [ -n "$line" ] && ZIPS+=("$line"); done < <(
  { find "$HERE/.." "$HOME" -maxdepth 4 -type f -iname '*.zip' 2>/dev/null \
      | grep -iE 'fcs.*(20260407|20260518)' \
      | while IFS= read -r f; do readlink -f "$f"; done; } | sort -u)
if [ "${#ZIPS[@]}" -eq 0 ]; then
  bad "no FCS zips found. Push them: cp ~/Downloads/*FCS*.zip ~/Downloads/indrani_code/nk92-ca-data/fcs_zips/"
  exit 1
fi
for z in "${ZIPS[@]}"; do ok "zip: $z"; done

# Indrani's original notebook
NB_ORIG="$(find "$HOME" -maxdepth 6 -name 'Ca_flow_analysis_mean_se.ipynb' 2>/dev/null | head -1)"
[ -n "$NB_ORIG" ] && ok "her notebook: $NB_ORIG" || warn "her original notebook not found (not fatal)"

NB_FIX="$HERE/Ca_flow_analysis_mean_se_fixed.ipynb"
[ -f "$NB_FIX" ] && ok "fixed notebook: $NB_FIX" || warn "fixed notebook missing"

# every Ca_NK92.csv on the system -- one of these is what the fits actually read
OLDCSV=()
while IFS= read -r line; do [ -n "$line" ] && OLDCSV+=("$line"); done < <(
  find "$HOME" -maxdepth 6 -name 'Ca_NK92.csv' -not -path "$WORK/*" 2>/dev/null | sort -u)
if [ "${#OLDCSV[@]}" -eq 0 ]; then
  warn "no existing Ca_NK92.csv found -- nothing to compare the new data against"
else
  for c in "${OLDCSV[@]}"; do ok "existing fit input: $c"; done
fi

has_flowkit=0
"$PY" -c 'import flowkit' 2>/dev/null && { has_flowkit=1; ok "flowkit available"; } \
  || warn "flowkit not importable -- will skip the notebook run, CLI parse still works"

# ---------------------------------------------------------------- parse
mkdir -p "$WORK"
NEW=()
for z in "${ZIPS[@]}"; do
  base="$(basename "$z" .zip)"
  case "$base" in
    *20260407*) tag=exp1_20260407 ;;
    *20260518*) tag=exp2_20260518 ;;
    *)          tag="$(echo "$base" | tr -cd 'A-Za-z0-9_')" ;;
  esac
  d="$WORK/$tag"
  say "parsing $tag"
  rm -rf "$d"; mkdir -p "$d"
  unzip -q -o "$z" -d "$d" || { bad "unzip failed: $z"; continue; }
  # flatten: put the four .fcs directly in $d so the notebook finds them too
  find "$d" -name '__MACOSX' -type d -exec rm -rf {} + 2>/dev/null
  find "$d" -mindepth 2 -name '*.fcs' ! -name '._*' -exec mv {} "$d"/ \; 2>/dev/null
  find "$d" -mindepth 1 -type d -empty -delete 2>/dev/null

  out="$WORK/Ca_NK92_$tag.csv"
  if "$PY" "$HERE/parse_ca_fcs.py" "$d" -o "$out" \
        --grid indrani --nbins "$NBINS" --schema indrani --trim; then
    NEW+=("$out"); ok "wrote $out"
  else
    bad "parse failed for $tag"; continue
  fi

  # --- her notebook, same folder, then diff
  if [ "$has_flowkit" = 1 ] && [ -f "$NB_FIX" ]; then
    cp "$NB_FIX" "$d/"
    ( cd "$d" && "$PY" -m jupyter nbconvert --to notebook --execute \
        Ca_flow_analysis_mean_se_fixed.ipynb --output executed.ipynb \
        --ExecutePreprocessor.timeout=1800 ) >"$WORK/${tag}_nb.log" 2>&1
    if [ -f "$d/Ca_NK92.csv" ]; then
      ok "notebook produced $d/Ca_NK92.csv"
      echo "   --- notebook vs CLI parser ---"
      "$PY" "$HERE/compare_ca_csv.py" "$d/Ca_NK92.csv" "$out" | sed 's/^/   /'
    else
      bad "notebook did not produce a CSV; see $WORK/${tag}_nb.log"
      tail -15 "$WORK/${tag}_nb.log" | sed 's/^/   | /'
    fi
  fi
done

# ---------------------------------------------------------------- the gate
EXP1="$WORK/Ca_NK92_exp1_20260407.csv"
if [ "${#OLDCSV[@]}" -gt 0 ] && [ -f "$EXP1" ]; then
  say "THE GATE -- new Exp 1 vs whatever the current fits were run on"
  for c in "${OLDCSV[@]}"; do
    echo "   --- $c ---"
    "$PY" "$HERE/compare_ca_csv.py" "$c" "$EXP1" | sed 's/^/   /'
  done
fi

# ---------------------------------------------------------------- picture
if [ "${#NEW[@]}" -gt 0 ]; then
  "$PY" "$HERE/plot_ca_data.py" "${NEW[@]}" -o "$WORK/ca_data.png" >/dev/null 2>&1 \
    && ok "plot: $WORK/ca_data.png   (bash ca_go.sh serve, then open /ca_check/ca_data.png)"
fi

say "summary"
echo "   new CSVs in $WORK:"
ls -1 "$WORK"/Ca_NK92_*.csv 2>/dev/null | sed 's/^/     /'
cat <<EOF

   If every diff above said IDENTICAL, the old SSRs still stand.
   If the GATE said DIFFERENT, the SSRs in the draft email are stale -- refit first.

   To put the new data in and submit:   bash $HERE/ca_go.sh swap
EOF

# ---------------------------------------------------------------- swap
if [ "$CMD" = swap ]; then
  say "swap"
  [ -f "$EXP1" ] || { bad "no $EXP1 -- run 'bash ca_go.sh' first"; exit 1; }
  if [ "${#OLDCSV[@]}" -eq 0 ]; then bad "no target Ca_NK92.csv found"; exit 1; fi
  if [ "${#OLDCSV[@]}" -gt 1 ]; then
    echo "   more than one Ca_NK92.csv exists:"
    printf '     %s\n' "${OLDCSV[@]}"
    echo "   set the one you want:  TARGET=<path> bash ca_go.sh swap"
  fi
  TARGET="${TARGET:-${OLDCSV[0]}}"
  cp -v "$TARGET" "$TARGET.bak.$(date +%F-%H%M)"
  cp -v "$EXP1" "$TARGET"
  FITDIR="$(dirname "$(dirname "$TARGET")")"
  ok "swapped. now submit with:"
  echo "       cd $FITDIR && bash ~/CD16_NK92_project/filesCC/finish.sh ca"
  echo "   (finish.sh ca -- NOT finish.sh submit, that deletes the 60 grid cells)"
fi
