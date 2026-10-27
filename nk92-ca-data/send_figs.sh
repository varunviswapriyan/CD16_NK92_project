#!/usr/bin/env bash
# send_figs.sh -- build the best Ca fit figures and push them, so you can pull
#                 them on your Mac and attach them to an email.
#
#   bash send_figs.sh
#
# Builds the 'orig' bootstrap band (Indrani's 5-parameter model) for Experiment 1
# and Experiment 2, drops both into the repo under figures/ca_fits/, commits and
# pushes. Prints the one command to run on your Mac.
set -uo pipefail

REPO="${REPO:-$HOME/CD16_NK92_project}"
F="${F:-$REPO/filesCC}"
BOOT="$F/bootstrap_ca2.py"
DEST="$REPO/figures/ca_fits"
VARIANT="${VARIANT:-orig}"

ok()   { printf '   \033[32m+\033[0m %s\n' "$*"; }
bad()  { printf '   \033[31m!\033[0m %s\n' "$*"; }
say()  { printf '\n\033[1m== %s\033[0m\n' "$*"; }

# ------------------------------------------------------------------ python env
[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh 2>/dev/null
module load Miniconda3/4.9.2 >/dev/null 2>&1
CONDA_SH=/gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh
[ -f "$CONDA_SH" ] && { . "$CONDA_SH" >/dev/null 2>&1; conda activate CD16_v2 >/dev/null 2>&1; }
PY="${PZ_PY:-}"
if [ -z "$PY" ]; then
  for c in "${CONDA_PREFIX:+$CONDA_PREFIX/bin/python}" python python3; do
    [ -n "$c" ] || continue
    q=$(command -v "$c" 2>/dev/null || { [ -x "$c" ] && echo "$c"; })
    [ -n "$q" ] || continue
    "$q" -c "import numpy,pandas,matplotlib" >/dev/null 2>&1 && { PY="$q"; break; }
  done
fi
[ -n "$PY" ] || { bad "no python with numpy/pandas/matplotlib"; exit 1; }
[ -f "$BOOT" ] || { bad "bootstrap_ca2.py not at $BOOT"; exit 1; }

mkdir -p "$DEST"

# ------------------------------------------------- find the two fit directories
# Experiment 1 = the directory whose ca_data/Ca_NK92.csv the original fits used.
# Experiment 2 = Ca_fit_exp2. Both are located, not assumed.
EXP1=""; best=0
while IFS= read -r d; do
  c=$(ls "$d/out_boot_ca2/$VARIANT"/boot_*.json 2>/dev/null | wc -l)
  [ "$c" -gt "$best" ] && { best="$c"; EXP1="$d"; }
done < <(find "$HOME" -maxdepth 2 -type d -name 'Ca_fit*' ! -name 'Ca_fit_exp2' 2>/dev/null | sort)
[ -n "$EXP1" ] && ok "Experiment 1 = $EXP1 ($best $VARIANT samples, the most of any Ca_fit* dir)"
EXP2="$HOME/Ca_fit_exp2"

say "building $VARIANT band plots"
n=0
for pair in "Exp1:$EXP1" "Exp2:$EXP2"; do
  tag="${pair%%:*}"; dir="${pair#*:}"
  [ -n "$dir" ] && [ -d "$dir/out_boot_ca2/$VARIANT" ] || { bad "$tag: no $VARIANT results"; continue; }
  ns=$(ls "$dir/out_boot_ca2/$VARIANT"/boot_*.json 2>/dev/null | wc -l)
  ( cd "$dir" && "$PY" "$BOOT" band "$VARIANT" >/dev/null 2>&1 )
  src="$dir/out_boot_ca2/band_$VARIANT.png"
  if [ -f "$src" ]; then
    cp "$src" "$DEST/Ca_fit_${tag}_${VARIANT}.png"
    ok "$tag -> figures/ca_fits/Ca_fit_${tag}_${VARIANT}.png  ($ns samples, from $dir)"
    n=$((n+1))
  else
    bad "$tag: band plot not produced"
  fi
done
[ "$n" -gt 0 ] || { bad "nothing built, stopping"; exit 1; }

# ------------------------------------------------------------- parameter table
( cd "${EXP1:-$HOME}" && "$PY" "$BOOT" report 2>/dev/null ) > "$DEST/params_Exp1.txt"
( cd "$EXP2" && "$PY" "$BOOT" report 2>/dev/null )          > "$DEST/params_Exp2.txt"
ok "parameter tables -> figures/ca_fits/params_Exp{1,2}.txt"

# --------------------------------------------------------------------- push it
say "pushing"
cd "$REPO" || exit 1
git add -f figures/ca_fits >/dev/null 2>&1
if git diff --cached --quiet; then
  ok "nothing changed since last push"
else
  git -c user.name="Varun Viswapriyan" -c user.email="varun.viswapriyan@nationwidechildrens.org" \
      commit -q -m "Ca fit figures ($VARIANT), Exp 1 and Exp 2" && ok "committed"
fi
PUSHLOG=$(git push 2>&1); PUSHRC=$?
echo "$PUSHLOG" | tail -3
if [ "$PUSHRC" -eq 0 ]; then
  ok "pushed"
  cat <<EOF

   On your Mac:

     cd ~/Downloads/indrani_code && git pull
     open nk92-ca-data/../figures/ca_fits   # or just: open figures/ca_fits

   Attach Ca_fit_Exp1_${VARIANT}.png and Ca_fit_Exp2_${VARIANT}.png

EOF
else
  bad "push failed (credentials?). The files are in $DEST --"
  bad "view them instead: bash $REPO/nk92-ca-data/plots.sh"
fi
