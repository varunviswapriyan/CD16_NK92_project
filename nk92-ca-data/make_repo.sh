#!/usr/bin/env bash
set -uo pipefail

DEST="${DEST:-$HOME/nk92-cd16-model}"
SRC_REPO="${SRC_REPO:-$HOME/CD16_NK92_project}"
PZFIT="${PZFIT:-$HOME/NK92_fit_v77/estimate_params_pzap_cleaned_up}"
PY="${PY:-python3}"

[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh 2>/dev/null
module load Miniconda3/4.9.2 >/dev/null 2>&1
CONDA_SH=/gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh
[ -f "$CONDA_SH" ] && { . "$CONDA_SH" >/dev/null 2>&1; conda activate CD16_v2 >/dev/null 2>&1; }
command -v python >/dev/null 2>&1 && PY=python

say() { printf '\n\033[1m%s\033[0m\n' "$*"; }
got() { printf '   %s\n' "$*"; }
miss() { printf '   \033[31mmissing\033[0m %s\n' "$*"; }

rm -rf "$DEST"
mkdir -p "$DEST"/{model,pzap,calcium,data,figures}

grab() {
  if [ -f "$1" ]; then cp "$1" "$2"; got "$(basename "$1") -> ${2#$DEST/}"
  else miss "$1"; fi
}

say "bngl models"
for f in "$PZFIT"/*.bngl; do [ -f "$f" ] && grab "$f" "$DEST/model/"; done
if ! ls "$DEST"/model/*.bngl >/dev/null 2>&1; then
  for f in $(find "$SRC_REPO" -name 'JJ_*.bngl' 2>/dev/null | head -5); do
    grab "$f" "$DEST/model/"
  done
fi

say "pzap"
grab "$SRC_REPO/filesCC/bootstrap_pzap.py" "$DEST/pzap/"
for n in calculate_SSR_zeta.py calculate_SSR_gamma.py calculate_SSR_mixed.py \
         calculate_pZAP.py interpolation.py estimate_k3_k4.py; do
  f=$(find "$SRC_REPO" -name "$n" 2>/dev/null | head -1)
  [ -n "$f" ] && grab "$f" "$DEST/pzap/"
done

say "calcium"
grab "$SRC_REPO/filesCC/bootstrap_ca2.py" "$DEST/calcium/"
f=$(find "$HOME" -maxdepth 5 -name 'ca_ODE.py' 2>/dev/null | head -1)
[ -n "$f" ] && grab "$f" "$DEST/calcium/"
f=$(find "$HOME" -maxdepth 5 -name 'estimate_ca_params.py' 2>/dev/null | head -1)
[ -n "$f" ] && grab "$f" "$DEST/calcium/"

say "data"
for f in "$HOME"/ca_check/Ca_NK92_exp*.csv; do [ -f "$f" ] && grab "$f" "$DEST/data/"; done
f=$(find "$HOME" -maxdepth 4 -name 'model_output_pzap_v77.csv' 2>/dev/null | head -1)
[ -n "$f" ] && grab "$f" "$DEST/data/"
f=$(find "$HOME" -maxdepth 6 -name 'pZAP70_Tyr493_mean.csv' 2>/dev/null | head -1)
[ -n "$f" ] && grab "$f" "$DEST/data/"

say "figures"
for n in clean_fig.py plot_pzap_fit.py parse_ca_fcs.py; do
  grab "$SRC_REPO/nk92-ca-data/$n" "$DEST/figures/"
done

say "stripping prose comments and docstrings from the helper scripts"
"$PY" - "$DEST" <<'PYSTRIP'
import os, sys, ast

root = sys.argv[1]
targets = ['figures/clean_fig.py', 'figures/plot_pzap_fit.py', 'figures/parse_ca_fcs.py']

def drop_docstrings(src):
    tree = ast.parse(src)
    for node in ast.walk(tree):
        if isinstance(node, (ast.Module, ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            if (node.body and isinstance(node.body[0], ast.Expr)
                    and isinstance(node.body[0].value, ast.Constant)
                    and isinstance(node.body[0].value.value, str)):
                node.body.pop(0)
                if not node.body:
                    node.body.append(ast.Pass())
    return ast.unparse(ast.fix_missing_locations(tree))

for rel in targets:
    p = os.path.join(root, rel)
    if not os.path.exists(p):
        continue
    src = open(p).read()
    try:
        new = drop_docstrings(src)
        compile(new, p, 'exec')
        open(p, 'w').write(new + "\n")
        print(f'   stripped {rel}')
    except Exception as e:
        print(f'   SKIPPED {rel}: {e}')
PYSTRIP

say "git"
cd "$DEST" || exit 1
printf '__pycache__/\n*.pyc\n*.png\n' > .gitignore
git init -q
git add -A
git -c user.name="Varun Viswapriyan" \
    -c user.email="varun.viswapriyan@nationwidechildrens.org" \
    commit -q -m "NK92 CD16 signalling model: BioNetGen pZAP70, calcium ODE, fitting"

say "contents"
git ls-files | sed 's/^/   /'

cat <<EOF

   Created $DEST

   To put it on GitHub:
     1. github.com/new  ->  name it  nk92-cd16-model  ->  Private  ->  no README
     2. cd $DEST
        git remote add origin git@github.com:varunviswapriyan/nk92-cd16-model.git
        git branch -M main
        git push -u origin main
     3. Settings -> Collaborators -> add  drsamalik

EOF
