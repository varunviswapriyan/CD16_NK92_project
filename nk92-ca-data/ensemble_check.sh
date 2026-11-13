#!/usr/bin/env bash
# ensemble_check.sh -- settle Saeed's question: is the pZAP70 time course that
# drives the Ca model an average over NFsim replicates, or one realization?
#
#   bash ensemble_check.sh
#
# Read-only. The repo has mixed_runs/analysis10, analysis11, ... which look like
# replicate runs. This finds out whether anything averages them.
set -uo pipefail

PY="${PY:-python3}"
REPO="${REPO:-$HOME/CD16_NK92_project}"

hdr() { printf '\n\033[1m===== %s =====\033[0m\n' "$*"; }
ok()  { printf '   \033[32m+\033[0m %s\n' "$*"; }
no()  { printf '   \033[31m-\033[0m %s\n' "$*"; }

[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh 2>/dev/null
module load Miniconda3/4.9.2 >/dev/null 2>&1
CONDA_SH=/gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh
[ -f "$CONDA_SH" ] && { . "$CONDA_SH" >/dev/null 2>&1; conda activate CD16_v2 >/dev/null 2>&1; }
if [ -z "${PZ_PY:-}" ]; then
  for c in "${CONDA_PREFIX:+$CONDA_PREFIX/bin/python}" python python3; do
    [ -n "$c" ] || continue
    q=$(command -v "$c" 2>/dev/null || { [ -x "$c" ] && echo "$c"; })
    [ -n "$q" ] || continue
    "$q" -c "import numpy,pandas" >/dev/null 2>&1 && { PY="$q"; break; }
  done
else
  PY="$PZ_PY"
fi

# find the directory that holds the *_runs/ trees
BASE=""
while IFS= read -r d; do BASE="$d"; break; done < <(
  find "$REPO" -maxdepth 6 -type d -name '*_runs' 2>/dev/null | head -1)
[ -n "$BASE" ] && BASE="$(dirname "$BASE")"
if [ -z "$BASE" ]; then
  no "no *_runs directory under $REPO"
  no "try:  find $REPO -name 'model.bngl' | head"
  exit 1
fi
ok "model directory: $BASE"

###############################################################################
hdr "1. how many replicate runs per condition?"
###############################################################################
for r in "$BASE"/*_runs; do
  [ -d "$r" ] || continue
  n=$(find "$r" -maxdepth 1 -type d -name 'analysis*' 2>/dev/null | wc -l)
  echo "   $(basename "$r"): $n analysis* director$([ "$n" = 1 ] && echo y || echo ies)"
  find "$r" -maxdepth 1 -type d -name 'analysis*' 2>/dev/null | sort | head -6 \
    | sed "s#$r/#     #"
done

###############################################################################
hdr "2. what is inside one of them?"
###############################################################################
ONE="$(find "$BASE" -maxdepth 2 -type d -name 'analysis1*' 2>/dev/null | sort | head -1)"
if [ -n "$ONE" ]; then
  ok "$ONE"
  ls -la "$ONE" | sed 's/^/     /' | head -15
else
  no "no analysis* directory found"
fi

###############################################################################
hdr "3. does any script average over them?"
###############################################################################
# This is the question. Look for a glob over analysis*/ plus a mean().
FOUND=0
for f in "$BASE"/*.py; do
  [ -f "$f" ] || continue
  if grep -qE 'analysis|_runs|glob' "$f" 2>/dev/null; then
    echo "   [$(basename "$f")]"
    grep -nE 'analysis|_runs|glob|\.mean\(|np\.mean|average|sum\(.*\)/|/ *len\(' "$f" \
      | head -14 | sed 's/^/     /'
    FOUND=1
  fi
done
[ "$FOUND" = 1 ] || no "no python here mentions analysis*/ or glob"

###############################################################################
hdr "4. the decisive test: average the runs and compare to the csv"
###############################################################################
PZCSV="$(find "$HOME" -maxdepth 4 -name 'model_output_pzap_v77.csv' 2>/dev/null | head -1)"
if [ -z "$PZCSV" ]; then
  no "model_output_pzap_v77.csv not found -- skipping"
else
  ok "comparing against $PZCSV"
  "$PY" - "$BASE" "$PZCSV" <<'PYCHK' 2>&1 | sed 's/^/   /'
import sys, glob, os
import numpy as np, pandas as pd

base, csvp = sys.argv[1], sys.argv[2]
out = pd.read_csv(csvp)
print('csv columns:', list(out.columns))

def read_gdat(p):
    with open(p) as fh:
        head = fh.readline().lstrip('#').split()
    a = np.loadtxt(p)
    return head, a

for cond, tag in (('mixed','hetero'), ('gamma','gamma'), ('zeta','zeta')):
    runs = sorted(glob.glob(os.path.join(base, f'{cond}_runs', 'analysis*')))
    if not runs:
        continue
    gd = []
    for r in runs:
        hits = glob.glob(os.path.join(r, '*.gdat'))
        if hits:
            gd.append(hits[0])
    print(f'\n{cond}: {len(runs)} analysis dirs, {len(gd)} with a .gdat')
    if len(gd) < 2:
        print('  fewer than 2 trajectory files -- cannot test averaging here')
        continue
    try:
        head, a0 = read_gdat(gd[0])
    except Exception as e:
        print('  could not read', gd[0], e); continue
    # find a pZAP-looking column
    j = next((i for i, h in enumerate(head) if 'zap' in h.lower() or 'pz' in h.lower()), None)
    if j is None:
        print('  columns:', head[:12]); print('  no obvious pZAP column'); continue
    print(f'  using column {head[j]!r}')
    stack = []
    for p in gd:
        try:
            _, a = read_gdat(p)
            if a.shape[0] == a0.shape[0]:
                stack.append(a[:, j])
        except Exception:
            pass
    if len(stack) < 2:
        print('  trajectories differ in length -- cannot stack'); continue
    S = np.array(stack)
    t = a0[:, 0]
    mean = S.mean(0)
    sd = S.std(0)
    print(f'  spread across runs: mean sd/mean = {np.nanmean(sd/np.maximum(mean,1e-12)):.3f}')
    col = next((c for c in out.columns if tag in c and 'mean' in c), None)
    if col is None:
        print(f'  no csv column for {tag}'); continue
    y = np.interp(t, out[out.columns[0]], out[col])
    d_avg = np.nanmax(np.abs(mean - y)) / max(np.nanmax(np.abs(y)), 1e-12)
    d_one = np.nanmax(np.abs(S[0]  - y)) / max(np.nanmax(np.abs(y)), 1e-12)
    print(f'  csv vs AVERAGE of runs : max rel diff {d_avg:.4g}')
    print(f'  csv vs SINGLE run      : max rel diff {d_one:.4g}')
    print('  ->', 'csv IS the ensemble average' if d_avg < d_one * 0.5
          else ('csv looks like ONE run' if d_one < d_avg * 0.5 else 'inconclusive'))
PYCHK
fi

hdr "what this means"
cat <<'EOF'
   If section 4 says "csv IS the ensemble average", then the pZAP input is
   already averaged over replicates -- tell Saeed the answer on the call was
   wrong, and no change is needed.

   If it says "csv looks like ONE run", he is right: we should average the
   runs before fitting the Ca model, and that is a real improvement.

   If sections 1-3 find the analysis* dirs but section 4 cannot read them,
   paste this output and we will work out the file format.
EOF
