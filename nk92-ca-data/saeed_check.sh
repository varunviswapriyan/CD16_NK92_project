#!/usr/bin/env bash
# saeed_check.sh -- answer the questions Saeed raised on the 13 Nov call.
#
#   bash saeed_check.sh
#
# Read-only. Changes nothing. Four sections:
#   A  is the pZAP70 input averaged over NFsim replicates, or one realization?
#   B  how slow is the Ca solve loop, really?
#   C  the solve/rhs code, to show him
#   D  what is and isn't in the GitHub repo
set -uo pipefail

PY="${PY:-python3}"
REPO="${REPO:-$HOME/CD16_NK92_project}"
BOOT="${BOOT:-$REPO/filesCC/bootstrap_ca2.py}"
PZFIT="${PZFIT:-$HOME/NK92_fit_v77/estimate_params_pzap_cleaned_up}"

hdr() { printf '\n\033[1m========== %s ==========\033[0m\n' "$*"; }
sub() { printf '\n\033[1m-- %s\033[0m\n' "$*"; }
ok()  { printf '   \033[32m+\033[0m %s\n' "$*"; }
no()  { printf '   \033[31m-\033[0m %s\n' "$*"; }

# ------------------------------------------------------------------ python env
[ -f /etc/profile.d/modules.sh ] && . /etc/profile.d/modules.sh 2>/dev/null
module load Miniconda3/4.9.2 >/dev/null 2>&1
CONDA_SH=/gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh
[ -f "$CONDA_SH" ] && { . "$CONDA_SH" >/dev/null 2>&1; conda activate CD16_v2 >/dev/null 2>&1; }
if [ -z "${PZ_PY:-}" ]; then
  for c in "${CONDA_PREFIX:+$CONDA_PREFIX/bin/python}" python python3; do
    [ -n "$c" ] || continue
    q=$(command -v "$c" 2>/dev/null || { [ -x "$c" ] && echo "$c"; })
    [ -n "$q" ] || continue
    "$q" -c "import numpy,pandas,scipy" >/dev/null 2>&1 && { PY="$q"; break; }
  done
else
  PY="$PZ_PY"
fi
echo "python: $PY"

###############################################################################
hdr "A.  Is the pZAP70 input an ensemble average or one run?"
###############################################################################
# Saeed: "are you averaging over some ensemble or is it just one random
# configuration?"  The columns are named mean_pZAP_*, so something is averaged
# -- the question is whether it is over NFsim replicates or over cells in one.

if [ ! -d "$PZFIT" ]; then
  no "pZAP fit directory not at $PZFIT  (set PZFIT=<dir>)"
else
  ok "pZAP fit directory: $PZFIT"

  sub "simulate() calls in the .bngl -- look for n_runs / seed / method=>nf"
  for f in "$PZFIT"/*.bngl; do
    [ -f "$f" ] || continue
    echo "   [$(basename "$f")]"
    grep -nE 'simulate|method *=>|seed|n_steps|t_end|suffix' "$f" | sed 's/^/     /' || echo "     (none)"
  done

  sub "NFsim trajectory files -- MANY per condition means replicates were run"
  for ext in gdat cdat; do
    n=$(find "$PZFIT" -maxdepth 3 -name "*.$ext" 2>/dev/null | wc -l)
    echo "   *.$ext : $n file(s)"
    find "$PZFIT" -maxdepth 3 -name "*.$ext" 2>/dev/null | head -5 | sed 's/^/     /'
  done

  sub "anything that looks like a per-run loop"
  grep -rnE 'for .*(seed|run|rep|i in)|--seed|-seed|nf .*-o' \
       "$PZFIT"/*.sh "$PZFIT"/*.py 2>/dev/null | head -12 | sed 's/^/     /' || echo "     (no .sh/.py found)"
fi

sub "the file the Ca model actually reads"
PZCSV="$(find "$HOME" -maxdepth 4 -name 'model_output_pzap_v77.csv' 2>/dev/null | head -1)"
if [ -z "$PZCSV" ]; then
  no "model_output_pzap_v77.csv not found"
else
  ok "$PZCSV"
  echo "   columns:"
  head -1 "$PZCSV" | tr ',' '\n' | sed 's/^/     /'
  echo "   rows: $(( $(wc -l < "$PZCSV") - 1 ))"
  echo
  echo "   READ THIS: a column named sd_/std_/se_ beside each mean_ would prove"
  echo "   averaging over replicates. Only mean_ columns = cannot tell from the"
  echo "   file alone; check the .bngl simulate() call above."
fi

###############################################################################
hdr "B.  How slow is the Ca solve loop?"
###############################################################################
# Saeed: "there is no problem with the code ... but I think it is a bit slow."
# He is right that solve() calls odeint once per time step in a Python loop.
# This measures it so the conversation has a number in it.

CADIR=""; best=0
while IFS= read -r d; do
  c=$(ls "$d/out_boot_ca2/orig"/boot_*.json 2>/dev/null | wc -l)
  [ "$c" -gt "$best" ] && { best="$c"; CADIR="$d"; }
done < <(find "$HOME" -maxdepth 2 -type d -name 'Ca_fit*' 2>/dev/null | sort)

if [ -z "$CADIR" ] || [ ! -f "$BOOT" ]; then
  no "need a Ca fit dir and $BOOT -- skipping"
else
  ok "timing in $CADIR"
  ( cd "$CADIR" && "$PY" - "$BOOT" <<'PYTIME'
import sys, runpy, time, json, os
import numpy as np
g = runpy.run_path(sys.argv[1], run_name='_timing')
cfg = g['VARIANTS']['orig']
jf = os.path.join(g.get('ROOT','out_boot_ca2'), 'orig', 'boot_000.json')
d  = json.load(open(jf))
p  = dict(g['BASEFIX']); p.update(g['N01']); p.update(cfg['fix'])
p.update(d.get('fixed', {})); p.update(d.get('params', {}))
data = g['resample'](g['load_base'](), seed=1, jitter=False)
for fn in ('solve','rhs'):
    g[fn].__globals__.update({'CFG':cfg,'PN':cfg['fit'],'DATA':data})

npts = {k: len(v['t']) for k, v in data.items()}
print('   time points per condition:', npts)

t0 = time.perf_counter()
N = 20
for _ in range(N):
    for k in ('zeta','gamma','hetero'):
        g['solve'](data[k], p, cfg)
dt = (time.perf_counter() - t0) / N
steps = sum(npts.values())
print(f'   one full cost evaluation (3 conditions): {dt*1000:.1f} ms')
print(f'   that is {steps} odeint calls -> {dt/steps*1e6:.0f} us per call')
print()
print('   A PSO fit is roughly particles x iterations x starts cost evaluations.')
print(f'   At 30 particles x 200 iterations x 4 starts that is ~{30*200*4*dt/60:.0f} min per sample,')
print(f'   and a 51-sample bootstrap is ~{51*30*200*4*dt/3600:.1f} h of CPU.')
print()
print('   Saeed\'s suggestion: one solve_ivp(RK45) call over the whole span with')
print('   the pZAP interpolator passed in, instead of a Python loop of odeint')
print('   calls. Same equations, same answers, fewer Python-level round trips.')
PYTIME
  ) 2>&1 | grep -v Warning
fi

###############################################################################
hdr "C.  The code he looked at"
###############################################################################
if [ -f "$BOOT" ]; then
  sub "rhs() -- the model (Indrani's two equations)"
  grep -n -A 12 "^def rhs" "$BOOT" | sed 's/^/   /'
  sub "solve() -- the integration loop he called slow"
  grep -n -A 12 "^def solve" "$BOOT" | sed 's/^/   /'
else
  no "$BOOT not found"
fi

###############################################################################
hdr "D.  What is in the GitHub repo (what Saeed can see once added)"
###############################################################################
if [ -d "$REPO/.git" ]; then
  ( cd "$REPO" && echo "   remote: $(git remote get-url origin 2>/dev/null)" \
    && echo "   tracked files:" && git ls-files | sed 's/^/     /' | head -40 )
  sub "is the BioNetGen model in the repo?"
  if ( cd "$REPO" && git ls-files | grep -q '\.bngl$' ); then
    ok "yes"
  else
    no "NO -- he asked 'where is the full model?'. To add it:"
    echo "     mkdir -p $REPO/bngl && cp $PZFIT/*.bngl $REPO/bngl/"
    echo "     cd $REPO && git add bngl && git commit -m 'add BioNetGen models' && git push"
  fi
else
  no "$REPO is not a git repo"
fi

hdr "summary of what to tell him"
cat <<'EOF'
   1. Ensemble: see section A. If the .bngl simulate() has no replicate loop and
      there is one trajectory file per condition, he is right and we average.
      If mean_ already spans replicates, correct what was said on the call.
   2. Speed: section B has the measured number. His solve_ivp suggestion is
      sound; it changes nothing about the results.
   3. Interpolation: already correct -- he confirmed this on screen.
   4. Repo: add drsamalik as a collaborator, and push the .bngl files.
EOF
