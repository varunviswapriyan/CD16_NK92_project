#!/bin/bash
# =============================================================================
# indrani_q.sh -- answers Indrani's questions 3 and 5 with results, not opinion.
#
#   bash ~/CD16_NK92_project/filesCC/indrani_q.sh submit    # Q3: the Ca ladder
#   bash ~/CD16_NK92_project/filesCC/indrani_q.sh flowkit   # Q5: find her code
#   bash ~/CD16_NK92_project/filesCC/indrani_q.sh flowprep  # Q5: install + read it
#   bash ~/CD16_NK92_project/filesCC/indrani_q.sh report    # Q3: the comparison
#   bash ~/CD16_NK92_project/filesCC/indrani_q.sh status
#
# Q3 -- "Did he change the Ca model?  I see a new parameter, S, being
#        estimated, but I do not see K3 and K4."
#
# Yes on all three counts, and this submits the ladder that shows whether any
# of it was necessary.  Four models, same data, same optimiser:
#
#   orig    C1 C2 g k3 k4        HER MODEL EXACTLY.  S = 1, both drive terms.
#                                Never run before -- even 'full' had S in it.
#   orig1   C1 C2 g k3           same, plus k4 = 1 to close the C1*k4 valley.
#   full1   C1 C2 g k3 S         her structure, S restored, k4 = 1.
#   lin1    C1 C2 g S            what we have been reporting.
#
# READ IT LIKE THIS
#   orig fits as well as lin1    -> S was unnecessary.  Report orig, and we
#                                   drop our addition entirely.
#   orig fits badly              -> S is a units conversion, not a free knob:
#                                   her k1 = k2 = 0.7 are micromolar while the
#                                   data is fluorescence 50-200, so with S = 1
#                                   the feedback term saturates at 1 and the
#                                   h-gate target collapses to 0.  Two of her
#                                   own mechanisms switch off.  That is the
#                                   justification for S, and the SSR is the
#                                   evidence.
#   full1 matches lin1 on SSR    -> k3 can stay in the model; only k4 is pinned,
#      with tight CIs               and nothing of hers is removed.
#
# Nothing here is a judgement call once the numbers are in.
# =============================================================================
set -o pipefail
MODE=${1:-submit}
F=~/CD16_NK92_project/filesCC
CADIR=$HOME/Ca_fit_c02
NCA=${NCA:-50}
VARIANTS=${VARIANTS:-"orig orig1 full1"}

echo "setting up the job environment..."
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
    "$q" -c "import numpy,pandas" >/dev/null 2>&1 && { PY="$q"; break; }
  done
fi
[ -n "$PY" ] || { echo "ERROR: no python with numpy/pandas"; exit 1; }
echo "using $PY"
echo

case "$MODE" in submit|flowkit|flowprep|report|status) ;; *)
  echo "use: submit | flowkit | flowprep | report | status"; exit 2 ;; esac

# ----------------------------------------------------------------- flowprep --
# Get ready for Oscar's .fcs files: make sure FlowKit is importable, then show
# what Indrani's notebook actually does so the new files can be run through the
# same pipeline rather than a reimplementation of it.
if [ "$MODE" = flowprep ]; then
  echo "=============== 1. FlowKit ==============="
  if "$PY" -c "import flowkit" >/dev/null 2>&1; then
    "$PY" -c "import flowkit; print('  already installed:', flowkit.__version__)"
  else
    echo "  not installed -- installing into $(dirname "$(dirname "$PY")")"
    "$PY" -m pip install --quiet flowkit 2>&1 | tail -5
    "$PY" -c "import flowkit; print('  now at', flowkit.__version__)" 2>/dev/null \
      || { echo "  INSTALL FAILED."; echo "  The cluster may block PyPI.  Try:"; \
           echo "    $PY -m pip install --user flowkit"; }
  fi

  NB=$(find "$HOME" -maxdepth 5 -name 'Ca_flow_analysis_mean_se.ipynb' 2>/dev/null | head -1)
  echo
  echo "=============== 2. Indrani's notebook ==============="
  if [ -z "$NB" ]; then echo "  not found"; exit 0; fi
  echo "  $NB"
  echo
  NB="$NB" "$PY" - <<'PYNB'
import json, os
nb = json.load(open(os.environ['NB']))
cells = [c for c in nb.get('cells', []) if c.get('cell_type') == 'code']
print('  %d code cells.  Source only, outputs stripped:' % len(cells))
print('  ' + '-' * 68)
for i, c in enumerate(cells, 1):
    src = ''.join(c.get('source', [])).rstrip()
    if not src.strip():
        continue
    print('\n  ----- cell %d -----' % i)
    for ln in src.splitlines():
        print('  ' + ln)
PYNB
  echo
  echo "=============== 3. .fcs files you already have ==============="
  find "$HOME" -maxdepth 5 -name '*.fcs' 2>/dev/null | xargs -r -n1 basename \
    | sort -u | sed 's|^|  |'
  echo
  echo "  Indrani's email says you still need, from Oscar:"
  echo "    Experiment 1, Replicate 2   (truncate at 280 s)"
  echo "    Experiment 2, Replicate 1   (fit separately, supplementary)"
  echo "  With Exp1 Rep2 added you have THREE replicates, which means the Ca"
  echo "  CIs can use the same day-resampling bootstrap as pZAP instead of the"
  echo "  parametric draw from mean +/- SE."
  exit 0
fi

# ------------------------------------------------------------------ flowkit --
if [ "$MODE" = flowkit ]; then
  echo "=============== Q5: the FlowKit code Indrani mentioned ==============="
  echo "she said: 'The code for that is shared in NCH cluster under ca_data folder.'"
  echo
  for d in "$CADIR/ca_data" /home/gddaslab/share/Varun_Indrani/ca_data \
           /home/gddaslab/share/ca_data; do
    [ -d "$d" ] || continue
    echo "--- $d"
    ls -la "$d" | sed 's|^|    |'
    echo
  done
  echo "--- anything else under a ca_data folder:"
  find "$HOME" /home/gddaslab/share -maxdepth 6 -ipath '*ca_data*' \
       \( -name '*.py' -o -name '*.ipynb' -o -name '*.fcs' \) 2>/dev/null \
    | head -30 | sed 's|^|    |'
  echo
  echo "--- is FlowKit installed in the env?"
  "$PY" -c "import flowkit; print('    flowkit', flowkit.__version__)" 2>/dev/null \
    || echo "    NOT installed -- pip install flowkit"
  exit 0
fi

# ------------------------------------------------------------------- status --
if [ "$MODE" = status ]; then
  for v in $VARIANTS lin1; do
    n=$(ls "$CADIR/out_boot_ca2/$v"/boot_*.json 2>/dev/null | wc -l)
    printf '  %-8s %s samples written\n' "$v" "$n"
  done
  echo; echo "queue: $(squeue -u "$USER" -h 2>/dev/null | wc -l) jobs"
  exit 0
fi

# ------------------------------------------------------------------- report --
if [ "$MODE" = report ]; then
  [ -d "$CADIR" ] || { echo "$CADIR not found"; exit 1; }
  cd "$CADIR" && bash "$F/ca_check.sh" 2>/dev/null || "$PY" "$F/bootstrap_ca2.py" report
  cat <<'EOF'

WHAT THE ANSWER TO INDRANI IS, BY OUTCOME
  orig  ref SSR ~= 5.21e4  -> S was not needed.  Drop it, report her model.
  orig  ref SSR much worse -> S is a units conversion; say so and give both SSRs.
  full1 ref SSR ~= 5.21e4
        and median w < 1   -> report full1: k3 kept, only k4 pinned, S retained
                              as a calibration constant.  Nothing of hers removed.
EOF
  exit 0
fi

# ------------------------------------------------------------------- submit --
[ -d "$CADIR" ] || { echo "ERROR: $CADIR not found"; exit 1; }
cd "$CADIR" || exit 1
echo "=============== Q3: submitting the Ca model ladder ==============="
for v in $VARIANTS; do
  if ! "$PY" -c "
import sys,re
src=open('$F/bootstrap_ca2.py').read()
ns={'N01':{'k3':12.81,'k4':0.1039,'S':31.6}}
exec(re.search(r'VARIANTS = \{.*?\n\}', src, re.S).group(0), ns)
sys.exit(0 if '$v' in ns['VARIANTS'] else 1)"; then
    echo "  '$v' is not in the deployed bootstrap_ca2.py -- sync/pull first, skipping"
    continue
  fi
  echo "--- $v  ($((NCA+1)) jobs)"
  "$PY" "$F/bootstrap_ca2.py" submit "$NCA" "$v" | tail -2
done
echo
echo "later:  bash \$0 status      then      bash \$0 report"
