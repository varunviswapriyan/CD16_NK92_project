#!/bin/bash
# run_ca_v77.sh -- re-run the Ca fits against the v77 pZAP curves (no gamma
# affinity fraction), so the Ca stage matches the pZAP version being sent.
#
#   bash ~/CD16_NK92_project/filesCC/run_ca_v77.sh          # export + submit
#   bash ~/CD16_NK92_project/filesCC/run_ca_v77.sh report   # results table
#
# Step 1 exports v77's pZAP curves from its BNG .gdat outputs (~1 min).
# Step 2 submits the 5 shared-parameter Ca variants in PARALLEL (~15-25 min):
#   m13 per-ITAM + saturation  (best on v69 inputs: 6.01e4)
#   m01 per-ITAM, her 5 params (7.61e4)
#   m12 per-ITAM + linear, 4 params
#   m02 per-ITAM with FITTED exponent  <- adapts if v77 inputs are spaced differently
#   c21 fully shared, no normalization (baseline for comparison)
# All run against the SAME v77 inputs, all parameters shared across adaptors.

set -e
MODE="${1:-run}"
V77="$HOME/NK92_fit_v77/estimate_params_pzap_cleaned_up"
CA="$HOME/Ca_fit_c02"
PZ_V77="$CA/optimized_model_pzap/model_output_pzap_v77.csv"
MODES="m13_itam_iff m01_itam m12_itam_lin m02_itam_alpha c21_shared"

ENV='module load Miniconda3/4.9.2; source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh; conda activate CD16_v2'

if [ "$MODE" = "report" ]; then
  cd "$CA"
  python3 - <<'PY'
import json, os, glob
print(f'{"mode":<16} {"#par":>4} {"raw SSR":>12} {"gamma gap":>11}')
print('-' * 50)
rows = []
for f in sorted(glob.glob('out_ca_v77/results_*.json')):
    r = json.load(open(f))
    if r.get('ok'):
        rows.append((r['mode'], r.get('n_params', '?'), r['raw_ssr'], r['gamma_peak_gap_in_SE']))
    else:
        print(f'{r.get("mode","?"):<16} FAILED {r.get("error","")[:40]}')
for m, n, s, g in sorted(rows, key=lambda t: t[2]):
    print(f'{m:<16} {n:>4} {s:>12.4e} {g:>8.1f} SE')
print('-' * 50)
print('v69-input reference:  m13 = 6.01e4 (4.9 SE),  m01 = 7.61e4 (6.0 SE)')
print('shared baseline:      8.40e4 (7.4 SE)')
PY
  exit 0
fi

# ---------- step 1: export v77 pZAP ----------
if [ ! -f "$PZ_V77" ]; then
  echo "Exporting v77 pZAP curves..."
  mkdir -p "$CA/optimized_model_pzap"
  cd "$CA"
  V77DIR="$V77" OUT="$PZ_V77" python3 - <<'PY'
import os, glob, numpy as np, pandas as pd
d = os.environ['V77DIR']; out = os.environ['OUT']
runs = {'zeta': 'zeta_runs', 'gamma': 'gamma_runs', 'hetero': 'mixed_runs'}
def read(path, obs):
    try:
        names = open(path).readline().lstrip('#').split()
        a = np.loadtxt(path, comments='#')
        if a.ndim != 2 or a.shape[1] != len(names): return None, None
        c = {n: a[:, i] for i, n in enumerate(names)}
        tkey = next((n for n in names if n.lower() == 'time'), None)
        if tkey is None: return None, None
        low = {n.lower(): n for n in names}
        for k in ('pzap_total', 'pzap_bound', 'pzap', 'pzap_free'):
            if k in low: return c[tkey], c[low[k]]
        hit = next((n for n in names if 'zap' in n.lower() and 'syk' not in n.lower()), None)
        if hit is None:
            print('    columns found:', names)
        return (c[tkey], c[hit]) if hit else (None, None)
    except Exception:
        return None, None
series, tg = {}, None
for cond, sub in runs.items():
    ad = sorted(glob.glob(os.path.join(d, sub, 'analysis*')), key=os.path.getmtime)[-10:]
    tr = []
    for a in ad:
        for g in glob.glob(os.path.join(a, '**', '*.gdat'), recursive=True):
            t, v = read(g, 'pZAP_total')
            if t is not None: tr.append((t, v))
    if not tr: raise SystemExit(f'ERROR: no pZAP gdat under {sub} -- run check_fit on v77 first')
    t0 = tr[0][0]
    series[cond] = np.mean([np.interp(t0, t, v) for t, v in tr], axis=0); tg = t0
    print(f'  {cond}: {len(tr)} traces')
pd.DataFrame({'time': tg, **{f'mean_pZAP_{c}': series[c] for c in series}}).to_csv(out, index=False)
print(f'  wrote {out}')
PY
fi

# ---------- step 2: submit Ca variants in parallel ----------
cd "$CA"; mkdir -p logs
SCRIPT="$HOME/CD16_NK92_project/filesCC/filesB.py"
for M in $MODES; do

  sbatch -J ca77_$M -N1 -n1 -c 32 --mem=16G -t 2:00:00 -o logs/ca77_$M.log --wrap \
    "$ENV; cd $CA; PZAP_OVERRIDE=$PZ_V77 OUT_OVERRIDE=out_ca_v77 python $SCRIPT $M"
done
echo
echo "Submitted. squeue -u \$USER"
echo "When done:  bash \$0 report"
