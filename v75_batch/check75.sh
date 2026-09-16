#!/bin/bash
# check75.sh -- inspect the v75..v96 batch. One command, three jobs:
#
#   bash ~/CD16_NK92_project/v75_batch/check75.sh          # ranked table + fitted values
#   bash ~/CD16_NK92_project/v75_batch/check75.sh plots    # submit check_fit for top 4
#   bash ~/CD16_NK92_project/v75_batch/check75.sh serve    # start web server for the PNGs
#
# The table shows cost, residue, worst single-point miss (if check_fit has run),
# and the fitted parameter values, so you can see at a glance whether anything
# is pinned at a bound. v69 (WITH the affinity fraction) is the benchmark row.

set -e
MODE="${1:-table}"
VARIANTS="75 76 77 78 79 80 81 82 83 84 85 86 87 88 89 90 91 92 93 94 95 96"
SUB="estimate_params_pzap_cleaned_up"
TOP="77 89 92 86"     # best FRAC-free candidates so far

if [ "$MODE" = "serve" ]; then
  IP=$(hostname -I | awk '{print $2}')
  echo "Collecting plots into ~/v75_plots/ ..."
  mkdir -p ~/v75_plots
  for V in $VARIANTS 69; do
    for P in ~/NK92_fit_v$V/$SUB/plot_N.png ~/NK92_fit_v$V/$SUB/plot_run*.png; do
      [ -f "$P" ] && cp "$P" ~/v75_plots/v${V}_$(basename $P) 2>/dev/null || true
    done
  done
  ls -1 ~/v75_plots/ | head -30
  echo
  echo "Open:  http://$IP:8000/"
  echo "(Ctrl+C here to stop the server when done)"
  cd ~/v75_plots && python -m http.server 8000
  exit 0
fi

if [ "$MODE" = "plots" ]; then
  echo "Submitting check_fit (10 replicates) for top candidates: $TOP"
  for V in $TOP; do
    D=$HOME/NK92_fit_v$V/$SUB
    [ -d "$D" ] || { echo "  v$V missing, skipping"; continue; }
    sbatch -J cf$V -N1 -n1 -c 32 --mem=32G -t 1:00:00 \
      -o $HOME/NK92_fit_v$V/checkfit_v$V.log --wrap \
      "module load Miniconda3/4.9.2; source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh; conda activate CD16_v2; cd $D; python check_fit.py 10"
  done
  echo
  echo "~10 min. Then:  bash \$0 table     (worst-miss column fills in)"
  echo "                bash \$0 serve     (view the plots)"
  exit 0
fi

# ---------------- table ----------------
python3 - "$HOME" "$SUB" $VARIANTS <<'PY'
import sys, os, json, re, glob
home, sub = sys.argv[1], sys.argv[2]
variants = sys.argv[3:]

def grab(v):
    d = os.path.join(home, f'NK92_fit_v{v}')
    logs = glob.glob(os.path.join(d, '*.log')) + glob.glob(os.path.join(d, sub, 'slurm*.out'))
    txt = ''
    for f in sorted(logs, key=os.path.getmtime):
        try: txt += open(f, errors='ignore').read()
        except Exception: pass
    out = {}
    m = re.findall(r'best cost:\s*([0-9.eE+-]+)', txt)
    if m: out['cost'] = float(m[-1])
    m = re.findall(r'residue \*\* \(1/4\) =\s*([0-9.eE+-]+)', txt)
    if m: out['res'] = float(m[-1])
    m = re.findall(r'best pos:\s*\[([^\]]+)\]', txt.replace('\n', ' '))
    if m: out['pos'] = [float(x) for x in m[-1].split()]
    # worst single-point miss printed by check_fit
    m = re.findall(r'worst[^0-9]*([0-9.]+)', txt, re.I)
    if m: out['worst'] = float(m[-1])
    cfg = os.path.join(d, sub, 'v_config.json')
    if os.path.exists(cfg):
        c = json.load(open(cfg))
        out['note'] = c.get('NOTE', '')
        out['params'] = c.get('PARAMS', [])
        out['fixed'] = c.get('FIXED', {})
        out['lb'] = c.get('LB', []); out['ub'] = c.get('UB', [])
    return out

rows = []
for v in variants:
    r = grab(v)
    if r.get('cost') is not None:
        rows.append((v, r))
rows.sort(key=lambda t: t[1]['cost'])

print(f'{"V":<5} {"cost":>10} {"res":>8} {"worst":>7} {"#fit":>4}  description')
print('-' * 84)
for v, r in rows:
    w = f'{r["worst"]:.3f}' if 'worst' in r else '  -  '
    print(f'v{v:<4} {r["cost"]:>10.5f} {r.get("res",0):>8.4f} {w:>7} '
          f'{len(r.get("params",[])):>4}  {r.get("note","")[:40]}')
print('-' * 84)
print(f'{"v69":<5} {0.00244:>10.5f} {0.2222:>8.4f} {0.043:>7} {7:>4}  BENCHMARK (with affinity fraction)')

print('\nFitted values for top 4 (flag = at a bound):')
for v, r in rows[:4]:
    if 'pos' not in r: continue
    print(f'\n  v{v}: {r.get("note","")}')
    for n, x, lo, hi in zip(r['params'], r['pos'], r['lb'], r['ub']):
        flag = ''
        if abs(x - lo) < 0.02: flag = '  <-- AT LOWER BOUND'
        if abs(x - hi) < 0.02: flag = '  <-- AT UPPER BOUND'
        print(f'      {n:<12} log10={x:+.4f}  value={10**x:.4g}{flag}')
    if r.get('fixed'):
        print(f'      fixed: {r["fixed"]}')
PY

echo
echo "Next:  bash \$0 plots    (run check_fit on the top candidates)"
echo "       bash \$0 serve    (view all plots in a browser)"
