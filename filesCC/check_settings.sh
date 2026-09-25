#!/bin/bash
# check_settings.sh -- show what every pZAP bootstrap run on disk is ACTUALLY
# configured with, and how far along it is.  Read-only; safe at any time.
#
#   bash ~/CD16_NK92_project/filesCC/check_settings.sh
set -e
PY=$(command -v python3 || command -v python)
SUB=estimate_params_pzap_cleaned_up

for d in ~/boot_pzap ~/boot_pzap_n16 ~/boot_pzap_t600; do
  [ -d "$d" ] || continue
  "$PY" - "$d" "$SUB" <<'PYEOF'
import os, sys, json, glob
root, sub = sys.argv[1], sys.argv[2]
name = os.path.basename(root)
cfg = None
for s in sorted(glob.glob(os.path.join(root, '*[0-9]'))):
    p = os.path.join(s, sub, 'v_config.json')
    if os.path.exists(p):
        cfg = json.load(open(p)); break
if cfg is None:
    print('%-16s (no samples built)' % name); raise SystemExit
FIT6 = {'lig0','kdl0','ZAP0','SYK0','KZP_MULT','KPR_MULT'}
def canon(k): return 'kdl0' if k.strip().lower() in ('kd10','kdl0') else k.strip()
done = 0
for s in sorted(glob.glob(os.path.join(root, '*[0-9]'))):
    f = os.path.join(s, sub, 'analysis_param_residue.dat')
    if not os.path.exists(f): continue
    for ln in open(f, errors='ignore'):
        if ln.strip().lower().startswith('linear'):
            got = {canon(t.split('=',1)[0]) for t in ln.split()[1:] if '=' in t}
            if FIT6 <= got: done += 1
            break
tot = len(glob.glob(os.path.join(root, '*[0-9]')))
box = {n: (round(10**lo, 4), round(10**hi, 4))
       for n, lo, hi in zip(cfg['PARAMS'], cfg['LB'], cfg['UB'])}
print('%-16s  N_REPS=%-3s FIT_TMAX=%-6s  %d of %d samples finished'
      % (name, cfg.get('N_REPS'), cfg.get('FIT_TMAX'), done, tot))
print('%-16s  fitted: %s' % ('', cfg['PARAMS']))
print('%-16s  fixed : %s' % ('', cfg.get('FIXED')))
print('%-16s  box   : %s' % ('', box))
print()
PYEOF
done

echo "jobs on the nodes, by output directory:"
for id in $(squeue -u "$USER" -h -o "%i" 2>/dev/null); do
  scontrol show job "$id" 2>/dev/null | tr ' ' '\n' | grep '^StdOut=' | sed 's|^StdOut=||'
done | sed 's|.*/\(boot_pzap[^/]*\)/.*|   \1|' | sort | uniq -c
