#!/usr/bin/env bash
# aggregate.sh -- read every replicate in runs_big/prod/<cond>/s*/log.spparks,
# report mean fold +/- SD, standard error, and the male/female ratio with a
# propagated standard error. Run after run_prod_array.sh jobs finish.
# Run:  cd ~/CD16_NK92_project/SPPARKS_Sex_Bias && bash aggregate.sh
set -uo pipefail
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || { echo "repo not found"; exit 1; }
source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh 2>/dev/null && conda activate CD16_v2 2>/dev/null

python3 - <<'PY'
import glob, math
from collections import defaultdict
folds=defaultdict(list)
for f in sorted(glob.glob("runs_big/prod/*/s*/log.spparks")):
    col=first=last=None
    for line in open(f):
        s=line.split()
        if not s: continue
        if "sum(i1)" in s:
            col=s.index("sum(i1)"); continue
        if col is not None and len(s)>col:
            try:
                float(s[0]); v=float(s[col])
                if first is None: first=v
                last=v
            except ValueError:
                pass
    if first and last:
        cond=f.split("/")[2]
        folds[cond].append(last/first)

def stats(v):
    n=len(v); m=sum(v)/n
    sd=math.sqrt(sum((x-m)**2 for x in v)/(n-1)) if n>1 else 0.0
    se=sd/math.sqrt(n) if n>1 else 0.0
    return n,m,sd,se

print(f"{'condition':12s} {'n':>3s} {'mean':>8s} {'SD':>8s} {'SE':>8s}")
S={}
for c in ("male","female","male_noseq"):
    if c in folds:
        n,m,sd,se=stats(folds[c]); S[c]=(m,se)
        print(f"{c:12s} {n:3d} {m:8.4f} {sd:8.4f} {se:8.4f}")

if "male" in S and "female" in S:
    (mm,mse),(fm,fse)=S["male"],S["female"]
    r=mm/fm
    rse=r*math.sqrt((mse/mm)**2+(fse/fm)**2) if mm and fm else 0.0
    print(f"\nmale/female ratio = {r:.4f} +/- {rse:.4f}  (experiment ~1.48)")
PY
