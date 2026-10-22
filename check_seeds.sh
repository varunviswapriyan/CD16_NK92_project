#!/usr/bin/env bash
# check_seeds.sh -- pull the start / day-21 cancer count and fold for every prod
# seed, so we can see whether the two seeds actually diverge (real replicates) or
# are identical (RNG seed not reaching the trajectory).
# Run:  cd ~/CD16_NK92_project/SPPARKS_Sex_Bias && bash check_seeds.sh
set -uo pipefail
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || { echo "repo not found"; exit 1; }
source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh 2>/dev/null && conda activate CD16_v2 2>/dev/null

python3 - <<'PY'
import glob
from collections import defaultdict
folds=defaultdict(list)
for f in sorted(glob.glob("runs_sweep/prod/*_s*/log.spparks")):
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
    name=f.split("/")[2]
    cond="_".join(name.split("_")[:-1])    # male / female / male_noseq
    if first and last:
        fold=last/first
        folds[cond].append(fold)
        print(f"{name:22s} start={first:9.0f} end={last:9.0f} fold={fold:.4f}")

print("\n--- per-condition mean fold (n seeds) ---")
m=folds.get("male",[]); fm=folds.get("female",[])
for c,v in folds.items():
    import statistics as st
    sd=st.pstdev(v) if len(v)>1 else 0.0
    print(f"{c:12s} mean={sum(v)/len(v):.4f}  sd={sd:.4f}  n={len(v)}")
if m and fm:
    print(f"\nmale/female ratio of means = {(sum(m)/len(m))/(sum(fm)/len(fm)):.4f}")
    if len(m)>1 or len(fm)>1:
        print("seeds DIVERGE -> real replicates, safe to scale up." if (len(set(round(x,6) for x in m))>1 or len(set(round(x,6) for x in fm))>1)
              else "seeds IDENTICAL to 6 dp -> RNG seed likely not threaded; fix before averaging.")
PY
