#!/usr/bin/env bash
# check_status.sh -- one command: see what's still running, count finished
# seeds, and aggregate/plot everything that's done (no-seq, CD4-off, spatial,
# diffusion). Safe to run anytime; it just reports whatever has finished.
# Run:  cd ~/CD16_NK92_project && git pull --no-edit && bash check_status.sh
set -uo pipefail
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || { echo "repo not found"; exit 1; }
source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh 2>/dev/null && conda activate CD16_v2 2>/dev/null

# locate the plot scripts wherever they landed (repo root or thesis_repro/)
PLOTSP=$(ls ~/CD16_NK92_project/plot_spatial.py ~/CD16_NK92_project/thesis_repro/plot_spatial.py 2>/dev/null | head -1 || true)
PLOTDF=$(ls ~/CD16_NK92_project/plot_diffusion.py ~/CD16_NK92_project/thesis_repro/plot_diffusion.py 2>/dev/null | head -1 || true)

echo "=============== STILL RUNNING (job : count) ==============="
squeue -u "$USER" -h -o "%j" 2>/dev/null | sort | uniq -c || echo "  (none / queue empty)"

echo
echo "=============== FINISHED SEEDS (expect 40 each) ==============="
for d in runs_big/noseq/male_noseq runs_big/noseq/female_noseq \
         runs_big/cd4off/male_nocd4 runs_big/cd4off/female_nocd4; do
  [ -d "$d" ] && printf "  %-34s %s\n" "$d" "$(ls $d/s*/log.spparks 2>/dev/null | wc -l)"
done
for c in male female; do
  [ -d runs_spatial_prod/$c ] && printf "  %-34s %s dumps (expect 3)\n" "runs_spatial_prod/$c" "$(ls runs_spatial_prod/$c/dump.2D_model.* 2>/dev/null | wc -l)"
done

echo
echo "=============== FOLD-CHANGE RESULTS ==============="
python3 - <<'PY'
import glob, math
def fold(f):
    # returns fold only if the run REACHED day 21 (Time >= 1000); else None
    col=first=last=lasttime=None
    for line in open(f):
        s=line.split()
        if not s: continue
        if "sum(i1)" in s: col=s.index("sum(i1)"); continue
        if col is not None and len(s)>col:
            try:
                t=float(s[0]); v=float(s[col])
                if first is None: first=v
                last=v; lasttime=t
            except ValueError: pass
    if first and last and lasttime is not None and lasttime>=1000:
        return last/first
    return None
groups={
 "no-seq (IFNGR blocking)": ["runs_big/noseq/male_noseq","runs_big/noseq/female_noseq"],
 "CD4 removed":             ["runs_big/cd4off/male_nocd4","runs_big/cd4off/female_nocd4"],
}
for title,dirs in groups.items():
    print(f"-- {title} --")
    for d in dirs:
        vs=[x for x in (fold(f) for f in sorted(glob.glob(f"{d}/s*/log.spparks"))) if x]
        name=d.split('/')[-1]
        if vs:
            m=sum(vs)/len(vs); sd=math.sqrt(sum((x-m)**2 for x in vs)/(len(vs)-1)) if len(vs)>1 else 0
            print(f"   {name:14s} n={len(vs):3d}  mean fold={m:.4f}  SD={sd:.4f}")
        else:
            print(f"   {name:14s} (no finished seeds yet)")
print("\n   compare standard:  male 2.3793   female 1.5456")
PY

echo
echo "=============== SPATIAL PLOTS ==============="
if [ -n "$PLOTSP" ] && ls runs_spatial_prod/male/dump.2D_model.* >/dev/null 2>&1 && ls runs_spatial_prod/female/dump.2D_model.* >/dev/null 2>&1; then
  python3 "$PLOTSP" runs_spatial_prod && echo "  -> runs_spatial_prod/Fig_spatial_cells.png , Fig_spatial_ifng.png"
else
  echo "  not ready (dumps missing) or plot_spatial.py not found"
fi

echo
echo "=============== EFFECTIVE DIFFUSION ==============="
if [ -n "$PLOTDF" ] && [ -f diffusion_front.csv ]; then
  python3 "$PLOTDF" diffusion_front.csv
else
  echo "  diffusion_front.csv or plot_diffusion.py not found"
fi

echo
echo "To view the PNGs:  python -m http.server 8000   (then open from the browser)"
