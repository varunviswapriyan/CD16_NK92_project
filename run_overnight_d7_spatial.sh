#!/usr/bin/env bash
# run_overnight_d7_spatial.sh -- two reviewer-motivated checks (Das, Nov 13).
#
#   bash run_overnight_d7_spatial.sh          find inputs + submit runs (minutes)
#   bash run_overnight_d7_spatial.sh analyze  (next morning) metrics + figures
#
# (A) DAY 7 -> 14. Same model and rates, started from day-7 slides instead of day 14.
#     The repo only has day-14 initial files, so step 1 SEARCHES the cluster (repo,
#     ~/giuseppe_ref, ~/rsgxg003_mount) for day-7 SPPARKS_init files. If found they are
#     run (n=20/sex); if not, the script lists any day-7 mIF data it finds so Claude can
#     write the converter. Nothing is invented.
# (B) DAY-21 SPATIAL FEATURES. Standard model, n=10/sex, with full-lattice dumps at
#     day 14 (= the mIF slide itself) and day 21. Metrics per sex, model day 21 vs day 14:
#       - densities (cells/mm^2) of AR- CD8, AR+ CD8, CD4 AR-, CD4 AR+, cancer
#       - T-cell density vs distance from the tumour edge (inside +, outside -)
#     Day-21 mIF files found in the search are listed for the data side of the comparison.
set -uo pipefail
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || exit 1
source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh 2>/dev/null && conda activate CD16_v2 2>/dev/null
SRC=rsc_out_2026-07-03; DECK="$PWD/$SRC/in.full"; BIN="$PWD/spk_sexbias"
JOB=$PWD/runs_big/_job_generic.sh
NSP="${NSP:-10}"; ND7="${ND7:-20}"; TIME="${TIME:-30:00:00}"
MS=$SRC/sex_bias_initialization_files/d14_M_4_SPPARKS_init
FS=$SRC/sex_bias_initialization_files/d14_F_4_SPPARKS_init
DUMPSED='s|^#dump .*|dump 1 sites 1008.0 dump.sp.* id x y i1 i2 i3 i4 i6 i7 i13 i14 d1|'
mkdir -p figures

if [ "${1:-}" = "analyze" ]; then
python3 - <<'PY'
import glob, os, numpy as np, csv
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
L=600; AREA_MM2=36.0; PIX_MM2=1e-4          # 10 um sites
def read(p):
    cols=None; data=[]; inv=False
    for line in open(p):
        s=line.split()
        if not s: continue
        if s[0]=="Values": inv=True; continue
        if cols is None and s[0]=="id": cols=[c for c in s if c!="values"]; continue
        if inv and cols and len(s)>=len(cols):
            try: data.append([float(x) for x in s[:len(cols)]])
            except ValueError: pass
    a=np.array(data); return {c:a[:,i] for i,c in enumerate(cols)}
def grid(d,keys):
    g=np.zeros((L,L)); ix=np.clip(d["x"].astype(int),0,L-1); iy=np.clip(d["y"].astype(int),0,L-1)
    for k in keys: np.add.at(g,(iy,ix),d[k])
    return g
def box_blur(a,k):
    c=np.cumsum(np.cumsum(np.pad(a,((k+1,k),(k+1,k))),0),1)
    return (c[2*k+1:,2*k+1:]-c[:-2*k-1,2*k+1:]-c[2*k+1:,:-2*k-1]+c[:-2*k-1,:-2*k-1])/(2*k+1)**2
def signed_dist(mask):
    # BFS-style distance (sites) to the mask edge; + inside, - outside
    def dist(m):
        d=np.where(m,np.inf,0.0); k=0
        while True:
            k+=1; P=np.pad(d,1,constant_values=0.0)
            nb=np.minimum.reduce([P[:-2,1:-1],P[2:,1:-1],P[1:-1,:-2],P[1:-1,2:]])+1
            new=np.where(m,np.minimum(d,nb),0.0)
            if np.array_equal(new,d) or k>L: return d
            d=new
    return dist(mask)-dist(~mask)
TYPES={"AR- CD8":["i2","i3","i4"],"AR+ CD8":["i6","i7"],"CD4 AR-":["i13"],"CD4 AR+":["i14"],"cancer":["i1"]}
BINS=np.arange(-500,801,100)                 # um from tumour edge
rows=[]; prof={}
for sex in ("male","female"):
    for tag,idx in (("day14",0),("day21",1008)):
        files=sorted(glob.glob(f"runs_spatial21/{sex}/s*/dump.sp.{idx}"))
        if not files: continue
        dens={k:[] for k in TYPES}; P=[]
        for f in files:
            d=read(f)
            can=grid(d,["i1"]); mask=box_blur((can>0).astype(float),15)>0.5
            sd=signed_dist(mask)*10.0
            T=grid(d,["i2","i3","i4","i6","i7"])
            for k,v in TYPES.items(): dens[k].append(grid(d,v).sum()/AREA_MM2)
            p=[]
            for a,b in zip(BINS[:-1],BINS[1:]):
                sel=(sd>=a)&(sd<b); p.append(T[sel].sum()/(sel.sum()*PIX_MM2) if sel.any() else np.nan)
            P.append(p)
        prof[(sex,tag)]=np.nanmean(np.array(P),0)
        for k in TYPES:
            v=np.array(dens[k]); rows.append([sex,tag,k,len(v),f"{v.mean():.2f}",f"{v.std(ddof=1)/np.sqrt(len(v)) if len(v)>1 else 0:.2f}"])
if rows:
    with open("spatial21_metrics.csv","w",newline="") as fh:
        w=csv.writer(fh); w.writerow(["sex","time","type","n","density_per_mm2","SE"]); w.writerows(rows)
    print("=== densities (cells/mm^2, whole 36 mm^2 box) ===")
    for r in rows: print("  %-6s %-5s %-8s n=%-3s %8s +/- %s"%tuple(r))
    fig,ax=plt.subplots(figsize=(6.4,4.8)); mid=(BINS[:-1]+BINS[1:])/2
    for (sex,tag),p in prof.items():
        ax.plot(mid,p,("-" if tag=="day21" else "--"),color=("#1f4db8" if sex=="male" else "#b81f5b"),lw=2,label=f"{sex}, {tag}")
    ax.axvline(0,color="grey",lw=1); ax.set_xlabel("Distance from tumour edge (µm; + = inside)")
    ax.set_ylabel("CD8+ T cells per mm²"); ax.set_title("CD8+ T-cell infiltration profile (model)")
    ax.legend(frameon=False); ax.grid(True,color="#eeeeee",lw=.6); fig.tight_layout()
    fig.savefig("figures/Infiltration_profile_model.png",dpi=200)
    print("wrote spatial21_metrics.csv, figures/Infiltration_profile_model.png")
else:
    print("(B) no spatial dumps yet")
# (A) day 7 -> 14
def fold(f):
    hdr=None; first=last=t=None
    for line in open(f):
        s=line.split()
        if "sum(i1)" in s:
            if hdr: break
            hdr=s; k=s.index("sum(i1)"); continue
        if hdr and len(s)==len(hdr):
            try: tt=float(s[0]); v=float(s[k])
            except ValueError: continue
            first=v if first is None else first; last=v; t=tt
    return last/first if first and t is not None and t>=1000 else None
res={}
for sex in ("male","female"):
    v=[x for x in (fold(f) for f in glob.glob(f"runs_d7/{sex}/s*/log.spparks")) if x]
    if v: res[sex]=np.mean(v); print(f"(A) day7->14 {sex:6s} n={len(v)} fold={np.mean(v):.3f}")
if len(res)==2: print(f"(A) day7->14 ratio male/female = {res['male']/res['female']:.2f}   (day14->21 standard: 1.54)")
elif not glob.glob("runs_d7/*"): print("(A) no day-7 runs (no day-7 initial files were found)")
PY
exit 0
fi

[ -x "$BIN" ] && [ -f "$JOB" ] || { echo "missing $BIN or $JOB"; exit 1; }
REPORT=overnight_d7_spatial_report.txt; : > $REPORT
say(){ echo "$@" | tee -a $REPORT; }

# ---------------- search for day-7 / day-21 inputs ----------------
say "=== searching for day-7 and day-21 slide files ==="
ROOTS=". $HOME/giuseppe_ref $HOME/rsgxg003_mount $HOME/CD16_NK92_project"
D7=$(find $ROOTS -maxdepth 7 -type f \( -iname "*d7*SPPARKS_init*" -o -iname "*day7*SPPARKS_init*" -o -iname "*d07*SPPARKS_init*" \) 2>/dev/null | grep -v runs_ | sort -u)
say "-- day-7 SPPARKS_init files:"; [ -n "$D7" ] && echo "$D7" | sed 's/^/   /' | tee -a $REPORT || say "   (none)"
say "-- other day-7 / day-21 slide data (mIF exports, csv, init builders):"
find $ROOTS -maxdepth 7 -type f \( -iname "*d7_*" -o -iname "*day7*" -o -iname "*day_7*" -o -iname "*d21_*" -o -iname "*day21*" -o -iname "*day_21*" -o -iname "*SPPARKS_init*.py" -o -iname "*make*init*" -o -iname "*slide*.py" \) 2>/dev/null \
  | grep -v -E "runs_|dump\.|log\.spparks" | sort -u | head -60 | sed 's/^/   /' | tee -a $REPORT

submit(){ # name base n slide fem extra-sed
  local name=$1 base=$2 n=$3 sl=$4 fem=$5 ex=$6
  [ -f $base/cfg.sh ] && { say "  $base already submitted"; return; }
  mkdir -p $base; printf '%s\n' "$ex" > $base/extra.sed
  cat > $base/cfg.sh <<EOF
BASE=$PWD/$base
SEEDBASE=9000
SLIDEPATH=$(readlink -f $sl)
SLIDENAME=$(basename $sl)
FEM=$fem
DECK=$DECK
BIN=$BIN
EXTRASED=$PWD/$base/extra.sed
EOF
  sbatch --job-name=$name --array=1-$n --time=$TIME --ntasks=1 --output=$PWD/$base/slurm_%A_%a.out \
         --export=ALL,CFG=$PWD/$base/cfg.sh "$JOB" >/dev/null && say "  submitted $name x$n -> $base"
}

# ---------------- (B) day-21 spatial runs ----------------
say "=== (B) standard runs with day-14 and day-21 lattice dumps ==="
sed "$DUMPSED" "$DECK" | grep -q "^dump 1 sites 1008.0" || { say "!! no '#dump' line in deck"; exit 1; }
submit sp21_m runs_spatial21/male   $NSP $MS 0 "$DUMPSED"
submit sp21_f runs_spatial21/female $NSP $FS 1 "$DUMPSED"

# ---------------- (A) day 7 -> 14 ----------------
say "=== (A) day 7 -> 14 ==="
M7=$(echo "$D7" | grep -i -E "_M_|male" | grep -v -i madef | head -1)
F7=$(echo "$D7" | grep -i -E "_F_|female" | head -1)
if [ -n "$M7" ] && [ -n "$F7" ]; then
  say "  male slide:   $M7"; say "  female slide: $F7"
  submit d7_m runs_d7/male   $ND7 "$M7" 0 ""
  submit d7_f runs_d7/female $ND7 "$F7" 1 ""
else
  say "  NOT submitted: no day-7 SPPARKS_init files for both sexes."
  say "  Send Claude this report; the day-7 slides must first be converted from the mIF data."
fi
say ""; say "Report: SPPARKS_Sex_Bias/$REPORT   Tomorrow: bash run_overnight_d7_spatial.sh analyze"
