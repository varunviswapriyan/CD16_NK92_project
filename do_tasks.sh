#!/usr/bin/env bash
# do_tasks.sh -- Nov-4 task list (Das / JD) in one script.
#
#   bash do_tasks.sh          -> does everything possible from EXISTING data now,
#                                and submits the new runs that are unavoidable
#   bash do_tasks.sh plot     -> (after the queue empties) makes the remaining plots
#
# Tasks covered
#   0  report: log columns + boundary-condition evidence (deck + SPPARKS source)
#   1  population kinetics (cancer, AR-/AR+ CD8, CD4 helpers, free IFN-g)
#      for all conditions  -> figures/Kinetics_*.png + table in tasks_report.txt
#   2  growth curves: Giuseppe's original (thesis) params vs ours, LINEAR x and y
#                         -> figures/Growth_thesis_vs_ours.png
#   3  IFN-g map with tumour outline + stretched colour scale; cell map with
#      smaller AR+ / bigger AR- dots -> figures/Figure3_overlay.png, Figure2_cells.png
#   4  true IFNGR blockade on ALL T cells (receptors = 0 on AR- and AR+ T cells,
#      cancer cells keep theirs), male+female n=40 -> runs_big/blockall/
#
# New simulations are needed ONLY for:
#   - task 4 (never run before), and
#   - task 1, IF the existing logs only record cancer (sum(i1)). The script checks
#     the log header first; if T-cell / IFN-g columns are already there it uses
#     the existing 40-seed runs and submits nothing for task 1.
set -uo pipefail
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || { echo "repo not found"; exit 1; }
source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh 2>/dev/null && conda activate CD16_v2 2>/dev/null
SRC=rsc_out_2026-07-03
DECK="$PWD/$SRC/in.full"
BIN="$PWD/spk_sexbias"
NKIN="${NKIN:-10}"        # seeds per condition for kinetics-only reruns
NBLOCK="${NBLOCK:-40}"    # seeds per sex for the full IFNGR blockade
TIME="${TIME:-24:00:00}"
PART="${PART:-}"
REPORT=tasks_report.txt
mkdir -p figures
DIAGCOLS="i2 sum i3 sum i4 sum i5 sum i6 sum i7 sum i8 sum i9 sum i10 sum i11 sum d1 sum"

# ---------------------------------------------------------------- plotting ---
make_plots() {
python3 - <<'PY'
import glob, os, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import PowerNorm
BLUE,PINK="#1f4db8","#b81f5b"
DAYS=np.linspace(14,21,141)
out=open("tasks_report.txt","a")
def say(*a):
    s=" ".join(str(x) for x in a); print(s); out.write(s+"\n")

def read_log(f):
    hdr=None; rows=[]
    for line in open(f):
        s=line.split()
        if not s: continue
        if "sum(i1)" in s:
            if rows: break
            hdr=s; continue
        if hdr and len(s)==len(hdr):
            try: v=[float(x) for x in s[:len(hdr)]]
            except ValueError: continue
            if rows and v[0]<=rows[-1][0]: break
            rows.append(v)
    if not hdr or len(rows)<3: return None
    a=np.array(rows); return {h:a[:,k] for k,h in enumerate(hdr)}

def group(dirglob, need_day21=True):
    runs=[]
    for f in sorted(glob.glob(dirglob)):
        r=read_log(f)
        if r is None: continue
        t=r.get("Time", list(r.values())[0])
        d=14+t/144.0
        if need_day21 and t[-1]<1000: continue
        runs.append({k:np.interp(DAYS,d,v) for k,v in r.items()})
    return runs

# ---------------- task 2: growth, thesis params vs ours, linear axes --------
fig,ax=plt.subplots(figsize=(6.4,5))
got=False
for lab,pat,ls in (("Giuseppe original params","runs_sweep/thesis/{}_s*/log.spparks","--"),
                   ("manuscript params","runs_big/prod/{}/s*/log.spparks","-")):
    for sex,col in (("male",BLUE),("female",PINK)):
        R=group(pat.format(sex))
        if not R: say(f"[task2] no finished runs for {lab} {sex} ({pat.format(sex)})"); continue
        F=np.array([r["sum(i1)"]/r["sum(i1)"][0] for r in R]); got=True
        ax.plot(DAYS,F.mean(0),ls,color=col,lw=2.2,label=f"{sex}, {lab} (n={len(F)})")
        say(f"[task2] {lab:26s} {sex:6s} n={len(F):3d} fold d21={F.mean(0)[-1]:.3f}")
ax.set_xlabel("Day"); ax.set_ylabel("Cancer cell fold change (linear)")
ax.set_title("Tumour growth: original vs manuscript parameters")
ax.legend(frameon=False,fontsize=8); ax.grid(True,color="#e6e6e6",lw=.6); fig.tight_layout()
if got: fig.savefig("figures/Growth_thesis_vs_ours.png",dpi=200); say("[task2] wrote figures/Growth_thesis_vs_ours.png")

# ---------------- task 1: kinetics ----------------------------------------
CONDS=[("standard","runs_big/prod/male","runs_big/prod/female"),
       ("AR+ internalization off","runs_big/noseq/male_noseq","runs_big/noseq/female_noseq"),
       ("CD4 helpers removed (corrected)","runs_big/cd4off_v2/male","runs_big/cd4off_v2/female"),
       ("IFNGR blocked, all T cells","runs_big/blockall/male","runs_big/blockall/female")]
PANELS=[("Cancer cells",["sum(i1)"]),
        ("AR- CD8",["sum(i2)","sum(i3)","sum(i4)"]),
        ("AR+ CD8",["sum(i5)","sum(i6)","sum(i7)"]),
        ("AR+ fraction of CD8",None),
        ("Total CD8",["sum(i2)","sum(i3)","sum(i4)","sum(i5)","sum(i6)","sum(i7)"]),
        ("CD4 helpers AR-",["sum(i13)"]),
        ("CD4 helpers AR+",["sum(i14)"]),
        ("Free IFN-γ",["sum(d1)"])]
def series(r,cols):
    if cols is None:
        a=sum(r[c] for c in ("sum(i5)","sum(i6)","sum(i7)"))
        b=sum(r[c] for c in ("sum(i2)","sum(i3)","sum(i4)"))
        return a/np.maximum(a+b,1)
    return sum(r[c] for c in cols)
free={}
for name,md,fd in CONDS:
    data={sex:group(f"{d}/s*/log.spparks") for sex,d in (("male",md),("female",fd))}
    data={k:[r for r in v if "sum(i2)" in r] for k,v in data.items()}
    if not any(data.values()): continue
    fig,axs=plt.subplots(2,4,figsize=(16,7.5))
    for ax,(title,cols) in zip(axs.flat,PANELS):
        for sex,col in (("male",BLUE),("female",PINK)):
            R=data[sex]
            if not R: continue
            try: S=np.array([series(r,cols) for r in R])
            except KeyError: continue
            m=S.mean(0); se=S.std(0,ddof=1)/math.sqrt(len(S)) if len(S)>1 else 0*m
            ax.plot(DAYS,m,color=col,lw=2,label=f"{sex} (n={len(S)})")
            ax.fill_between(DAYS,m-se,m+se,color=col,alpha=.2,lw=0)
            if title.startswith("Free"): free[(name,sex)]=m
        ax.set_title(title,fontsize=11); ax.set_xlabel("Day"); ax.grid(True,color="#eee",lw=.6)
    axs.flat[0].legend(frameon=False,fontsize=9)
    fig.suptitle(f"Population kinetics — {name}",fontsize=14); fig.tight_layout(rect=[0,0,1,.95])
    fn="figures/Kinetics_"+name.replace(" ","_").replace("+","p").replace("(","").replace(")","").replace(",","")+".png"
    fig.savefig(fn,dpi=150); say(f"[task1] wrote {fn}")
    say(f"[task1] {name}: values at day 14 / 17.5 / 21 (mean over seeds)")
    for sex in ("male","female"):
        R=data[sex]
        if not R: continue
        for title,cols in PANELS:
            try: m=np.array([series(r,cols) for r in R]).mean(0)
            except KeyError: continue
            i=[0,70,140]
            fmt=(lambda x:f"{x:7.1%}") if cols is None else (lambda x:f"{x:12.4g}")
            say(f"      {sex:6s} {title:26s} "+"  ".join(fmt(m[k]) for k in i))
# free IFN-g comparison across conditions (answers Fig 3 / Fig 8 questions)
if free:
    fig,axs=plt.subplots(1,2,figsize=(12,4.6),sharey=False)
    for ax,sex in zip(axs,("male","female")):
        for (name,s),m in free.items():
            if s==sex: ax.plot(DAYS,m,lw=2,label=name)
        ax.set_yscale("log"); ax.set_ylim(bottom=1e5)
        ax.set_title(f"Free IFN-γ — {sex}"); ax.set_xlabel("Day"); ax.set_ylabel("total free IFN-γ (log)")
        ax.legend(frameon=False,fontsize=8)
        ax.grid(True,color="#eee",lw=.6)
    fig.tight_layout(); fig.savefig("figures/Kinetics_free_IFNg_by_condition.png",dpi=150)
    say("[task1] wrote figures/Kinetics_free_IFNg_by_condition.png")

# ---------------- task 3: spatial overlays + dot sizes -----------------------
L=600; TIMES=[(0,"Day 14"),(1,"Day 17.5"),(2,"Day 21")]
def read_dump(p):
    cols=None; data=[]; inv=False
    for line in open(p):
        s=line.split()
        if not s: continue
        if s[0]=="Values": inv=True; continue
        if cols is None and s[0]=="id": cols=[c for c in s if c!="values"]; continue
        if inv and cols and len(s)>=len(cols):
            try: data.append([float(x) for x in s[:len(cols)]])
            except ValueError: pass
    if not cols or not data: return None
    a=np.array(data); return {c:a[:,k] for k,c in enumerate(cols)}
def grid(d,f):
    g=np.zeros((L,L)); ix=np.clip(d["x"].astype(int),0,L-1); iy=np.clip(d["y"].astype(int),0,L-1)
    g[iy,ix]=d[f]; return g
try:
    from scipy.ndimage import uniform_filter
except ImportError:
    def uniform_filter(a,k):
        from numpy.lib.stride_tricks import sliding_window_view
        p=np.pad(a,k//2,mode='edge'); return sliding_window_view(p,(k,k)).mean((-1,-2))
cache={}
for sex in ("male","female"):
    for t,_ in TIMES:
        c=glob.glob(f"runs_spatial_prod/{sex}/dump.2D_model.{t}")
        cache[(sex,t)]=read_dump(c[0]) if c else None
if any(v is not None for v in cache.values()):
    vmax=max(np.percentile(grid(d,"d1"),99.9) for d in cache.values() if d is not None) or 1
    fig,axs=plt.subplots(2,3,figsize=(13,8.2)); im=None
    for r,sex in enumerate(("male","female")):
        for c,(t,lab) in enumerate(TIMES):
            ax=axs[r,c]; ax.set_xticks([]); ax.set_yticks([]); d=cache[(sex,t)]
            if d is None: continue
            im=ax.imshow(grid(d,"d1"),origin="lower",cmap="magma",norm=PowerNorm(0.35,vmin=0,vmax=vmax))
            can=grid(d,"i1")
            ax.contour(uniform_filter((can>0).astype(float),31),levels=[0.5],colors="white",linewidths=1.0)
            ax.set_title(f"{sex} — {lab}",fontsize=11)
    fig.suptitle("Free IFN-γ and tumour edge (white line)",fontsize=15)
    if im is not None: fig.colorbar(im,ax=axs,shrink=.72,label="free IFN-γ per site (nonlinear scale)")
    fig.savefig("figures/Figure3_overlay.png",dpi=150); say("[task3] wrote figures/Figure3_overlay.png")

    from matplotlib.patches import Patch
    from matplotlib.lines import Line2D
    fig,axs=plt.subplots(2,3,figsize=(13,8.6))
    for r,sex in enumerate(("male","female")):
        for c,(t,lab) in enumerate(TIMES):
            ax=axs[r,c]; ax.set_xticks([]); ax.set_yticks([]); d=cache[(sex,t)]
            if d is None: continue
            ix=np.clip(d["x"].astype(int),0,L-1); iy=np.clip(d["y"].astype(int),0,L-1)
            arm=d["i2"]+d["i3"]+d["i4"]; arp=d["i5"]+d["i6"]+d["i7"]
            ax.imshow(grid(d,"i1"),origin="lower",cmap="Reds",vmin=0,vmax=2,alpha=.7)
            m=arm>0; ax.scatter(ix[m],iy[m],s=4.0,c="#1565ff",linewidths=0,alpha=.7)
            m=arp>0; ax.scatter(ix[m],iy[m],s=4.0,c="#00b000",linewidths=0,alpha=.95)
            ax.set_title(f"{sex} — {lab}  (AR+ {int(arp.sum())}, AR- {int(arm.sum())})",fontsize=10)
    h=[Patch(facecolor="#c0392b",label="cancer"),
       Line2D([0],[0],marker='o',color='none',markerfacecolor="#1565ff",markersize=7,label="AR- CD8"),
       Line2D([0],[0],marker='o',color='none',markerfacecolor="#00b000",markersize=7,label="AR+ CD8")]
    fig.suptitle("Spatial cell positions",fontsize=16,y=.99)
    fig.legend(handles=h,loc="upper center",ncol=3,frameon=False,fontsize=12,bbox_to_anchor=(.5,.955))
    fig.tight_layout(rect=[0,0,1,.92]); fig.savefig("figures/Figure2_cells.png",dpi=150)
    say("[task3] wrote figures/Figure2_cells.png")
else:
    say("[task3] no dumps in runs_spatial_prod/")

# ---------------- task 4 summary --------------------------------------------
def fold(f):
    r=read_log(f)
    if r is None: return None
    t=r.get("Time",list(r.values())[0])
    return r["sum(i1)"][-1]/r["sum(i1)"][0] if t[-1]>=1000 else None
for sex in ("male","female"):
    v=[x for x in (fold(f) for f in glob.glob(f"runs_big/blockall/{sex}/s*/log.spparks")) if x]
    if v:
        m=np.mean(v); se=np.std(v,ddof=1)/math.sqrt(len(v)) if len(v)>1 else 0
        say(f"[task4] IFNGR blocked on all T cells  {sex:6s} n={len(v):3d} fold={m:.4f} ± {se:.4f} (SE)")
say("[task4] compare: standard male 2.3793 female 1.5456 | AR+ internalization off male 1.7689 female 1.5022")
out.close()
PY
}

if [ "${1:-}" = "plot" ]; then
  rm -f figures/Kinetics_CD4_helpers_removed.png figures/Kinetics_standard_kin.png figures/Kinetics_OLD_cd4off_run_=_recruitment_markers_removed.png
  echo "=== queue ==="; squeue -u "$USER" -h -o "%j" 2>/dev/null | sort | uniq -c
  echo "=== plots (log appended to $REPORT) ==="
  { echo; echo "######## plot pass $(date)"; } >> $REPORT
  make_plots
  echo; echo "View:  cd ~/CD16_NK92_project/SPPARKS_Sex_Bias && python -m http.server 8001"
  exit 0
fi

if [ "${1:-}" = "fix" ]; then
  echo "=== cancelling the unneeded kinetics reruns and the first blockade submission ==="
  scancel -u "$USER" -n kin_std_m,kin_std_f,kin_nsq_m,kin_nsq_f,kin_cd4_m,kin_cd4_f,blk_m,blk_f 2>/dev/null
  sleep 5; rm -rf runs_big/kin runs_big/blockall; echo "   removed runs_big/kin and runs_big/blockall"
fi

# --------------------------------------------------------------- task 0 -----
{
echo "######## do_tasks.sh  $(date)"
echo "=== log header (what the existing runs record) ==="
LOG1=$(ls runs_big/prod/male/s*/log.spparks 2>/dev/null | head -1)
HDR=$(grep -m1 "sum(i1)" "$LOG1" 2>/dev/null || true)
echo "  $LOG1"; echo "  $HDR"
echo "=== deck: boundary / lattice / diag lines ==="
grep -niE "boundary|lattice|region|create_box|diag_style|^dump|stats" "$DECK" | sed 's/^/  /'
echo "=== SPPARKS app source: periodic / boundary-crossing code ==="
for f in $(find . -path ./runs_big -prune -o -name "app_rxn_diff*.cpp" -print -o -name "app_lattice.cpp" -print 2>/dev/null | head -5); do
  echo "  -- $f"
  grep -nE "&& *false|periodic|xprd|boundary|nx *< *0|>= *nx" "$f" | head -25 | sed 's/^/     /'
done
} | tee $REPORT

LOG1=$(ls runs_big/prod/male/s*/log.spparks 2>/dev/null | head -1)
HDR=$(grep -m1 "sum(i1)" "$LOG1" 2>/dev/null || true)
HAVE_KIN=0; echo "$HDR" | grep -q "sum(i2)" && HAVE_KIN=1
echo "  kinetics columns already in logs: $HAVE_KIN" | tee -a $REPORT
grep -q "^diag_style *array" "$DECK" && CANDIAG=1 || CANDIAG=0

# ------------------------------------------------------- job submission -----
[ -x "$BIN" ] || { echo "binary $BIN missing -- build first"; exit 1; }
JOB=$PWD/runs_big/_job_generic.sh
mkdir -p runs_big
cat > "$JOB" <<'EOF'
#!/usr/bin/env bash
source "$CFG"
module load foss/2023a
seed=$((SEEDBASE + SLURM_ARRAY_TASK_ID))
d=$BASE/s$seed
mkdir -p "$d/sex_bias_initialization_files"
cp -f "$SLIDEPATH" "$d/sex_bias_initialization_files/$SLIDENAME"
sed -e "/^seed/ s/[0-9][0-9]*/$seed/" \
    -e "s|read_sites.*|read_sites sex_bias_initialization_files/$SLIDENAME|" \
    -f "$EXTRASED" "$DECK" > "$d/in.simulation"
if [ "$FEM" = 1 ]; then
  awk '{ if ($0 ~ /recruitment/ && $0 !~ /app_style/){ if($0 ~ /AR\+/)$6=0; else if($0 ~ /Tp recruitment/||$0 ~ /Tm recruitment/)$6=$6/0.63 }
         if ($0 ~ /Tp differentiation/ && $0 !~ /app_style/){ if($0 ~ /AR\+ Tm/)$6=0; else $6=$6*2 } print }' \
      "$d/in.simulation" > "$d/in.tmp" && mv "$d/in.tmp" "$d/in.simulation"
fi
cd "$d"; mpirun -np 1 "$BIN" < in.simulation > log.spparks 2>&1
EOF
chmod +x "$JOB"

# submit COND BASE N SLIDE FEM SEDLINES...
submit() {
  local name=$1 base=$2 n=$3 slidepath=$4 fem=$5; shift 5
  if [ -f "$base/cfg.sh" ]; then echo "   $name already submitted ($base) -- skipping"; return; fi
  mkdir -p "$base"
  printf '%s\n' "$@" > "$base/extra.sed"
  cat > "$base/cfg.sh" <<EOF
BASE=$PWD/$base
SEEDBASE=9000
SLIDEPATH=$PWD/$slidepath
SLIDENAME=$(basename "$slidepath")
FEM=$fem
DECK=$DECK
BIN=$BIN
EXTRASED=$PWD/$base/extra.sed
EOF
  # dry-run check that every edit actually hits the deck
  sed -f "$base/extra.sed" "$DECK" > "$base/deck_check.txt"
  sbatch --job-name="$name" --array=1-$n --time=$TIME --ntasks=1 \
         --output="$PWD/$base/slurm_%A_%a.out" ${PART:+--partition=$PART} \
         --export=ALL,CFG="$PWD/$base/cfg.sh" "$JOB" >/dev/null \
    && echo "   submitted $name x$n -> $base" | tee -a $REPORT
}
MSLIDE=$SRC/sex_bias_initialization_files/d14_M_4_SPPARKS_init
FSLIDE=$SRC/sex_bias_initialization_files/d14_F_4_SPPARKS_init

# ---- task 4: IFNGR blocked on all T cells (receptor args = AR- T, AR+ T, cancer)
BLOCK='/app_style/ s/ 2\.0 1000 1000 1000 / 2.0 0 0 1000 /'
if sed "$BLOCK" "$DECK" | grep -q " 2\.0 0 0 1000 "; then
  echo "=== task 4: submitting IFNGR-block (all T cells) ===" | tee -a $REPORT
  submit blk_m runs_big/blockall/male   $NBLOCK $MSLIDE 0 "$BLOCK"
  submit blk_f runs_big/blockall/female $NBLOCK $FSLIDE 1 "$BLOCK"
else
  echo "!!! task 4 NOT submitted: could not find ' 2.0 1000 1000 1000 ' on the app_style line." | tee -a $REPORT
  grep -n "app_style" "$DECK" | tee -a $REPORT
fi

# ---- corrected CD4-helper removal: helpers are slide columns i13 (Th) and i14 (ARTh)
#      = awk fields 16 and 17. (The old run zeroed fields 10/11 = i8 recruitment
#      marker and i9 error signal, i.e. it removed T-cell recruitment, not CD4.)
mkdir -p runs_big/cd4off_v2/slides
NOCD4M=runs_big/cd4off_v2/slides/d14_M_4_SPPARKS_init_noCD4
NOCD4F=runs_big/cd4off_v2/slides/d14_F_4_SPPARKS_init_noCD4
for pair in "$MSLIDE:$NOCD4M" "$FSLIDE:$NOCD4F"; do
  a=${pair%%:*}; b=${pair##*:}
  awk 'f && NF>15 {$16=0; $17=0} /^Values/{f=1} {print}' "$a" > "$b"
  echo "   $(basename $b): CD4 entries before=$(awk 'f&&NF>15{c+=($16!=0)+($17!=0)} /^Values/{f=1} END{print c+0}' "$a") after=$(awk 'f&&NF>15{c+=($16!=0)+($17!=0)} /^Values/{f=1} END{print c+0}' "$b")  recruit markers kept=$(awk 'f&&NF>15{c+=($10!=0)} /^Values/{f=1} END{print c+0}' "$b")" | tee -a $REPORT
done
submit cd4v2_m runs_big/cd4off_v2/male   $NBLOCK $NOCD4M 0 ""
submit cd4v2_f runs_big/cd4off_v2/female $NBLOCK $NOCD4F 1 ""

echo | tee -a $REPORT
echo "=== plots possible right now (tasks 2, 3, and 1 if columns exist) ===" | tee -a $REPORT
make_plots

echo
echo "Done. Report: SPPARKS_Sex_Bias/$REPORT   Figures: SPPARKS_Sex_Bias/figures/"
echo "When 'squeue -u \$USER' is empty:   cd ~/CD16_NK92_project && bash do_tasks.sh plot"
