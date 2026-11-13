#!/usr/bin/env bash
# plot_fboff_fold.sh -- tumour growth fold change vs time, male/female, with IFNGR
# blocked on all T cells (solid) over the standard model (dashed). Same style and
# axes as Figure4_giuseppe_style.png (time in hours, log y, ratio bracket).
# Also makes the same plot for the other conditions (for the SI), in case needed.
# Run: cd ~/CD16_NK92_project && git pull --no-edit && bash plot_fboff_fold.sh
cd ~/CD16_NK92_project/SPPARKS_Sex_Bias || exit 1
source /gpfs0/scratch/miniforge3/24.11.2/etc/profile.d/conda.sh 2>/dev/null && conda activate CD16_v2 2>/dev/null
mkdir -p figures
python3 - <<'PY'
import glob, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullLocator, FuncFormatter
H=np.linspace(0,168,337)
def runs(d):
    out=[]
    for f in sorted(glob.glob(f"{d}/s*/log.spparks")):
        hdr=None; t=[]; c=[]
        for line in open(f):
            s=line.split()
            if not s: continue
            if "sum(i1)" in s:
                if t: break
                hdr=s; k=s.index("sum(i1)"); continue
            if hdr and len(s)==len(hdr):
                try: tv=float(s[0]); cv=float(s[k])
                except ValueError: continue
                if t and tv<=t[-1]: break
                t.append(tv); c.append(cv)
        if len(t)<3 or t[-1]<1000: continue
        out.append(np.interp(H,np.array(t)/6.0,np.array(c)/c[0]))
    return np.array(out)
STD={"male":runs("runs_big/prod/male"),"female":runs("runs_big/prod/female")}
CONDS=[("feedback off","runs_big/fboff","Fold_feedback_off.png")]
SUB={"runs_big/noseq":("male_noseq","female_noseq")}
for lab,base,fn in CONDS:
    ms,fs=SUB.get(base,("male","female"))
    M=runs(f"{base}/{ms}"); F=runs(f"{base}/{fs}")
    if len(M)==0 or len(F)==0: print("no finished runs for",lab); continue
    fig,ax=plt.subplots(figsize=(8.5,7))
    for A,col,sex in ((STD["male"],"#0000ff","male"),(STD["female"],"#ff00ff","female")):
        ax.plot(H,A.mean(0),"--",color=col,lw=1.6,alpha=.6,label=f"{sex} mouse, standard (n={len(A)})")
    for A,col,sex in ((M,"#0000ff","male"),(F,"#ff00ff","female")):
        for r in A: ax.plot(H,r,color=col,alpha=.06,lw=1)
        ax.plot(H,A.mean(0),color=col,lw=2.4,label=f"{sex} mouse, {lab} (n={len(A)})")
    ax.axhline(1,ls="--",color="grey",lw=1.2)
    ax.set_yscale("log")
    top=max(M.mean(0)[-1],F.mean(0)[-1],STD["male"].mean(0)[-1])
    ticks=np.round(np.arange(1.0,top+0.2,0.2),1)
    ax.yaxis.set_major_locator(FixedLocator(ticks)); ax.yaxis.set_minor_locator(NullLocator())
    ax.yaxis.set_major_formatter(FuncFormatter(lambda v,_: r"$10^0$" if abs(v-1)<1e-9 else rf"${v:.1f}\times10^0$"))
    ax.set_ylim(0.95,top*1.05); ax.set_xlim(-8,178); ax.set_xticks([0,25,50,75,100,125,150,175])
    ax.set_xlabel("Time (hours)",fontsize=14,fontweight="bold")
    ax.set_ylabel("Cancer Cell\nPopulation Fold Change",fontsize=14,fontweight="bold")
    ax.tick_params(labelsize=12)
    ym,yf=M.mean(0)[-1],F.mean(0)[-1]; ratio=ym/yf
    ax.plot([170,173,173,170],[ym,ym,yf,yf],color="#3d6a85",lw=3,clip_on=False)
    ax.text(176,np.sqrt(ym*yf)*(1.03 if abs(ym-yf)<0.15 else 1),f"x{ratio:.2f}",fontsize=20,fontweight="bold",va="center",clip_on=False)
    ax.legend(loc="upper left",fontsize=10.5)
    fig.tight_layout(); fig.subplots_adjust(right=0.84)
    fig.savefig(f"figures/{fn}",dpi=200)
    print(f"wrote figures/{fn}   {lab}: male {ym:.3f}  female {yf:.3f}  ratio x{ratio:.2f}")
PY
echo "View: cd ~/CD16_NK92_project/SPPARKS_Sex_Bias && python -m http.server 8001  ->  /figures/Fold_feedback_off.png"
