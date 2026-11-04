#!/usr/bin/env bash
# plot_growth_gstyle.sh -- tumour growth (manuscript params, n=40) drawn the way
# Giuseppe drew thesis Fig 4.x: x = time in hours (0-168), LOG y-axis with ticks
# 1.0, 1.2, ... , all individual seeds faint + mean bold, dashed line at 1,
# bracket with the male/female ratio at the end.
# Run: cd ~/CD16_NK92_project && git pull --no-edit && bash plot_growth_gstyle.sh
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
        hrs=np.array(t)/6.0                      # 1 time unit = 10 min
        out.append(np.interp(H,hrs,np.array(c)/c[0]))
    return np.array(out)
M=runs("runs_big/prod/male"); F=runs("runs_big/prod/female")
fig,ax=plt.subplots(figsize=(8.5,7))
for A,col,lab in ((M,"#0000ff","male mouse"),(F,"#ff00ff","female mouse")):
    for r in A: ax.plot(H,r,color=col,alpha=.08,lw=1)
    ax.plot(H,A.mean(0),color=col,lw=2.2,label=f"{lab} (n={len(A)})")
ax.axhline(1,ls="--",color="grey",lw=1.5)
ax.set_yscale("log")
top=max(M.mean(0)[-1],F.mean(0)[-1])
ticks=np.round(np.arange(1.0,top+0.2,0.2),1)
ax.yaxis.set_major_locator(FixedLocator(ticks)); ax.yaxis.set_minor_locator(NullLocator())
ax.yaxis.set_major_formatter(FuncFormatter(lambda v,_: r"$10^0$" if abs(v-1)<1e-9 else rf"${v:.1f}\times10^0$"))
ax.set_ylim(0.95,top*1.04); ax.set_xlim(-8,178)
ax.set_xticks([0,25,50,75,100,125,150,175])
ax.set_xlabel("Time (hours)",fontsize=14,fontweight="bold")
ax.set_ylabel("Cancer Cell\nPopulation Fold Change",fontsize=14,fontweight="bold")
ax.tick_params(labelsize=12)
ratio=M.mean(0)[-1]/F.mean(0)[-1]; ym,yf=M.mean(0)[-1],F.mean(0)[-1]
ax.plot([170,173,173,170],[ym,ym,yf,yf],color="#3d6a85",lw=3,clip_on=False)
ax.text(176,np.sqrt(ym*yf),f"x{ratio:.2f}",fontsize=22,fontweight="bold",va="center",clip_on=False)
ax.legend(loc="upper left",fontsize=13)
fig.tight_layout(); fig.subplots_adjust(right=0.84)
fig.savefig("figures/Figure4_giuseppe_style.png",dpi=200)
print(f"wrote figures/Figure4_giuseppe_style.png   male {ym:.3f}  female {yf:.3f}  ratio x{ratio:.2f}")
PY
