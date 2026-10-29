#!/usr/bin/env python3
"""plot_spatial.py RESDIR  -> figures/Figure1.png (cells), figures/Figure2.png (IFN-g)."""
import sys, os, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.lines import Line2D

L = 600
TIMES = [(0, "Day 14"), (1, "Day 17.5"), (2, "Day 21")]
CONDS = ["male", "female"]
OUT = "figures"; os.makedirs(OUT, exist_ok=True)

def read_dump(path):
    cols=None; data=[]; inv=False
    with open(path) as f:
        for line in f:
            s=line.split()
            if not s: continue
            if s[0]=="Values": inv=True; continue
            if cols is None and s[0]=="id":
                cols=[c for c in s if c!="values"]; continue
            if inv and cols is not None and len(s)>=len(cols):
                try: data.append([float(x) for x in s[:len(cols)]])
                except ValueError: pass
    if cols is None or not data: return None
    arr=np.array(data)
    return {c: arr[:,k] for k,c in enumerate(cols)}

def grid(d, field):
    g=np.zeros((L,L))
    ix=np.clip(d["x"].astype(int),0,L-1); iy=np.clip(d["y"].astype(int),0,L-1)
    g[iy,ix]=d[field]; return g

def main():
    resdir=sys.argv[1] if len(sys.argv)>1 else "runs_spatial_prod"
    cache={}; ifng_max=0.0
    for cond in CONDS:
        for t,_ in TIMES:
            cand=glob.glob(os.path.join(resdir,cond,f"dump.2D_model.{t}")) or \
                 glob.glob(os.path.join(resdir,cond,f"dump*{t}"))
            d=read_dump(cand[0]) if cand else None
            cache[(cond,t)]=d
            if d is not None: ifng_max=max(ifng_max,float(np.percentile(grid(d,"d1"),99.9)))
    if ifng_max<=0: ifng_max=1.0

    # ---- Figure 1: cells ----
    figc,axc=plt.subplots(2,3,figsize=(13,8.6))
    for r,cond in enumerate(CONDS):
        for c,(t,lab) in enumerate(TIMES):
            ax=axc[r,c]; ax.set_xticks([]); ax.set_yticks([])
            d=cache[(cond,t)]
            if d is None: ax.text(.5,.5,"no dump",ha="center",va="center"); continue
            cancer=grid(d,"i1")
            ix=np.clip(d["x"].astype(int),0,L-1); iy=np.clip(d["y"].astype(int),0,L-1)
            arm=d["i2"]+d["i3"]+d["i4"]; arp=d["i5"]+d["i6"]+d["i7"]
            ax.imshow(cancer,origin="lower",cmap="Reds",vmin=0,vmax=2,alpha=0.70)
            m=arm>0
            if m.any(): ax.scatter(ix[m],iy[m],s=2.5,c="#1565ff",linewidths=0,alpha=.75)
            m=arp>0
            if m.any(): ax.scatter(ix[m],iy[m],s=14,c="#00c000",edgecolors="#0a3d0a",linewidths=.3,alpha=.95)
            ax.set_title(f"{cond} — {lab}",fontsize=11)
    handles=[Patch(facecolor="#c0392b",label="cancer"),
             Line2D([0],[0],marker='o',color='none',markerfacecolor="#1565ff",markersize=8,label="AR- CD8"),
             Line2D([0],[0],marker='o',color='none',markerfacecolor="#00c000",markeredgecolor="#0a3d0a",markersize=9,label="AR+ CD8")]
    figc.suptitle("Spatial cell positions", fontsize=16, y=0.99)
    figc.legend(handles=handles, loc="upper center", ncol=3, frameon=False, fontsize=12, bbox_to_anchor=(0.5,0.955))
    figc.tight_layout(rect=[0,0,1,0.92])
    figc.savefig(os.path.join(OUT,"Figure1.png"),dpi=150); print("wrote figures/Figure1.png")

    # ---- Figure 2: IFN-g (shared scale, no scale-note in title) ----
    figi,axi=plt.subplots(2,3,figsize=(13,8.0)); im=None
    for r,cond in enumerate(CONDS):
        for c,(t,lab) in enumerate(TIMES):
            ax=axi[r,c]; ax.set_xticks([]); ax.set_yticks([])
            d=cache[(cond,t)]
            if d is None: ax.text(.5,.5,"no dump",ha="center",va="center"); continue
            im=ax.imshow(grid(d,"d1"),origin="lower",cmap="viridis",vmin=0,vmax=ifng_max)
            ax.set_title(f"{cond} — {lab}",fontsize=11)
    figi.suptitle("Free IFN-γ concentration", fontsize=15)
    if im is not None: figi.colorbar(im,ax=axi,shrink=0.72,label="free IFN-γ (per site)")
    figi.savefig(os.path.join(OUT,"Figure2.png"),dpi=150); print("wrote figures/Figure2.png")

if __name__=="__main__": main()
