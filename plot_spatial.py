#!/usr/bin/env python3
"""
plot_spatial.py RESDIR
Spatial figures for the prod runs. Two outputs:
  Fig_spatial_cells.png -- cancer density (red) with AR- CD8 (blue) and AR+ CD8
                           (green, enlarged so the key rare population is visible)
  Fig_spatial_ifng.png  -- free IFN-g, SHARED colour scale across all 6 panels
Columns in the dumps: id x y i1 i2 i3 i4 i5 i6 i7 d1
  i1=cancer  i2,i3,i4 = AR- CD8 (Tp,Tm,Te)  i5,i6,i7 = AR+ CD8  d1 = free IFN-g
"""
import sys, os, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

L = 600
TIMES = [(0, "Day 14"), (1, "Day 17.5"), (2, "Day 21")]
CONDS = ["male", "female"]

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
    ix=np.clip(d["x"].astype(int),0,L-1)
    iy=np.clip(d["y"].astype(int),0,L-1)
    g[iy,ix]=d[field]
    return g

def main():
    resdir=sys.argv[1] if len(sys.argv)>1 else "runs_spatial_prod"
    # first pass: cache dumps + find a shared IFN-g scale (99.9th pct, robust to spikes)
    cache={}; ifng_max=0.0
    for cond in CONDS:
        for t,_ in TIMES:
            cand=glob.glob(os.path.join(resdir,cond,f"dump.2D_model.{t}")) or \
                 glob.glob(os.path.join(resdir,cond,f"dump*{t}"))
            d=read_dump(cand[0]) if cand else None
            cache[(cond,t)]=d
            if d is not None:
                ifng_max=max(ifng_max, float(np.percentile(grid(d,"d1"),99.9)))
    if ifng_max<=0: ifng_max=1.0

    # ---- cells figure ----
    figc,axc=plt.subplots(2,3,figsize=(13,8.6))
    for r,cond in enumerate(CONDS):
        for c,(t,lab) in enumerate(TIMES):
            ax=axc[r,c]; ax.set_xticks([]); ax.set_yticks([])
            d=cache[(cond,t)]
            if d is None:
                ax.text(.5,.5,"no dump",ha="center",va="center"); continue
            cancer=grid(d,"i1")
            ix=np.clip(d["x"].astype(int),0,L-1); iy=np.clip(d["y"].astype(int),0,L-1)
            arm=d["i2"]+d["i3"]+d["i4"]            # AR- CD8
            arp=d["i5"]+d["i6"]+d["i7"]            # AR+ CD8
            ax.imshow(cancer,origin="lower",cmap="Reds",vmin=0,vmax=2,alpha=0.70)
            m=arm>0
            if m.any(): ax.scatter(ix[m],iy[m],s=2.5,c="#1565ff",linewidths=0,alpha=.75)
            m=arp>0                                 # AR+ : big, bright, drawn on top
            if m.any(): ax.scatter(ix[m],iy[m],s=14,c="#00e000",edgecolors="#0a3d0a",
                                   linewidths=.3,alpha=.95)
            ax.set_title(f"{cond} — {lab}",fontsize=11)
    figc.suptitle("Spatial cell positions   (red = cancer density,  blue = AR− CD8,  green = AR+ CD8)",
                  fontsize=13)
    figc.tight_layout()
    figc.savefig(os.path.join(resdir,"Fig_spatial_cells.png"),dpi=150)
    print("wrote Fig_spatial_cells.png")

    # ---- IFN-g figure, shared scale + one colorbar ----
    figi,axi=plt.subplots(2,3,figsize=(13,8.0))
    im=None
    for r,cond in enumerate(CONDS):
        for c,(t,lab) in enumerate(TIMES):
            ax=axi[r,c]; ax.set_xticks([]); ax.set_yticks([])
            d=cache[(cond,t)]
            if d is None:
                ax.text(.5,.5,"no dump",ha="center",va="center"); continue
            im=ax.imshow(grid(d,"d1"),origin="lower",cmap="viridis",vmin=0,vmax=ifng_max)
            ax.set_title(f"{cond} — {lab}",fontsize=11)
    figi.suptitle("Free IFN-γ concentration   (shared colour scale)",fontsize=13)
    if im is not None:
        figi.colorbar(im,ax=axi,shrink=0.72,label="free IFN-γ (per site)")
    figi.savefig(os.path.join(resdir,"Fig_spatial_ifng.png"),dpi=150)
    print("wrote Fig_spatial_ifng.png")

if __name__=="__main__":
    main()
