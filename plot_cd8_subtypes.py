#!/usr/bin/env python3
"""
plot_cd8_subtypes.py -- why does total CD8 dip until ~day 16 and then rise?
Splits CD8 into progenitor / mature / exhausted (AR- and AR+) over time,
mean of all finished seeds, standard model, male and female.
Columns (from the deck): i2 AR- progenitor, i3 AR- mature, i4 AR- exhausted,
                         i5 AR+ progenitor (always 0), i6 AR+ mature, i7 AR+ exhausted.
Writes figures/CD8_subtypes_standard.png and prints a table (day 14/16/18/21).
Run: cd ~/CD16_NK92_project/SPPARKS_Sex_Bias && python3 ~/CD16_NK92_project/plot_cd8_subtypes.py
"""
import glob, os, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
DAYS=np.linspace(14,21,141)
COLS=["sum(i2)","sum(i3)","sum(i4)","sum(i6)","sum(i7)"]
NAMES={"sum(i2)":"AR- progenitor","sum(i3)":"AR- mature","sum(i4)":"AR- exhausted",
       "sum(i6)":"AR+ mature","sum(i7)":"AR+ exhausted"}
COL={"sum(i2)":"#2ca02c","sum(i3)":"#1f77b4","sum(i4)":"#d62728",
     "sum(i6)":"#17becf","sum(i7)":"#ff7f0e"}

def read(f):
    hdr=None; rows=[]
    for line in open(f):
        s=line.split()
        if not s: continue
        if "sum(i1)" in s:
            if rows: break
            hdr=s; continue
        if hdr and len(s)==len(hdr):
            try: v=[float(x) for x in s]
            except ValueError: continue
            if rows and v[0]<=rows[-1][0]: break
            rows.append(v)
    if not hdr or len(rows)<3 or rows[-1][0]<1000: return None
    a=np.array(rows); d=14+a[:,0]/144.0
    return {c:np.interp(DAYS,d,a[:,hdr.index(c)]) for c in COLS if c in hdr}

os.makedirs("figures",exist_ok=True)
fig,axs=plt.subplots(1,2,figsize=(13,5),sharey=False)
for ax,sex in zip(axs,("male","female")):
    R=[r for r in (read(f) for f in sorted(glob.glob(f"runs_big/prod/{sex}/s*/log.spparks"))) if r]
    if not R: print("no finished runs for",sex); continue
    tot=np.zeros_like(DAYS)
    print(f"\n{sex} (n={len(R)})   day:      14       16       18       21")
    for c in COLS:
        m=np.mean([r[c] for r in R],axis=0); tot+=m
        ax.plot(DAYS,m,color=COL[c],lw=2,label=NAMES[c])
        print(f"  {NAMES[c]:16s}        "+"  ".join(f"{m[i]:7.0f}" for i in (0,40,80,140)))
    ax.plot(DAYS,tot,"k--",lw=2,label="total CD8")
    print(f"  {'total CD8':16s}        "+"  ".join(f"{tot[i]:7.0f}" for i in (0,40,80,140)))
    ax.set_title(f"CD8 subtypes, standard model — {sex} (n={len(R)})")
    ax.set_xlabel("Day"); ax.set_ylabel("cells"); ax.grid(True,color="#eee",lw=.6)
    ax.legend(frameon=False,fontsize=9)
fig.tight_layout(); fig.savefig("figures/CD8_subtypes_standard.png",dpi=170)
print("\nwrote figures/CD8_subtypes_standard.png")
