#!/usr/bin/env python3
"""
plot_folds.py -- fold-change-vs-day figures for the sensitivity runs, each
drawn against the standard prod male/female (dashed reference). One PNG each:
  IFNGR blocking (no-seq), equal AR+/AR- recruitment, CD4+ removed.
Run:  cd ~/CD16_NK92_project/SPPARKS_Sex_Bias && python3 ~/CD16_NK92_project/plot_folds.py
"""
import glob, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

DAYGRID = np.linspace(14.0, 21.0, 141)
BLUE, PINK = "#1f4db8", "#b81f5b"
BLUE_L, PINK_L = "#9bb4e8", "#e8a7c2"

def read_log(path):
    tcol=ccol=None; t=[]; c=[]
    for line in open(path):
        s=line.split()
        if not s: continue
        if "sum(i1)" in s:
            if t: break                       # first stats block only
            ccol=s.index("sum(i1)"); tcol=s.index("Time") if "Time" in s else 0
            continue
        if ccol is not None and len(s)>ccol:
            try: tv=float(s[tcol]); cv=float(s[ccol])
            except ValueError: continue
            if t and tv<=t[-1]: break
            t.append(tv); c.append(cv)
    if len(c)<3: return None
    return np.array(t), np.array(c)

def curve(dirpath):
    folds=[]
    for f in sorted(glob.glob(os.path.join(dirpath,"s*","log.spparks"))):
        r=read_log(f)
        if r is None: continue
        t,c=r; days=14.0+t/144.0
        if days[-1]<20.9: continue            # finished runs only
        folds.append(np.interp(DAYGRID, days, c/c[0]))
    if not folds: return None
    F=np.array(folds)
    return F.mean(0), F.std(0,ddof=1)/np.sqrt(len(F)), len(F)

std={lab:curve(d) for lab,d in (("male","runs_big/prod/male"),
                                ("female","runs_big/prod/female"))}

EXPTS=[
 ("IFNGR blocking (no sequestration)",
   [("male","runs_big/noseq/male_noseq",BLUE),("female","runs_big/noseq/female_noseq",PINK)],
   "runs_big/noseq/Fig_folds_ifngr_block.png"),
 ("Equal AR+/AR- recruitment",
   [("male","runs_big/equalrec/male_eq",BLUE),("female","runs_big/equalrec/female",PINK)],
   "runs_big/equalrec/Fig_folds_equalrec.png"),
 ("CD4+ T cells removed",
   [("male","runs_big/cd4off/male_nocd4",BLUE),("female","runs_big/cd4off/female_nocd4",PINK)],
   "runs_big/cd4off/Fig_folds_cd4off.png"),
]

for title, conds, out in EXPTS:
    fig,ax=plt.subplots(figsize=(6,5))
    for lab,col in (("male",BLUE_L),("female",PINK_L)):
        if std.get(lab):
            m,se,n=std[lab]; ax.plot(DAYGRID,m,"--",color=col,lw=1.5,label=f"{lab} — standard")
    got=False
    for lab,sub,col in conds:
        r=curve(sub)
        if r is None: continue
        got=True; m,se,n=r
        ax.plot(DAYGRID,m,"-",color=col,lw=2.3,label=f"{lab} — this condition (n={n})")
        ax.fill_between(DAYGRID,m-se,m+se,color=col,alpha=.25,lw=0)
    ax.set(xlabel="Day",ylabel="Cancer cell population fold change",title=title)
    ax.legend(frameon=False,fontsize=8); ax.grid(True,color="#e6e6e6",lw=.6)
    fig.tight_layout()
    if got: fig.savefig(out,dpi=200); print("wrote",out)
    else:   print("no finished data for",title)
