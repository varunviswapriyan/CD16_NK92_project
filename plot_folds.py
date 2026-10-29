#!/usr/bin/env python3
"""
Fold-change figures -> figures/Figure3.png (growth m/f),
Figure4.png (IFNGR blocking), Figure6.png (equal recruitment),
Figure7.png (CD4 removed). Solid = this experiment's condition, named
explicitly; dashed = the baseline. Change REF_LABEL to rename the baseline.
Run:  cd ~/CD16_NK92_project/SPPARKS_Sex_Bias && python3 ~/CD16_NK92_project/plot_folds.py
"""
import glob, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT="figures"; os.makedirs(OUT, exist_ok=True)
DAYGRID=np.linspace(14.0,21.0,141)
BLUE,PINK="#1f4db8","#b81f5b"
BLUE_L,PINK_L="#9bb4e8","#e8a7c2"
REF_LABEL="standard"          # <-- dashed baseline label; e.g. "baseline" / "full model"

def read_log(path):
    tcol=ccol=None; t=[]; c=[]
    for line in open(path):
        s=line.split()
        if not s: continue
        if "sum(i1)" in s:
            if t: break
            ccol=s.index("sum(i1)"); tcol=s.index("Time") if "Time" in s else 0
            continue
        if ccol is not None and len(s)>ccol:
            try: tv=float(s[tcol]); cv=float(s[ccol])
            except ValueError: continue
            if t and tv<=t[-1]: break
            t.append(tv); c.append(cv)
    if len(c)<3: return None
    return np.array(t),np.array(c)

def curve(d):
    folds=[]
    for f in sorted(glob.glob(os.path.join(d,"s*","log.spparks"))):
        r=read_log(f)
        if r is None: continue
        t,c=r; days=14.0+t/144.0
        if days[-1]<20.9: continue
        folds.append(np.interp(DAYGRID,days,c/c[0]))
    if not folds: return None
    F=np.array(folds)
    return F.mean(0),F.std(0,ddof=1)/np.sqrt(len(F)),len(F)

std={"male":curve("runs_big/prod/male"),"female":curve("runs_big/prod/female")}

def fold_fig(title, conds, out, reference=True):
    fig,ax=plt.subplots(figsize=(6,5))
    if reference:
        for lab,col in (("male",BLUE_L),("female",PINK_L)):
            if std.get(lab):
                m,se,n=std[lab]; ax.plot(DAYGRID,m,"--",color=col,lw=1.6,label=f"{lab} ({REF_LABEL})")
    got=False
    for d,col,clab in conds:
        r=curve(d)
        if r is None: continue
        got=True; m,se,n=r
        ax.plot(DAYGRID,m,"-",color=col,lw=2.4,label=clab)
        ax.fill_between(DAYGRID,m-se,m+se,color=col,alpha=.22,lw=0)
    ax.set_xlabel("Day"); ax.set_ylabel("Cancer cell population fold change")
    ax.set_title(title); ax.legend(frameon=False,fontsize=9); ax.grid(True,color="#e6e6e6",lw=.6)
    fig.tight_layout()
    if got:
        p=os.path.join(OUT,out); fig.savefig(p,dpi=200); print("wrote",p)
    else:
        print("no finished data for",title)

# Figure 3 -- main tumour growth, male vs female (no baseline overlay)
fold_fig("Tumour growth (fold change)",
         [("runs_big/prod/male",BLUE,"male (n=40)"),
          ("runs_big/prod/female",PINK,"female (n=40)")],
         "Figure3.png", reference=False)

# Figure 4 -- IFNGR blocking
fold_fig("IFNGR blocking",
         [("runs_big/noseq/male_noseq",BLUE,"male (IFNGR blocked)"),
          ("runs_big/noseq/female_noseq",PINK,"female (IFNGR blocked)")],
         "Figure4.png")

# Figure 6 -- equal AR+/AR- recruitment (SI)
fold_fig("Equal AR+/AR- recruitment",
         [("runs_big/equalrec/male_eq",BLUE,"male (equal recruitment)"),
          ("runs_big/equalrec/female",PINK,"female (equal recruitment)")],
         "Figure6.png")

# Figure 7 -- CD4+ removed (SI)
fold_fig("CD4+ T cells removed",
         [("runs_big/cd4off/male_nocd4",BLUE,"male (CD4 removed)"),
          ("runs_big/cd4off/female_nocd4",PINK,"female (CD4 removed)")],
         "Figure7.png")
