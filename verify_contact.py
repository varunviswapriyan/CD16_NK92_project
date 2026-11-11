#!/usr/bin/env python3
"""
verify_contact.py -- test the explanation for the early CD8 dip:
"at day 14 many mature CD8 T cells are already in contact with tumour cells,
 so the cancer-contact exhaustion rule converts them quickly."

Uses the day-14 state (dump.2D_model.0 of the spatial runs = the initial slide).
"In contact" = same site or one of the 4 nearest neighbours has >=1 cancer cell
(the model's own definition). Reports, per sex:
  - fraction of AR- mature CD8 (i3) in contact with cancer
  - mean number of cancer cells in contact
  - expected time to exhaustion from the deck rate (0.0106 per 10-min step per
    contacting cancer cell, local + neighbour reactions)
Run: cd ~/CD16_NK92_project/SPPARKS_Sex_Bias && python3 ~/CD16_NK92_project/verify_contact.py
"""
import glob, numpy as np
L=600; RATE=0.0106; STEP_H=1/6.0          # one time unit = 10 min
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
    a=np.array(data); return {c:a[:,k] for k,c in enumerate(cols)}
for sex in ("male","female"):
    f=glob.glob(f"runs_spatial_prod/{sex}/dump.2D_model.0")
    if not f: print(sex,": no day-14 dump found"); continue
    d=read_dump(f[0])
    can=np.zeros((L,L)); mat=np.zeros((L,L))
    ix=np.clip(d["x"].astype(int),0,L-1); iy=np.clip(d["y"].astype(int),0,L-1)
    can[iy,ix]=d["i1"]; mat[iy,ix]=d["i3"]
    # cancer cells touching each site: own site + 4 neighbours (no wrap: edges blocked)
    P=np.pad(can,1)
    touch=P[1:-1,1:-1]+P[:-2,1:-1]+P[2:,1:-1]+P[1:-1,:-2]+P[1:-1,2:]
    n=mat.sum(); m=mat>0
    inc=(mat*(touch>0)).sum()
    w=touch[m]; ww=mat[m]
    mean_t=(w*ww).sum()/n
    rate=RATE*w                                    # per step, per mature cell
    with_c=rate>0
    tmed=np.median(1/rate[with_c])*STEP_H if with_c.any() else float("nan")
    print(f"\n{sex}: AR- mature CD8 at day 14 = {n:.0f}")
    print(f"  in contact with cancer       : {inc:.0f}  ({inc/n:.0%})")
    print(f"  mean cancer cells in contact : {mean_t:.1f}")
    print(f"  median expected time to exhaustion (contacting cells): {tmed:.1f} h")
print("\nIf most mature cells are in contact and the expected time is a few hours,"
      "\nthe early mature->exhausted conversion (and the CD8 dip) is explained.")
