#!/usr/bin/env python3
"""
plot_big.py -- growth figure from the n=40 array runs in runs_big/prod/.
HARDENED: parses only the FIRST monotonic stats block of each log (ignores any
repeated header / time reset), and interpolates every seed onto a common day
grid before averaging, so differing log lengths can't create phantom curves.
Writes runs_big/prod/Fig_growth_big.png.
Run:  cd ~/CD16_NK92_project && python3 plot_big.py
"""
import glob, os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RES = os.path.expanduser("~/CD16_NK92_project/SPPARKS_Sex_Bias/runs_big/prod")
COND = [("male", "#1f4db8", "male"),
        ("female", "#b81f5b", "female"),
        ("male_noseq", "#5aa469", "male (no sequestration)")]
DAYGRID = np.linspace(14.0, 21.0, 141)     # common axis for all seeds

def read_log(path):
    """First monotonic (Time, sum(i1)) block only."""
    tcol = ccol = None
    t, c = [], []
    for line in open(path):
        s = line.split()
        if not s:
            continue
        if "sum(i1)" in s:
            if ccol is not None:          # a SECOND header -> stop, keep block 1
                break
            ccol = s.index("sum(i1)")
            tcol = s.index("Time") if "Time" in s else 0
            continue
        if ccol is not None and len(s) > ccol:
            try:
                tv = float(s[tcol]); cv = float(s[ccol])
            except ValueError:
                continue
            if t and tv <= t[-1]:          # time reset/non-monotonic -> stop
                break
            t.append(tv); c.append(cv)
    if len(c) < 3:
        return None
    return np.array(t), np.array(c)

def collect(cond):
    folds, absol = [], []
    for f in sorted(glob.glob(os.path.join(RES, cond, "s*", "log.spparks"))):
        r = read_log(f)
        if r is None:
            continue
        t, c = r
        days = 14.0 + t / 144.0            # 144 native units per day
        if days[-1] < 20.9:                # unfinished run -> skip
            continue
        fold = c / c[0]
        folds.append(np.interp(DAYGRID, days, fold))
        absol.append(np.interp(DAYGRID, days, c))
    if not folds:
        return None
    return np.array(folds), np.array(absol)

fig, ax = plt.subplots(1, 2, figsize=(12, 5))
for cond, col, lab in COND:
    d = collect(cond)
    if d is None:
        print("no data for", cond); continue
    folds, absoluteC = d
    nseed = folds.shape[0]
    for axi, Y, scale in ((ax[0], folds, 1.0), (ax[1], absoluteC, 1e3)):
        m = Y.mean(0) / scale
        se = Y.std(0, ddof=1) / np.sqrt(nseed) / scale
        axi.plot(DAYGRID, m, color=col, lw=2, label=f"{lab} (n={nseed})")
        axi.fill_between(DAYGRID, m - se, m + se, color=col, alpha=0.25, lw=0)

ax[0].set(xlabel="Day", ylabel="Cancer cell population fold change",
          title="Tumour growth (fold change)")
ax[1].set(xlabel="Day", ylabel="Tumour size (thousands of cancer cells)",
          title="Absolute tumour size")
for a in ax:
    a.legend(frameon=False); a.grid(True, color="#e6e6e6", lw=.6)
fig.tight_layout()
out = os.path.join(RES, "Fig_growth_big.png")
fig.savefig(out, dpi=200)
print("wrote", out)
