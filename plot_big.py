#!/usr/bin/env python3
"""
plot_big.py -- growth figure from the n=40 array runs in runs_big/prod/.
Reads every runs_big/prod/<cond>/s*/log.spparks, builds the cancer fold-change
curve per seed (sum(i1)/sum(i1)[0]), averages across seeds, and draws mean +/-
SE bands for male / female / male_noseq. Two panels: fold change and absolute.
Writes runs_big/prod/Fig_growth_big.png.
Run:  cd ~/CD16_NK92_project/SPPARKS_Sex_Bias && python3 ~/CD16_NK92_project/plot_big.py
      (or wherever plot_big.py landed; it uses an absolute results path)
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

def read_log(path):
    """return (time_array, cancer_array) from a SPPARKS stats log."""
    tcol = ccol = None
    t, c = [], []
    for line in open(path):
        s = line.split()
        if not s:
            continue
        if "sum(i1)" in s:
            ccol = s.index("sum(i1)")
            tcol = s.index("Time") if "Time" in s else 0
            continue
        if ccol is not None and len(s) > ccol:
            try:
                tv = float(s[tcol]); cv = float(s[ccol])
                t.append(tv); c.append(cv)
            except ValueError:
                pass
    if not c:
        return None
    return np.array(t), np.array(c)

def collect(cond):
    folds, absolutes, tgrid = [], [], None
    for f in sorted(glob.glob(os.path.join(RES, cond, "s*", "log.spparks"))):
        r = read_log(f)
        if r is None:
            continue
        t, c = r
        if tgrid is None:
            tgrid = t
        n = min(len(tgrid), len(c))
        folds.append(c[:n] / c[0])
        absolutes.append(c[:n])
    if not folds:
        return None
    n = min(len(x) for x in folds)
    folds = np.array([x[:n] for x in folds])
    absolutes = np.array([x[:n] for x in absolutes])
    days = 14 + tgrid[:n] / 144.0          # 144 native units per day
    return days, folds, absolutes

fig, ax = plt.subplots(1, 2, figsize=(12, 5))
for cond, col, lab in COND:
    d = collect(cond)
    if d is None:
        continue
    days, folds, absol = d
    for axi, Y, scale in ((ax[0], folds, 1.0), (ax[1], absol, 1e3)):
        m = Y.mean(0) / scale
        se = Y.std(0, ddof=1) / np.sqrt(Y.shape[0]) / scale
        axi.plot(days, m, color=col, lw=2, label=f"{lab} (n={Y.shape[0]})")
        axi.fill_between(days, m - se, m + se, color=col, alpha=0.2, lw=0)

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
