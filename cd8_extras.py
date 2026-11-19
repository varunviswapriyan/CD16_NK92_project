#!/usr/bin/env python3
"""
cd8_extras.py -- the two follow-ups from the Nov 19 meeting, data vs model, as bar graphs
(no lines between time points). Run from ~/CD16_NK92_project/SPPARKS_Sex_Bias:
    python3 ~/CD16_NK92_project/cd8_extras.py zihai_data/Processed_Zihai_Data

1. CD8 change between time points: ratio day 14 / day 7 and day 21 / day 14, male and female.
   Data: mean over mice of CD8 per slide, and CD8 per mm2 of tissue (30 um squares, Giuseppe's area).
   Model: total CD8 at the end / start of the day 7->14 runs and of the day 14->21 runs (mean over runs).
   -> figures/CD8_ratio_bars.png
2. CD8 density in the tumour periphery vs core (Das's suggestion: a fixed distance from the boundary).
   Edge = boundary of the cancer mask (same definition as the infiltration figure);
   periphery = within PERI um of the edge on either side; core = more than PERI um inside.
   Data: day 14 and 21 slides. Model: day 14 (start) and day 21 snapshots of runs_spatial21.
   -> figures/CD8_periphery_core_bars.png
Numbers -> cd8_extras_summary.txt
"""
import sys, os, glob, numpy as np
from convert_zihai import read_slide, types, bins
from zihai_analysis import read_log, read_dump, box_blur, edt, UM
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

PERI = 200.0                      # um, periphery half-width around the tumour edge
COL = {"male": "#1f4db8", "female": "#b81f5b"}
out = []


def regions(cancer, cd8):
    mask = box_blur((cancer > 0).astype(float), 15) > 0.5
    sd = (edt(mask) - edt(~mask)) * UM
    peri = np.abs(sd) <= PERI; core = sd > PERI
    a = lambda m: m.sum() * (UM / 1000) ** 2
    return (cd8[peri].sum() / a(peri) if a(peri) > 0 else np.nan, cd8[core].sum() / a(core) if a(core) > 0.3 else np.nan)


# ------------------------------------------------------------------ data
ddir = sys.argv[1] if len(sys.argv) > 1 else "zihai_data/Processed_Zihai_Data"
D = {}                                                    # (day, sex) -> list of dicts
for f in sorted(glob.glob(os.path.join(ddir, "d*_*_*.csv"))):
    day, sx, _ = os.path.basename(f)[:-4].split("_"); day = int(day[1:]); sex = "male" if sx == "M" else "female"
    s = read_slide(f); T = types(s); X, Y = bins(s); H, W = Y.max() + 1, X.max() + 1
    g = lambda m: np.bincount(Y[m] * W + X[m], minlength=H * W).reshape(H, W).astype(float)
    cd8 = T["i2"] | T["i3"] | T["i4"] | T["i6"] | T["i7"]
    area30 = len(set(zip((s["x"] // 30).astype(int), (s["y"] // 30).astype(int)))) * 9e-4
    rng = np.random.default_rng(1)
    p, c = regions(g(~s["cd3"] & (rng.random(len(cd8)) < 0.7)), g(cd8))
    D.setdefault((day, sex), []).append(dict(n=cd8.sum(), dens=cd8.sum() / area30, peri=p, core=c))
    print(f"  {os.path.basename(f):12s} CD8 {cd8.sum():6d}  area {area30:5.1f} mm2  periphery {p:6.0f}/mm2  core {c:6.0f}/mm2", flush=True)
m = lambda day, sex, k: np.nanmean([r[k] for r in D[(day, sex)]])

# ------------------------------------------------------------------ model
def model_ratio(pattern):
    R = []
    for f in glob.glob(pattern):
        d = read_log(f)
        if d is None: continue
        t = d["sum(i2)"] + d["sum(i3)"] + d["sum(i4)"] + d["sum(i6)"] + d["sum(i7)"]
        R.append(t[-1] / t[0])
    return np.array(R)

fig, axs = plt.subplots(1, 3, figsize=(13, 4.3), sharey=True)
labels = ["day 14 / day 7", "day 21 / day 14"]
for ax, kind in zip(axs, ("data: CD8 count", "data: CD8 per mm²", "model: CD8 count")):
    for i, sex in enumerate(("male", "female")):
        if kind.startswith("data"):
            k = "n" if "count" in kind else "dens"
            v = [m(14, sex, k) / m(7, sex, k), m(21, sex, k) / m(14, sex, k)]; e = [0, 0]
        else:
            r1 = model_ratio(f"runs_d7cd4/{sex}/s*/log.spparks"); r2 = model_ratio(f"runs_big/prod/{sex}/s*/log.spparks")
            v = [r1.mean() if len(r1) else np.nan, r2.mean() if len(r2) else np.nan]
            e = [r.std(ddof=1) / np.sqrt(len(r)) if len(r) > 1 else 0 for r in (r1, r2)]
        ax.bar(np.arange(2) + (i - .5) * .36, v, .36, yerr=e, capsize=3, color=COL[sex], label=sex)
        out.append(f"{kind:18s} {sex:6s}  day14/day7 {v[0]:.2f}   day21/day14 {v[1]:.2f}")
    ax.axhline(1, color="grey", lw=1, ls="--"); ax.set_xticks(range(2)); ax.set_xticklabels(labels); ax.set_title(kind)
    ax.grid(True, axis="y", color="#e6e6e6", lw=.6)
axs[0].set_ylabel("CD8+ T cells, fold change"); axs[0].legend(frameon=False)
fig.tight_layout(); os.makedirs("figures", exist_ok=True); fig.savefig("figures/CD8_ratio_bars.png", dpi=200); plt.close(fig)

# ------------------------------------------------------------------ periphery vs core
MP = {}
for sex in ("male", "female"):
    for day in (14, 21):
        vals = []
        for dd in sorted(glob.glob(f"runs_spatial21/{sex}/s*")):
            fs = sorted(glob.glob(f"{dd}/dump.sp.*"), key=lambda p: float(p.rsplit(".", 1)[1]))
            if len(fs) < 2: continue
            gd = read_dump(fs[0] if day == 14 else fs[-1])
            vals.append(regions(gd(["i1"]), gd(["i2", "i3", "i4", "i6", "i7"])))
        MP[(day, sex)] = np.nanmean(vals, 0) if vals else (np.nan, np.nan)
fig, axs = plt.subplots(1, 2, figsize=(10.5, 4.3), sharey=True)
for ax, sex in zip(axs, ("male", "female")):
    groups = [("data d14", m(14, sex, "peri"), m(14, sex, "core")), ("model d14", *MP[(14, sex)]),
              ("data d21", m(21, sex, "peri"), m(21, sex, "core")), ("model d21", *MP[(21, sex)])]
    x = np.arange(len(groups))
    ax.bar(x - .18, [g[1] for g in groups], .36, color=COL[sex], label=f"periphery (±{PERI:.0f} µm of edge)")
    ax.bar(x + .18, [g[2] for g in groups], .36, color=COL[sex], alpha=.4, label=f"core (> {PERI:.0f} µm inside)")
    ax.set_xticks(x); ax.set_xticklabels([g[0] for g in groups]); ax.set_title(sex); ax.grid(True, axis="y", color="#e6e6e6", lw=.6)
    ax.legend(frameon=False, fontsize=8)
    for g in groups: out.append(f"{sex:6s} {g[0]:9s} periphery {g[1]:6.0f}/mm2  core {g[2]:6.0f}/mm2")
axs[0].set_ylabel("CD8+ T cells per mm²"); fig.tight_layout()
fig.savefig("figures/CD8_periphery_core_bars.png", dpi=200); plt.close(fig)
open("cd8_extras_summary.txt", "w").write("\n".join(out) + "\n"); print("\n".join(out))
print("wrote figures/CD8_ratio_bars.png, figures/CD8_periphery_core_bars.png, cd8_extras_summary.txt")
