#!/usr/bin/env python3
"""
zihai_analysis.py -- compare the model with the Zihai mIF slides (day 7, 14, 21; male and female).
Run from ~/CD16_NK92_project/SPPARKS_Sex_Bias.

  python3 zihai_analysis.py data DATA_DIR      per-slide numbers + infiltration profiles from the slides
  python3 zihai_analysis.py figures            model vs data figures + summary (needs `data` first)

Outputs (in SPPARKS_Sex_Bias/):
  zihai_slide_metrics.csv            one row per slide: counts, AR+ % of CD8, CD4/CD8, densities
  zihai_profiles.npz                 slide infiltration profiles (cached)
  figures/ARp_fraction_model_vs_data.png
  figures/Infiltration_model_vs_data.png
  figures/Day7_to_14_fold.png
  zihai_summary.txt
Same definitions on both sides: cancer/tissue mask = 31x31-site (310 um) box average of
"site has a CD3- cell (data) / cancer cell (model)" > 0.5; distance = Euclidean distance to the mask
edge, + inside, - outside; CD8 density = CD8 cells per mm^2 of lattice area in each 100 um band.
"""
import sys, os, glob, csv
import numpy as np
from convert_zihai import read_slide, types, bins

UM = 10.0; PIX_MM2 = (UM / 1000) ** 2
BANDS = np.arange(-500, 801, 100)
try:
    from scipy.ndimage import distance_transform_edt as _edt
except Exception:
    _edt = None


def box_blur(a, k):
    c = np.cumsum(np.cumsum(np.pad(a, ((k + 1, k), (k + 1, k))), 0), 1)
    return (c[2 * k + 1:, 2 * k + 1:] - c[:-2 * k - 1, 2 * k + 1:] - c[2 * k + 1:, :-2 * k - 1] + c[:-2 * k - 1, :-2 * k - 1]) / (2 * k + 1) ** 2


def edt(m):
    if _edt is not None: return _edt(m)
    d = np.where(m, np.inf, 0.0)                     # fallback: 4-neighbour distance
    while True:
        P = np.pad(d, 1, constant_values=0.0)
        n = np.minimum(np.minimum(np.minimum(P[:-2, 1:-1], P[2:, 1:-1]), P[1:-1, :-2]), P[1:-1, 2:]) + 1
        new = np.where(m, np.minimum(d, n), 0.0)
        if np.array_equal(new, d): return d
        d = new


def profile(cancer_grid, cd8_grid):
    mask = box_blur((cancer_grid > 0).astype(float), 15) > 0.5
    sd = (edt(mask) - edt(~mask)) * UM
    out = []
    for a, b in zip(BANDS[:-1], BANDS[1:]):
        sel = (sd >= a) & (sd < b)
        out.append(cd8_grid[sel].sum() / (sel.sum() * PIX_MM2) if sel.sum() > 0 else np.nan)
    return np.array(out), mask.sum() * PIX_MM2


def cmd_data(ddir):
    files = sorted(glob.glob(os.path.join(ddir, "d*_*_*.csv")))
    if not files: sys.exit(f"no slide csv files in {ddir}")
    rows, prof = [], {}
    for f in files:
        name = os.path.basename(f)[:-4]; day, sex, rep = name.split("_"); day = int(day[1:])
        s = read_slide(f); T = types(s); X, Y = bins(s)
        H, W = Y.max() + 1, X.max() + 1
        g = lambda m: np.bincount(Y[m] * W + X[m], minlength=H * W).reshape(H, W).astype(float)
        cd8 = T["i2"] | T["i3"] | T["i4"] | T["i6"] | T["i7"]; arp = T["i6"] | T["i7"]; cd4 = T["i13"] | T["i14"]
        p, mask_mm2 = profile(g(~s["cd3"]), g(cd8))
        prof[name] = p
        tissue_mm2 = (box_blur((g(np.ones(len(cd8), bool)) > 0).astype(float), 5) > 0.2).sum() * PIX_MM2
        rows.append(dict(slide=name, day=day, sex="male" if sex == "M" else "female", cells=len(cd8),
                         CD8=int(cd8.sum()), ARp_CD8=int(arp.sum()), CD4=int(cd4.sum()), ARp_CD4=int(T["i14"].sum()),
                         ARp_pct_CD8=round(100 * arp.sum() / max(cd8.sum(), 1), 2),
                         CD4_per_CD8=round(cd4.sum() / max(cd8.sum(), 1), 3),
                         tissue_mm2=round(tissue_mm2, 2), tumour_mask_mm2=round(mask_mm2, 2),
                         CD8_per_mm2=round(cd8.sum() / tissue_mm2, 1), ARp_CD8_per_mm2=round(arp.sum() / tissue_mm2, 2)))
        print(f"  {name:8s} CD8 {cd8.sum():6d}  AR+ {100 * arp.sum() / max(cd8.sum(), 1):5.1f}%  CD4/CD8 {cd4.sum() / max(cd8.sum(), 1):4.2f}  "
              f"tissue {tissue_mm2:6.1f} mm2  CD8 {cd8.sum() / tissue_mm2:6.1f}/mm2", flush=True)
    with open("zihai_slide_metrics.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    np.savez("zihai_profiles.npz", **prof)
    print("wrote zihai_slide_metrics.csv, zihai_profiles.npz")


# ----------------------------------------------------------------- model side
def read_log(f):
    hdr, rows = None, []
    for line in open(f):
        p = line.split()
        if "sum(i1)" in p:
            if rows: break
            hdr = p; continue
        if hdr and len(p) == len(hdr):
            try: v = [float(x) for x in p]
            except ValueError: continue
            if rows and v[0] <= rows[-1][0]: break
            rows.append(v)
    if not hdr or len(rows) < 3 or rows[-1][0] < 1000: return None
    a = np.array(rows); return {k: a[:, i] for i, k in enumerate(hdr)}


def model_arp(pattern, day0):
    days = np.linspace(day0, day0 + 7, 141); R, F = [], []
    for f in sorted(glob.glob(pattern)):
        d = read_log(f)
        if d is None: continue
        arp = d["sum(i5)"] + d["sum(i6)"] + d["sum(i7)"]; arm = d["sum(i2)"] + d["sum(i3)"] + d["sum(i4)"]
        R.append(np.interp(days, day0 + d["Time"] / 144.0 if "Time" in d else day0 + d[list(d)[0]] / 144.0,
                           100 * arp / np.maximum(arp + arm, 1)))
        F.append(d["sum(i1)"][-1] / d["sum(i1)"][0])
    return days, np.array(R), np.array(F)


def read_dump(p):
    cols, data, inv = None, [], False
    for line in open(p):
        s = line.split()
        if not s: continue
        if s[0] == "Values": inv = True; continue
        if cols is None and s[0] == "id": cols = [c for c in s if c != "values"]; continue
        if inv and cols and len(s) >= len(cols):
            try: data.append([float(x) for x in s[:len(cols)]])
            except ValueError: pass
    a = np.array(data); d = {c: a[:, i] for i, c in enumerate(cols)}
    n = int(round(np.sqrt(len(a)))); ix = np.clip(d["x"].astype(int), 0, n - 1); iy = np.clip(d["y"].astype(int), 0, n - 1)
    def g(keys):
        G = np.zeros((n, n))
        for k in keys: np.add.at(G, (iy, ix), d[k])
        return G
    return g


def cmd_figures():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    os.makedirs("figures", exist_ok=True)
    rows = list(csv.DictReader(open("zihai_slide_metrics.csv"))); prof = dict(np.load("zihai_profiles.npz"))
    col = {"male": "#1f4db8", "female": "#b81f5b"}; lines = []
    # ---- AR+ fraction
    fig, ax = plt.subplots(figsize=(6.6, 4.8))
    for sex in ("male", "female"):
        for pat, d0, ls in ((f"runs_d7cd4/{sex}/s*/log.spparks", 7, "--"), (f"runs_big/prod/{sex}/s*/log.spparks", 14, "-")):
            days, R, F = model_arp(pat, d0)
            if len(R):
                m = R.mean(0); ax.plot(days, m, ls, color=col[sex], lw=2, label=f"model {sex}, day {d0}-{d0 + 7} (n={len(R)})")
                lines.append(f"model {sex} day {d0}->{d0 + 7}: AR+ % of CD8 {m[0]:.1f} -> {m[-1]:.1f} (n={len(R)})")
        for day in (7, 14, 21):
            v = [float(r["ARp_pct_CD8"]) for r in rows if r["sex"] == sex and int(r["day"]) == day]
            if v:
                ax.errorbar(day + (0.15 if sex == "female" else -0.15), np.mean(v), yerr=np.std(v, ddof=1) / np.sqrt(len(v)) if len(v) > 1 else 0,
                            fmt="o", color=col[sex], ms=7, capsize=4, label=f"data {sex}" if day == 7 else None)
                ax.scatter([day + (0.15 if sex == "female" else -0.15)] * len(v), v, s=10, color=col[sex], alpha=.35)
                lines.append(f"data  {sex} day {day}: AR+ % of CD8 {np.mean(v):.1f} (slides {', '.join(f'{x:.1f}' for x in v)})")
    ax.set_xlabel("Day"); ax.set_ylabel("AR+ cells (% of CD8+ T cells)"); ax.set_title("AR+ fraction of CD8+ T cells: model vs mIF")
    ax.set_xticks([7, 14, 21]); ax.set_ylim(bottom=0); ax.grid(True, color="#e6e6e6", lw=.6); ax.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig("figures/ARp_fraction_model_vs_data.png", dpi=200); plt.close(fig)
    # ---- infiltration profiles
    mid = (BANDS[:-1] + BANDS[1:]) / 2
    fig, axs = plt.subplots(1, 2, figsize=(11, 4.4), sharey=True)
    for ax, sex in zip(axs, ("male", "female")):
        for day, ls in ((14, "--"), (21, "-")):
            P = [prof[r["slide"]] for r in rows if r["sex"] == sex and int(r["day"]) == day]
            if P:
                ax.plot(mid, np.nanmean(P, 0), ls, color="k", lw=2, label=f"data day {day} (n={len(P)})")
        for day, ls in ((14, "--"), (21, "-")):
            P = []
            for dd in sorted(glob.glob(f"runs_spatial21/{sex}/s*")):
                fs = sorted(glob.glob(f"{dd}/dump.sp.*"), key=lambda p: float(p.rsplit(".", 1)[1]))
                if len(fs) < 2: continue
                g = read_dump(fs[0] if day == 14 else fs[-1])
                P.append(profile(g(["i1"]), g(["i2", "i3", "i4", "i6", "i7"]))[0])
            if P:
                ax.plot(mid, np.nanmean(P, 0), ls, color=col[sex], lw=2, label=f"model day {day} (n={len(P)})")
                lines.append(f"model {sex} day {day} CD8/mm2 by 100-um band (-500..800 um): " + ", ".join(f"{v:.0f}" for v in np.nanmean(P, 0)))
            Pd = [prof[r["slide"]] for r in rows if r["sex"] == sex and int(r["day"]) == day]
            if Pd: lines.append(f"data  {sex} day {day} CD8/mm2 by 100-um band (-500..800 um): " + ", ".join(f"{v:.0f}" for v in np.nanmean(Pd, 0)))
        ax.axvline(0, color="grey", lw=1); ax.set_title(sex); ax.set_xlabel("Distance from tumour edge (µm; + = inside)")
        ax.grid(True, color="#e6e6e6", lw=.6); ax.legend(frameon=False, fontsize=8)
    axs[0].set_ylabel("CD8+ T cells per mm²"); fig.tight_layout()
    fig.savefig("figures/Infiltration_model_vs_data.png", dpi=200); plt.close(fig)
    # ---- day 7 -> 14 fold
    res = {}
    for sex in ("male", "female"):
        _, _, F = model_arp(f"runs_d7cd4/{sex}/s*/log.spparks", 7); _, _, F2 = model_arp(f"runs_big/prod/{sex}/s*/log.spparks", 14)
        res[sex] = (F, F2)
        if len(F): lines.append(f"model {sex} cancer fold day 7->14: {F.mean():.3f} +/- {F.std(ddof=1) / np.sqrt(len(F)) if len(F) > 1 else 0:.3f} (n={len(F)})")
        if len(F2): lines.append(f"model {sex} cancer fold day 14->21: {F2.mean():.3f} (n={len(F2)})")
    if len(res["male"][0]) and len(res["female"][0]):
        r7 = res["male"][0].mean() / res["female"][0].mean(); lines.append(f"male/female fold ratio day 7->14: {r7:.2f}")
        if len(res["male"][1]) and len(res["female"][1]):
            lines.append(f"male/female fold ratio day 14->21: {res['male'][1].mean() / res['female'][1].mean():.2f}")
        fig, ax = plt.subplots(figsize=(5.4, 4.4)); xs = np.arange(2); w = .35
        for i, sex in enumerate(("male", "female")):
            vals = [res[sex][0].mean(), res[sex][1].mean() if len(res[sex][1]) else np.nan]
            ax.bar(xs + (i - .5) * w, vals, w, color=col[sex], label=sex)
        ax.set_xticks(xs); ax.set_xticklabels(["day 7 → 14", "day 14 → 21"]); ax.set_ylabel("Cancer cell fold change (model)")
        ax.legend(frameon=False); ax.grid(True, axis="y", color="#e6e6e6", lw=.6); fig.tight_layout()
        fig.savefig("figures/Day7_to_14_fold.png", dpi=200); plt.close(fig)
    open("zihai_summary.txt", "w").write("\n".join(lines) + "\n"); print("\n".join(lines))
    print("wrote figures/ARp_fraction_model_vs_data.png, figures/Infiltration_model_vs_data.png, figures/Day7_to_14_fold.png, zihai_summary.txt")


if __name__ == "__main__":
    if len(sys.argv) >= 3 and sys.argv[1] == "data": cmd_data(sys.argv[2])
    elif len(sys.argv) >= 2 and sys.argv[1] == "figures": cmd_figures()
    else: print(__doc__)
