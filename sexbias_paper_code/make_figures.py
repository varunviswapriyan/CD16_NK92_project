#!/usr/bin/env python3
"""
Make the manuscript figures from sex_bias_sim_v7 / ifng_effective_diffusion output.

  Fig_cancer_timecourse.png  cancer fold change (and free IFNG) vs time,
                             male vs female, mean +/- SD over seeds
                             (+ male without sequestration if present)
  Fig_spatial_cells.png      cell positions, male vs female, day 14 / 17.5 / 21
  Fig_spatial_ifng.png       free IFNG maps at the same times
  Fig_effective_diffusion.png IFNG profiles + front^2 vs t (thesis Fig 4.8 remake)

Usage:  python3 make_figures.py RESULTS_DIR
RESULTS_DIR must contain runs/<condition>_s<seed>_timecourse.csv,
runs/<condition>_s<seed>_snap_tNNNN.txt and diffusion/ifng_diffusion_*.csv
(run_all.sh lays it out this way).
"""
import glob
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

RES = sys.argv[1] if len(sys.argv) > 1 else "results"
RUNS = os.path.join(RES, "runs")
OUT = os.path.join(RES, "figures")
os.makedirs(OUT, exist_ok=True)

COND = {  # label, colour, linestyle
    "male": ("Male", "#2a5fbf", "-"),
    "female": ("Female", "#c2338a", "-"),
    "male_noseq": ("Male, no AR+ internalization", "#2a5fbf", "--"),
}
plt.rcParams.update({"font.size": 10, "axes.spines.top": False,
                     "axes.spines.right": False, "axes.grid": True,
                     "grid.color": "#e3e3e3", "grid.linewidth": 0.6})


def load_tc(path):
    return np.genfromtxt(path, delimiter=",", names=True)


# ---------------------------------------------------------------- time course
def fig_timecourse():
    data = {}
    for cond in COND:
        files = sorted(glob.glob(os.path.join(RUNS, f"{cond}_s*_timecourse.csv")))
        if not files:
            continue
        tcs = [load_tc(f) for f in files]
        n = min(len(t) for t in tcs)
        days = 14 + tcs[0]["time_days"][:n]
        fold = np.array([t["cancer"][:n] / t["cancer"][0] for t in tcs])
        ifng = np.array([t["free_ifng"][:n] for t in tcs])
        data[cond] = (days, fold, ifng, len(files))
    if not data:
        print("no time-course files found"); return

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for cond, (days, fold, ifng, nseed) in data.items():
        lab, col, ls = COND[cond]
        for ax, y in ((axes[0], fold), (axes[1], ifng)):
            m, s = y.mean(0), y.std(0)
            ax.plot(days, m, color=col, ls=ls, lw=2, label=f"{lab} (n={nseed})")
            if nseed > 1:
                ax.fill_between(days, m - s, m + s, color=col, alpha=0.15, lw=0)
    axes[0].axhline(1, color="#999", lw=0.8)
    axes[0].set(xlabel="Day after tumour implantation", ylabel="Cancer cells (fold change from day 14)",
                title="Tumour growth")
    axes[1].set(xlabel="Day after tumour implantation", ylabel="Total free IFN-γ (molecules)",
                title="Free IFN-γ in the tissue")
    axes[0].legend(frameon=False, fontsize=9)
    for ax in axes:
        ax.set_xlim(14, 21)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "Fig_cancer_timecourse.png"), dpi=300)
    plt.close(fig)

    # summary table for the text
    with open(os.path.join(OUT, "summary_day21.txt"), "w") as fh:
        fh.write("condition  n_seeds  fold_change_day21 (mean +/- SD)  free_IFNG_day21 (mean)\n")
        for cond, (days, fold, ifng, nseed) in data.items():
            fh.write(f"{cond:11s} {nseed:3d}   {fold[:, -1].mean():.4f} +/- {fold[:, -1].std():.4f}"
                     f"   {ifng[:, -1].mean():.4g}\n")
        if "male" in data and "female" in data:
            r = data["male"][1][:, -1].mean() / data["female"][1][:, -1].mean()
            fh.write(f"\nmale/female fold-change ratio at day 21: {r:.4f}\n")
        if "male_noseq" in data and "female" in data:
            r = data["male_noseq"][1][:, -1].mean() / data["female"][1][:, -1].mean()
            fh.write(f"male(no internalization)/female ratio at day 21: {r:.4f}\n")
    print(open(os.path.join(OUT, "summary_day21.txt")).read())


# ------------------------------------------------------------------- spatial
def load_snap(path):
    a = np.loadtxt(path, comments="#")
    return a  # x y cancer Tp Tm Te ARPE ARTe CD4 ARCD4 free_ifng bound_frac


def snap_files(cond):
    first = sorted(glob.glob(os.path.join(RUNS, f"{cond}_s*_snap_t*.txt")))
    if not first:
        return []
    seed_tag = os.path.basename(first[0]).split("_snap_")[0]
    return sorted(glob.glob(os.path.join(RUNS, f"{seed_tag}_snap_t*.txt")))


def fig_spatial():
    conds = [c for c in ("male", "female") if snap_files(c)]
    if not conds:
        print("no snapshot files found"); return
    ntime = max(len(snap_files(c)) for c in conds)

    # --- cells
    fig, axes = plt.subplots(len(conds), ntime, figsize=(3.3 * ntime, 3.4 * len(conds)),
                             squeeze=False)
    layers = [  # column(s), label, colour, size
        ((2,), "Cancer", "#d9d9d9", 0.15),
        ((3, 4), "AR− CD8 (progenitor + mature)", "#2a5fbf", 0.8),
        ((6,), "AR+ progenitor-exhausted CD8", "#e8740c", 1.6),
        ((5, 7), "Exhausted CD8", "#222222", 0.8),
    ]
    for r, cond in enumerate(conds):
        for c, f in enumerate(snap_files(cond)):
            a = load_snap(f)
            t_units = int(os.path.basename(f).split("_t")[-1].split(".")[0])
            ax = axes[r][c]
            for cols, lab, col, size in layers:
                m = a[:, list(cols)].sum(1) > 0
                ax.scatter(a[m, 0] * 0.01, a[m, 1] * 0.01, s=size, c=col, lw=0,
                           label=lab, rasterized=True)
            ax.set_aspect("equal"); ax.set_xlim(0, 6); ax.set_ylim(0, 6)
            ax.grid(False)
            ax.set_title(f"{COND[cond][0]}, day {14 + t_units / 144:.1f}")
            ax.set_xlabel("mm"); ax.set_ylabel("mm")
    h, l = axes[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=4, frameon=False, markerscale=8)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(os.path.join(OUT, "Fig_spatial_cells.png"), dpi=300)
    plt.close(fig)

    # --- IFNG (coarse 30 um pixels, shared log colour scale)
    grids, vmax = {}, 0.0
    for cond in conds:
        for f in snap_files(cond):
            a = load_snap(f)
            g = np.zeros((200, 200))
            xi, yi = (a[:, 0] // 3).astype(int), (a[:, 1] // 3).astype(int)
            g[yi, xi] = np.maximum(g[yi, xi], a[:, 10])
            grids[(cond, f)] = g
            vmax = max(vmax, g.max())
    if vmax <= 0:
        print("IFNG all zero in snapshots"); return
    norm = LogNorm(vmin=vmax * 1e-4, vmax=vmax)
    fig, axes = plt.subplots(len(conds), ntime, figsize=(3.3 * ntime, 3.4 * len(conds)),
                             squeeze=False)
    for r, cond in enumerate(conds):
        for c, f in enumerate(snap_files(cond)):
            t_units = int(os.path.basename(f).split("_t")[-1].split(".")[0])
            g = np.ma.masked_less_equal(grids[(cond, f)], vmax * 1e-4)
            ax = axes[r][c]
            ax.set_facecolor("white")
            im = ax.imshow(g, origin="lower", extent=(0, 6, 0, 6), cmap="Purples", norm=norm)
            ax.grid(False)
            ax.set_title(f"{COND[cond][0]}, day {14 + t_units / 144:.1f}")
            ax.set_xlabel("mm"); ax.set_ylabel("mm")
    cb = fig.colorbar(im, ax=axes, shrink=0.8)
    cb.set_label("Free IFN-γ per 10 µm site (molecules)")
    fig.savefig(os.path.join(OUT, "Fig_spatial_ifng.png"), dpi=300, bbox_inches="tight")
    plt.close(fig)


# ------------------------------------------------------- effective diffusion
def fig_diffusion():
    pf = os.path.join(RES, "diffusion", "ifng_diffusion_profiles.csv")
    ff = os.path.join(RES, "diffusion", "ifng_diffusion_front.csv")
    if not (os.path.exists(pf) and os.path.exists(ff)):
        print("no diffusion output found"); return
    prof = np.genfromtxt(pf, delimiter=",", names=True)
    fr = np.genfromtxt(ff, delimiter=",", names=True, dtype=None, encoding=None)

    def fit(run, lo=300.0, hi=2500.0):
        sel = (fr["run"] == run) & (fr["front_um"] >= lo) & (fr["front_um"] <= hi)
        t = fr["time_h"][sel] * 3600.0
        x2 = fr["front_um2"][sel]
        slope, icpt = np.polyfit(t, x2, 1)
        ll = np.polyfit(np.log(t), np.log(x2), 1)[0]
        return t, x2, slope, icpt, ll

    tb, xb, sb, ib, llb = fit("bound")
    t0, x0, s0, i0, ll0 = fit("bare")
    D = 20.0
    D_eff = D * sb / s0                       # calibrated against bare diffusion
    D_pred = D / (1 + 9000.0 / 100.0)         # D / (1 + R0/K_D)

    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    times = np.unique(prof["time_h"])
    cmap = plt.get_cmap("Purples")
    for k, th in enumerate(times):
        m = (prof["time_h"] == th) & (prof["x_um"] >= 0)
        axes[0].plot(prof["x_um"][m], prof["free_ifng_per_site"][m],
                     color=cmap(0.35 + 0.65 * k / max(1, len(times) - 1)), lw=2,
                     label=f"{th / 24:.0f} d")
    axes[0].set(xlim=(0, 1500), xlabel="Distance from IFN-γ source (µm)",
                ylabel="Free IFN-γ per 10 µm site", title="IFN-γ spreading into a uniform tumour nest")
    axes[0].legend(frameon=False, title="time", fontsize=8)

    ax = axes[1]
    ax.plot(tb / 3600, xb / 1e6, "o", ms=3, color="#5b3c9e", label="simulation")
    tt = np.linspace(tb.min(), tb.max(), 50)
    ax.plot(tt / 3600, (sb * tt + ib) / 1e6, "-", color="#222", lw=1.5, label="linear fit")
    ax.set(xlabel="Time (h)", ylabel="Front position² (mm²)", title="IFN-γ front (1% of maximum)")
    ax.text(0.03, 0.97,
            f"log–log slope = {llb:.2f}  (1 = diffusive)\n"
            f"$D_{{eff}}$ = {D_eff:.2f} µm²/s  ({D_eff / D:.3f} × D)\n"
            f"$D/(1+R_0/K_D)$ = {D_pred:.2f} µm²/s  ({D_pred / D:.3f} × D)",
            transform=ax.transAxes, va="top", fontsize=9)
    ax.legend(frameon=False, loc="lower right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "Fig_effective_diffusion.png"), dpi=300)
    plt.close(fig)
    with open(os.path.join(OUT, "effective_diffusion.txt"), "w") as fh:
        fh.write(f"bound run: x^2 slope = {sb:.4f} um^2/s, log-log slope = {llb:.3f}\n")
        fh.write(f"bare  run: x^2 slope = {s0:.4f} um^2/s, log-log slope = {ll0:.3f} (D = {D})\n")
        fh.write(f"D_eff = D * slope_bound/slope_bare = {D_eff:.4f} um^2/s = {D_eff / D:.4f} D\n")
        fh.write(f"analytic D/(1+R0/K_D) = {D_pred:.4f} um^2/s = {D_pred / D:.4f} D\n")
    print(open(os.path.join(OUT, "effective_diffusion.txt")).read())


if __name__ == "__main__":
    fig_timecourse()
    fig_spatial()
    fig_diffusion()
    print("figures written to", OUT)
