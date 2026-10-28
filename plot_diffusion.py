#!/usr/bin/env python3
"""
plot_diffusion.py FRONT.csv
Plots IFN-g front position^2 vs time (the MSD-style plot, thesis Fig 4.8b) for
the bound and bare-diffusion runs, fits each with a line, and reports the
effective diffusion constant D_eff = D_bare * (slope_bound / slope_bare).
front.csv columns: run,time_h,front_um,front_um2   (run = 'bound' or 'bare')
"""
import sys, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

path = sys.argv[1] if len(sys.argv) > 1 else "diffusion_front.csv"
D_BARE = 20.0  # um^2/s, microscopic IFN-g diffusion constant

rows = {}
with open(path) as f:
    for r in csv.DictReader(f):
        rows.setdefault(r["run"], []).append((float(r["time_h"]), float(r["front_um2"])))

fig, ax = plt.subplots(figsize=(6, 5))
slopes = {}
colors = {"bound": "#b81f5b", "bare": "#1f4db8"}
for run, pts in rows.items():
    pts = np.array(sorted(pts))
    t, x2 = pts[:, 0], pts[:, 1]
    sl, icpt = np.polyfit(t, x2, 1)
    slopes[run] = sl
    ax.plot(t, x2, "o", ms=4, color=colors.get(run, "#555"), label=f"{run} (slope={sl:.3g})")
    ax.plot(t, sl * t + icpt, "-", color=colors.get(run, "#555"), lw=1.5)

ax.set(xlabel="Time (hours)", ylabel=r"Front position$^2$ ($\mu m^2$)",
       title="IFN-$\\gamma$ spread into the tumour (effective diffusion)")
ax.legend(frameon=False); ax.grid(True, color="#e6e6e6", lw=.6)

if "bound" in slopes and "bare" in slopes and slopes["bare"]:
    d_eff = D_BARE * slopes["bound"] / slopes["bare"]
    ax.text(.05, .95, f"D_eff ≈ {d_eff:.2f} µm²/s\n(bare D = {D_BARE:.0f})",
            transform=ax.transAxes, va="top",
            bbox=dict(boxstyle="round", fc="white", ec="#ccc"))
    print(f"slope_bound={slopes['bound']:.4g}  slope_bare={slopes['bare']:.4g}")
    print(f"D_eff = {d_eff:.3f} um^2/s   (thesis old model gave 0.76)")

fig.tight_layout()
fig.savefig("Fig_effective_diffusion.png", dpi=200)
print("wrote Fig_effective_diffusion.png")
