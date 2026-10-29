#!/usr/bin/env python3
"""
plot_diffusion.py FRONT.csv
Effective-diffusion figure. The bare-diffusion front races out and SATURATES
against the box edge, so a straight line through the whole bare curve is wrong
(that was the "weird" fit). Here each run is fit ONLY over its linear diffusive
regime (front^2 between 5% and 55% of its own max), shown as two clean panels.
D_eff = D_bare * slope_bound / slope_bare.
front.csv columns: run,time_h,front_um,front_um2   (run = 'bound' or 'bare')
"""
import sys, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

path = sys.argv[1] if len(sys.argv) > 1 else "diffusion_front.csv"
D_BARE = 20.0

rows = {}
with open(path) as f:
    for r in csv.DictReader(f):
        rows.setdefault(r["run"], []).append((float(r["time_h"]), float(r["front_um2"])))

def linfit(t, x2):
    """fit the linear diffusive regime: 5%..55% of the run's own max front^2."""
    mx = x2.max()
    m = (x2 > 0.05*mx) & (x2 < 0.55*mx)
    if m.sum() < 3:
        m = np.ones_like(x2, bool)
    sl, ic = np.polyfit(t[m], x2[m], 1)
    return sl, ic, m

fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.6))
res = {}
for run, ax, col in (("bound", a1, "#b81f5b"), ("bare", a2, "#1f4db8")):
    if run not in rows:
        ax.text(.5,.5,f"no '{run}' data",ha="center",va="center",transform=ax.transAxes); continue
    pts = np.array(sorted(rows[run])); t, x2 = pts[:,0], pts[:,1]
    sl, ic, m = linfit(t, x2); res[run] = sl
    ax.plot(t, x2, "o", ms=3, color=col, alpha=.45, label="data")
    ax.plot(t[m], sl*t[m]+ic, "-", color="k", lw=1.8, label="linear-regime fit")
    ax.set(xlabel="Time (hours)", ylabel=r"Front position$^2$ ($\mu m^2$)")
    ax.legend(frameon=False, fontsize=9); ax.grid(True, color="#e6e6e6", lw=.6)

a1.set_title("Bound (IFN-γ binds receptors)")
a2.set_title("Bare diffusion (D = 20, calibration)")

if "bound" in res and "bare" in res and res["bare"]:
    d_eff = D_BARE * res["bound"] / res["bare"]
    fig.suptitle(f"IFN-γ spread into the tumour — effective diffusion "
                 f"D_eff ≈ {d_eff:.2f} µm²/s   (bare D = 20 µm²/s; reduced {res['bare']/res['bound']:.0f}×)",
                 fontsize=12)
    print(f"slope_bound={res['bound']:.4g}  slope_bare={res['bare']:.4g}  D_eff={d_eff:.3f}")
else:
    fig.suptitle("IFN-γ spread into the tumour — effective diffusion", fontsize=12)

fig.tight_layout()
fig.savefig("Fig_effective_diffusion.png", dpi=200)
print("wrote Fig_effective_diffusion.png")
