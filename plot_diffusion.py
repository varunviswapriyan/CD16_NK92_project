#!/usr/bin/env python3
"""
plot_diffusion.py FRONT.csv  ->  figures/Figure5.png
Single-panel MSD figure (matches Giuseppe thesis Fig 4.8): the IFN-g front's
squared position vs time for the with-receptor run. A straight line
(front^2 proportional to t) is the signature of diffusion; its slope gives the
effective diffusion constant D_eff. The receptor-free run is read only to
calibrate slope -> D (not plotted); D_eff is printed for the caption.
front.csv columns: run,time_h,front_um,front_um2   (run = 'bound' or 'bare')
"""
import sys, csv, os
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
    mx = x2.max(); m = (x2 > 0.05*mx) & (x2 < 0.55*mx)
    if m.sum() < 3: m = np.ones_like(x2, bool)
    sl, ic = np.polyfit(t[m], x2[m], 1); return sl, ic

b = np.array(sorted(rows["bound"])); bt, bx = b[:, 0], b[:, 1]
bsl, bic = linfit(bt, bx)
d_eff = None
if "bare" in rows:
    r = np.array(sorted(rows["bare"])); rsl, _ = linfit(r[:, 0], r[:, 1])
    if rsl: d_eff = D_BARE * bsl / rsl

fig, ax = plt.subplots(figsize=(6.2, 5))
ax.plot(bt, bx, "o", ms=3, color="#b81f5b", alpha=.5, label="simulation")
ax.plot(bt, bsl*bt + bic, "-", color="k", lw=2, label="linear fit")
ax.set_xlabel("Time (hours)")
ax.set_ylabel(r"Front position$^2$, $x^2$ ($\mu m^2$)")
ax.set_title("IFN-γ spread into the tumour")
ax.legend(frameon=False); ax.grid(True, color="#e6e6e6", lw=.6)
fig.tight_layout()
os.makedirs("figures", exist_ok=True)
fig.savefig("figures/Figure5.png", dpi=200)
print("wrote figures/Figure5.png")
if d_eff is not None:
    print(f"D_eff = {d_eff:.2f} um^2/s  (for the caption; bare D = {D_BARE:.0f})")
