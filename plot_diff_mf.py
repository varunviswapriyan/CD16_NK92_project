#!/usr/bin/env python3
"""
plot_diff_mf.py -- same style as the earlier diffusion figure (front^2 vs time +
linear fit), now female-like (no AR+ cells) vs male-like (AR+ cells internalizing).
Uses the simplified-model outputs already made by run_fb_diff.sh (diff_arp/).
Run: cd ~/CD16_NK92_project/SPPARKS_Sex_Bias && python3 ~/CD16_NK92_project/plot_diff_mf.py [male_density]
"""
import sys, csv, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
MALE = sys.argv[1] if len(sys.argv) > 1 else "0.03"   # AR+ cells per IFN-g pixel for "male"
def load(a):
    R = {}
    for r in csv.DictReader(open(f"diff_arp/arp_{a}_front.csv")):
        R.setdefault(r["run"], []).append((float(r["time_h"]), float(r["front_um2"])))
    return np.array(sorted(R["bound"])), np.array(sorted(R["bare"]))
def window(x):
    mx = x.max(); m = (x > 0.05*mx) & (x < 0.55*mx)
    return m if m.sum() >= 3 else np.ones_like(x, bool)
fig, ax = plt.subplots(figsize=(6.4, 5))
for a, col, lab in (("0", "#b81f5b", "female (no AR+ T cells)"),
                    (MALE, "#1f4db8", f"male (AR+ T cells, {MALE} per 30 µm pixel)")):
    b, z = load(a)
    t, x = b[:, 0], b[:, 1]; m = window(x)
    sl, ic = np.polyfit(t[m], x[m], 1)
    zw = window(z[:, 1]); zsl = np.polyfit(z[zw, 0], z[zw, 1], 1)[0]
    d = 20.0 * sl / zsl
    if a == "0":
        ax.plot(t, x, "-", color=col, lw=1.2, alpha=.45)
        ax.plot(t, sl*t + ic, "--", color=col, lw=1.5, label=f"{lab}: linear fit, D_eff = {d:.2f} µm²/s")
        ax.plot([],[],"-",color=col,lw=1.2,alpha=.45,label=f"{lab}: simulation")
        print(f"{lab:40s} D_eff = {d:.3f} um^2/s")
    else:  # internalization is a sink: front stalls, a single D_eff does not describe it
        ax.plot(t, x, "-", color=col, lw=2, label=f"{lab}: front stalls")
        e = t <= 24; se = np.polyfit(t[e], x[e], 1)[0]
        print(f"{lab:40s} front^2 at end = {x[-1]:.0f} um^2 vs female {fem_end:.0f}  (early slope {se:.0f} um^2/h)")
    if a == "0": fem_end = x[-1]
ax.set_xlabel("Time (hours)"); ax.set_ylabel(r"Front position$^2$, $x^2$ ($\mu m^2$)")
ax.set_title("IFN-γ spread into the tumour"); ax.legend(frameon=False, fontsize=9)
ax.grid(True, color="#e6e6e6", lw=.6); fig.tight_layout()
fig.savefig("figures/Diffusion_male_vs_female.png", dpi=200)
print("wrote figures/Diffusion_male_vs_female.png")
