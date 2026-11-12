#!/usr/bin/env python3
"""
plot_arp_fraction.py -- fraction of AR+ cells among CD8+ T cells vs time,
standard model, male and female (mean of all finished seeds, SE band).
AR+ CD8 = i6 (mature) + i7 (terminally exhausted); AR- CD8 = i2 + i3 + i4.
Run: cd ~/CD16_NK92_project/SPPARKS_Sex_Bias && python3 ~/CD16_NK92_project/plot_arp_fraction.py
"""
import glob, os, numpy as np, matplotlib
matplotlib.use("Agg"); import matplotlib.pyplot as plt
DAYS = np.linspace(14, 21, 141)
def read(f):
    hdr = None; rows = []
    for line in open(f):
        s = line.split()
        if "sum(i1)" in s:
            if rows: break
            hdr = s; continue
        if hdr and len(s) == len(hdr):
            try: v = [float(x) for x in s]
            except ValueError: continue
            if rows and v[0] <= rows[-1][0]: break
            rows.append(v)
    if not hdr or len(rows) < 3 or rows[-1][0] < 1000: return None
    a = np.array(rows); c = lambda k: a[:, hdr.index(k)]
    arp = c("sum(i6)") + c("sum(i7)") + c("sum(i5)")
    arm = c("sum(i2)") + c("sum(i3)") + c("sum(i4)")
    return np.interp(DAYS, 14 + a[:, 0] / 144.0, 100 * arp / np.maximum(arp + arm, 1))
os.makedirs("figures", exist_ok=True)
fig, ax = plt.subplots(figsize=(6.2, 4.8))
for sex, col in (("male", "#1f4db8"), ("female", "#b81f5b")):
    R = np.array([r for r in (read(f) for f in sorted(glob.glob(f"runs_big/prod/{sex}/s*/log.spparks"))) if r is not None])
    if len(R) == 0: print("no runs for", sex); continue
    m = R.mean(0); se = R.std(0, ddof=1) / np.sqrt(len(R))
    ax.plot(DAYS, m, color=col, lw=2.4, label=sex)
    ax.fill_between(DAYS, m - se, m + se, color=col, alpha=.2, lw=0)
    print(f"{sex:6s} n={len(R)}  AR+ % of CD8 at day 14/16/18/21: " + "  ".join(f"{m[i]:.1f}" for i in (0, 40, 80, 140)))
ax.set_xlabel("Day"); ax.set_ylabel("AR+ cells (% of CD8+ T cells)")
ax.set_title("AR+ fraction of CD8+ T cells"); ax.set_ylim(bottom=0)
ax.legend(frameon=False); ax.grid(True, color="#e6e6e6", lw=.6); fig.tight_layout()
fig.savefig("figures/ARp_fraction_CD8.png", dpi=200); print("wrote figures/ARp_fraction_CD8.png")
