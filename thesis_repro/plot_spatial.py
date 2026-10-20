#!/usr/bin/env python3
"""
plot_spatial.py RESDIR

Reads the SPPARKS spatial dumps in RESDIR/<cond>/dump.2D_model.<t> (written by
run_spatial.sh) and draws two figures, male on top / female below, at
day 14 (t=0), 17.5 (t=504), 21 (t=1008):

  Fig_spatial_cells.png  -- cancer density (red) with CD8 T cells overlaid
                            (AR- = blue, AR+ = green)
  Fig_spatial_ifng.png   -- free IFN-g concentration map

Dump columns (from run_spatial.sh): id x y i1 i2 i3 i4 i5 i6 i7 d1
  i1=cancer  i2,i3,i4 = AR- CD8 (Tp,Tm,Te)  i5,i6,i7 = AR+ CD8   d1 = free IFN-g
"""
import sys, os, glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

L = 600
# SPPARKS names the dumps by index: .0 = day 14, .1 = day 17.5, .2 = day 21
TIMES = [(0, "Day 14"), (1, "Day 17.5"), (2, "Day 21")]
CONDS = ["male", "female"]


def read_dump(path):
    """Parse a SPPARKS native 'dump sites' file:
       ... header ...
       id x y i1 i2 i3 i4 i5 i6 i7 d1 values
       ... box bounds ...
       Values
       <blank>
       <id x y i1 ... d1>   one row per site
    """
    cols = None
    data = []
    in_values = False
    with open(path) as f:
        for line in f:
            s = line.split()
            if not s:
                continue
            if s[0] == "Values":
                in_values = True
                continue
            if cols is None and s[0] == "id":
                cols = [c for c in s if c != "values"]
                continue
            if in_values and cols is not None and len(s) >= len(cols):
                try:
                    data.append([float(x) for x in s[:len(cols)]])
                except ValueError:
                    pass
    if cols is None or not data:
        return None
    arr = np.array(data)
    return {c: arr[:, k] for k, c in enumerate(cols)}


def grid(d, field):
    g = np.zeros((L, L))
    ix = np.clip(d["x"].astype(int), 0, L - 1)
    iy = np.clip(d["y"].astype(int), 0, L - 1)
    g[iy, ix] = d[field]
    return g


def main():
    resdir = sys.argv[1] if len(sys.argv) > 1 else "runs_spatial"

    figc, axc = plt.subplots(2, 3, figsize=(12, 8))
    figi, axi = plt.subplots(2, 3, figsize=(12, 8))

    for r, cond in enumerate(CONDS):
        for c, (t, lab) in enumerate(TIMES):
            # dump filename may be dump.2D_model.<t> (exact) — find it
            cand = glob.glob(os.path.join(resdir, cond, f"dump.2D_model.{t}"))
            if not cand:
                cand = glob.glob(os.path.join(resdir, cond, f"dump*{t}"))
            ax1, ax2 = axc[r, c], axi[r, c]
            if not cand:
                for ax in (ax1, ax2):
                    ax.text(.5, .5, f"no dump\n{cond} t={t}", ha="center", va="center")
                    ax.set_xticks([]); ax.set_yticks([])
                continue
            d = read_dump(cand[0])
            if d is None:
                continue
            cancer = grid(d, "i1")
            ifng = grid(d, "d1")
            arm = d["i2"] + d["i3"] + d["i4"]          # AR- CD8
            arp = d["i5"] + d["i6"] + d["i7"]          # AR+ CD8

            # cells panel: cancer density + T cell dots
            ax1.imshow(cancer, origin="lower", cmap="Reds", vmin=0, vmax=2)
            ix = np.clip(d["x"].astype(int), 0, L - 1)
            iy = np.clip(d["y"].astype(int), 0, L - 1)
            for vals, col in ((arm, "#1f6feb"), (arp, "#2ca02c")):
                m = vals > 0
                if m.any():
                    ax1.scatter(ix[m], iy[m], s=.4, c=col, linewidths=0)
            ax1.set_title(f"{cond} — {lab}"); ax1.set_xticks([]); ax1.set_yticks([])

            # ifng panel
            im = ax2.imshow(ifng, origin="lower", cmap="viridis")
            ax2.set_title(f"{cond} — {lab}"); ax2.set_xticks([]); ax2.set_yticks([])

    figc.suptitle("Spatial cell positions  (red = cancer density, blue = AR- CD8, green = AR+ CD8)")
    figi.suptitle("Free IFN-γ concentration")
    for fig, name in ((figc, "Fig_spatial_cells.png"), (figi, "Fig_spatial_ifng.png")):
        fig.tight_layout()
        out = os.path.join(resdir, name)
        fig.savefig(out, dpi=140)
        print("wrote", out)


if __name__ == "__main__":
    main()
