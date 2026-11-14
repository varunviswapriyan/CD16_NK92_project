#!/usr/bin/env python3
"""
convert_zihai.py -- turn one Zihai mIF slide (Processed_Zihai_Data/*.csv: one row per cell with
X/Y position in um and CD3/CD8/PD1/Tim3/TCF1/AR calls) into a SPPARKS starting file in the
17-column CD4 format our model reads:
    id i1 i2 i3 i4 i5 i6 i7 d1 i8 i9 i10 i11 i12 d2 i13 i14

Rules (reverse-engineered from Giuseppe's d14_M_4 / d14_F_4 files and checked against them):
  lattice      600 x 600 sites of 10 um (6 x 6 mm window of the slide); site = floor(um / 10)
  i1  cancer   a random 70% of CD3- cells (thesis); at most 2 per site, at most 1 where a T cell sits
  i2  AR- CD8 progenitor   CD3+ CD8+ AR- TCF1+ Tim3-
  i3  AR- CD8 mature       CD3+ CD8+ AR- TCF1- Tim3-
  i4  AR- CD8 exhausted    CD3+ CD8+ AR- Tim3+
  i5  AR+ progenitor       always 0 (new model has no AR+ progenitor)
  i6  AR+ CD8 mature       CD3+ CD8+ AR+ Tim3-
  i7  AR+ CD8 exhausted    CD3+ CD8+ AR+ Tim3+
  i13 CD4 AR-              CD3+ CD8- AR-
  i14 CD4 AR+              CD3+ CD8- AR+
  i8  recruitment trigger  1 at site id 1 only
  i11 lattice edge         1 on the outer ring of sites
  i12 recruitment allowed  1 in every 500 x 500 um block that contains any cell, else 0
  d1, d2, i9, i10          0
Window: if a reference starting file is given (--ref) the window is found by matching CD8 positions
to it (so the same piece of tissue Giuseppe used); otherwise the 6 x 6 mm window holding most cells.

Usage:
  python3 convert_zihai.py SLIDE.csv OUT_init [--ref GIUSEPPE_init] [--seed 1]
  python3 convert_zihai.py --check SLIDE.csv REF_init      (rebuild + compare, writes nothing)
"""
import sys, csv, argparse
import numpy as np

L, UM, BLOCK = 600, 10.0, 50
COLS = "id i1 i2 i3 i4 i5 i6 i7 d1 i8 i9 i10 i11 i12 d2 i13 i14".split()


def read_slide(path):
    with open(path, newline="") as fh:
        r = csv.reader(fh); h = next(r); ix = {c: i for i, c in enumerate(h)}
        X, Y, cd3, cd8, ar, tcf, tim, tis = [], [], [], [], [], [], [], []
        for row in r:
            X.append(float(row[ix["X Position"]])); Y.append(float(row[ix["Y Position"]]))
            cd3.append(row[ix["CD3"]] == "CD3+"); cd8.append(row[ix["CD8"]] == "CD8+")
            ar.append(row[ix["AR"]] == "AR+"); tcf.append(row[ix["TCF1"]] == "TCF1+")
            tim.append(row[ix["Tim3"]] == "Tim3+"); tis.append(row[ix["Tissue Category"]])
    a = lambda v, t=bool: np.array(v, dtype=t)
    return dict(x=a(X, float), y=a(Y, float), cd3=a(cd3), cd8=a(cd8), ar=a(ar), tcf=a(tcf), tim=a(tim),
                tumor=np.array([t == "Tumor" for t in tis]))


def types(s):
    c8 = s["cd3"] & s["cd8"]; c4 = s["cd3"] & ~s["cd8"]; ar, tcf, tim = s["ar"], s["tcf"], s["tim"]
    return {"i2": c8 & ~ar & tcf & ~tim, "i3": c8 & ~ar & ~tcf & ~tim, "i4": c8 & ~ar & tim,
            "i6": c8 & ar & ~tim, "i7": c8 & ar & tim, "i13": c4 & ~ar, "i14": c4 & ar}


def read_init(path):
    cols, rows = None, []
    for line in open(path):
        p = line.split()
        if p and p[0] == "id": cols = [c for c in p if c != "values"]; continue
        if cols and len(p) == len(cols):
            try: rows.append([float(v) for v in p])
            except ValueError: pass
    a = np.array(rows); n = int(round(np.sqrt(len(a))))
    return {c: a[:, i].reshape(n, n) for i, c in enumerate(cols)}


def bins(s, fx=0.0, fy=0.0):
    return (np.floor((s["x"] - s["x"].min() + fx) / UM).astype(int),
            np.floor((s["y"] - s["y"].min() + fy) / UM).astype(int))


def grid(X, Y, m, c0, r0):
    k = m & (X >= c0) & (X < c0 + L) & (Y >= r0) & (Y < r0 + L)
    g = np.zeros((L, L)); np.add.at(g, (Y[k] - r0, X[k] - c0), 1); return g


def find_window(s, ref=None):
    """(col0, row0, fx, fy). With ref: match CD8 positions; else window with most cells."""
    X, Y = bins(s)
    if ref is None:
        H, W = Y.max() + 1, X.max() + 1
        A = np.zeros((max(H, L), max(W, L))); np.add.at(A, (Y, X), 1)
        C = np.pad(A.cumsum(0).cumsum(1), ((1, 0), (1, 0)))
        S = C[L:, L:] - C[:-L, L:] - C[L:, :-L] + C[:-L, :-L]
        r, c = np.unravel_index(S.argmax(), S.shape); return int(c), int(r), 0.0, 0.0
    T = types(s); c8 = s["cd3"] & s["cd8"]
    R8 = sum(ref[k] for k in ("i2", "i3", "i4", "i5", "i6", "i7") if k in ref)
    N = max(Y.max(), X.max()) + 1 + 2 * L
    P = np.zeros((N, N)); np.add.at(P, (Y[c8] + L, X[c8] + L), 1)          # pad so negative offsets work
    K = np.zeros((N, N)); K[:L, :L] = R8
    cc = np.real(np.fft.ifft2(np.fft.fft2(P) * np.conj(np.fft.fft2(K))))
    r, c = np.unravel_index(cc.argmax(), cc.shape); r0, c0 = r - L, c - L
    best = (-1, c0, r0, 0.0, 0.0)
    for fx in range(10):                                                    # sub-site shift, 1 um steps
        for fy in range(10):
            X2, Y2 = bins(s, fx, fy)
            for dc in (-1, 0, 1):
                for dr in (-1, 0, 1):
                    g = grid(X2, Y2, c8, c0 + dc, r0 + dr)
                    sc = ((g == R8) & (R8 > 0)).sum()
                    if sc > best[0]: best = (sc, c0 + dc, r0 + dr, float(fx), float(fy))
    return best[1], best[2], best[3], best[4]


def build(s, c0, r0, fx=0.0, fy=0.0, seed=1):
    X, Y = bins(s, fx, fy); T = types(s); out = {c: np.zeros((L, L)) for c in COLS[1:]}
    for k, m in T.items(): out[k] = grid(X, Y, m, c0, r0)
    rng = np.random.default_rng(seed)
    neg = ~s["cd3"]; pick = neg & (rng.random(len(neg)) < 0.7)
    can = grid(X, Y, pick, c0, r0)
    tcell = sum(out[k] for k in ("i2", "i3", "i4", "i6", "i7", "i13", "i14"))
    out["i1"] = np.where(tcell > 0, np.minimum(can, 1), np.minimum(can, 2))
    allc = grid(X, Y, np.ones(len(neg), bool), c0, r0)
    blk = allc.reshape(L // BLOCK, BLOCK, L // BLOCK, BLOCK).sum((1, 3)) > 0
    out["i12"] = np.kron(blk, np.ones((BLOCK, BLOCK)))
    out["i11"][0, :] = out["i11"][-1, :] = out["i11"][:, 0] = out["i11"][:, -1] = 1
    out["i8"][0, 0] = 1
    return out


def write(out, path):
    with open(path, "w") as fh:
        fh.write("Site file written by dump sites 1 command\n\n" + " ".join(COLS) + " values\n\nValues\n\n")
        M = np.column_stack([np.arange(1, L * L + 1)] + [out[c].ravel().astype(int) for c in COLS[1:]])
        np.savetxt(fh, M, fmt="%d")


def compare(out, ref, label=""):
    """Per-column totals and how many CD8/cancer sites agree exactly. Returns dict of checks."""
    print(f"   {'column':6s} {'ours':>9s} {'Giuseppe':>9s}  diff")
    res = {}
    for c in ("i1", "i2", "i3", "i4", "i5", "i6", "i7", "i11", "i12", "i13", "i14"):
        if c not in ref: continue
        a, b = out[c].sum(), ref[c].sum()
        print(f"   {c:6s} {a:9.0f} {b:9.0f}  {100 * (a - b) / b if b else 0:+.1f}%")
    o8 = sum(out[k] for k in ("i2", "i3", "i4", "i6", "i7"))
    r8 = sum(ref[k] for k in ("i2", "i3", "i4", "i5", "i6", "i7") if k in ref)
    res["cd8_total"] = (o8.sum() - r8.sum()) / r8.sum()
    res["cd8_site_match"] = ((o8 == r8) & (r8 > 0)).sum() / max((r8 > 0).sum(), 1)
    res["cancer_total"] = (out["i1"].sum() - ref["i1"].sum()) / ref["i1"].sum()
    res["cancer_site_match"] = (((out["i1"] > 0) == (ref["i1"] > 0))).mean()
    if "i13" in ref:
        o4, r4 = out["i13"] + out["i14"], ref["i13"] + ref["i14"]
        res["cd4_total"] = (o4.sum() - r4.sum()) / r4.sum()
    arp_o = out["i6"].sum() + out["i7"].sum(); arp_r = sum(ref[k].sum() for k in ("i5", "i6", "i7") if k in ref)
    res["arp_total"] = (arp_o - arp_r) / arp_r if arp_r else 0
    print(f"   CD8 total {100 * res['cd8_total']:+.1f}%, CD8 sites identical {100 * res['cd8_site_match']:.0f}%, "
          f"cancer total {100 * res['cancer_total']:+.1f}%, cancer occupancy identical {100 * res['cancer_site_match']:.0f}%"
          + (f", CD4 total {100 * res['cd4_total']:+.1f}%" if "cd4_total" in res else "")
          + f", AR+ CD8 total {100 * res['arp_total']:+.1f}%")
    return res


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("slide"); ap.add_argument("out"); ap.add_argument("--ref"); ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--check", action="store_true", help="out is a reference file: rebuild and compare only")
    a = ap.parse_args()
    s = read_slide(a.slide)
    ref = read_init(a.out if a.check else a.ref) if (a.check or a.ref) else None
    c0, r0, fx, fy = find_window(s, ref)
    print(f"{a.slide}: {len(s['x'])} cells, window starts at ({c0 * 10:.0f}, {r0 * 10:.0f}) um, shift ({fx:.0f}, {fy:.0f}) um")
    out = build(s, c0, r0, fx, fy, a.seed)
    if ref is not None: compare(out, ref)
    if not a.check: write(out, a.out); print("   wrote", a.out)
