#!/usr/bin/env python3
"""
parse_and_plot.py RESDIR

Reads the SPPARKS stats table from each RESDIR/<cond>_s<seed>/log.spparks,
pulls the cancer-cell count over time (the 'i1' column of diag_style array),
writes one timecourse CSV per run, and draws the male-vs-female tumour-growth
figure (absolute size + fold change), mean +/- SD over seeds.

SPPARKS prints a stats header line that starts with 'Time'. The diag_style
array columns are named i1,i2,... ; i1 = cancer (species A). If the header
uses different labels on your build, the script prints what it found so you
can point it at the right column with  --cancer-col NAME_or_INDEX.
"""
import sys, os, glob, argparse
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

UNITS_PER_DAY = 144.0        # 1 native unit = 10 min
DAY0 = 14.0                  # runs start at day 14


def parse_log(path, cancer_col):
    """Return (time[], cancer[]) from a SPPARKS log, or (None,None)."""
    header = None
    rows = []
    with open(path, errors="replace") as f:
        for line in f:
            t = line.split()
            if not t:
                continue
            if header is None:
                if t[0].lower() == "time":
                    header = t
                continue
            # data rows: all numeric, same width as header
            try:
                vals = [float(x) for x in t]
            except ValueError:
                if rows:        # table ended
                    break
                continue
            if len(vals) == len(header):
                rows.append(vals)
    if header is None or not rows:
        return None, None, header
    arr = np.array(rows)
    tcol = 0
    # choose cancer column
    if cancer_col is None:
        ccol = None
        for cand in ("sum(i1)", "i1"):     # SPPARKS labels the cancer column sum(i1)
            if cand in header:
                ccol = header.index(cand); break
        if ccol is None:
            ccol = 6            # Time Naccept Nreject Nsweeps CPU Propnsty <cancer>...
    elif cancer_col.isdigit():
        ccol = int(cancer_col)
    else:
        ccol = header.index(cancer_col)
    return arr[:, tcol], arr[:, ccol], header


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("resdir")
    ap.add_argument("--cancer-col", default=None,
                    help="header name or integer index of the cancer column")
    a = ap.parse_args()

    runs = {}   # cond -> list of (time, cancer)
    printed_header = False
    for log in sorted(glob.glob(os.path.join(a.resdir, "*_s*", "log.spparks"))):
        run = os.path.basename(os.path.dirname(log))          # cond_sSEED
        cond = run.rsplit("_s", 1)[0]
        t, c, header = parse_log(log, a.cancer_col)
        if t is None:
            print("  !! no stats table in %s" % log); continue
        if not printed_header:
            print("stats header detected:", header)
            print("using cancer column ->",
                  (a.cancer_col or ("i1" if "i1" in header else "index 1")))
            printed_header = True
        # per-run csv
        np.savetxt(log.replace("log.spparks", "timecourse.csv"),
                   np.column_stack([DAY0 + t / UNITS_PER_DAY, c]),
                   delimiter=",", header="day,cancer", comments="")
        runs.setdefault(cond, []).append((t, c))

    if not runs:
        sys.exit("no parseable logs found in %s" % a.resdir)

    # aggregate + plot
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    colors = {"male": "#1f4e8c", "female": "#b4245a",
              "male_noseq": "#6a9a3a", "female_real": "#d08a00"}
    summ = []
    for cond, series in sorted(runs.items()):
        # interpolate all seeds onto a common time grid
        tmax = min(s[0][-1] for s in series)
        grid = np.linspace(0, tmax, 200)
        C = np.array([np.interp(grid, t, c) for t, c in series])
        day = DAY0 + grid / UNITS_PER_DAY
        mean, sd = C.mean(0), C.std(0)
        fold = C / C[:, :1]
        fmean, fsd = fold.mean(0), fold.std(0)
        col = colors.get(cond, "#555")
        axes[0].plot(day, mean / 1e3, color=col, label=cond)
        axes[0].fill_between(day, (mean - sd) / 1e3, (mean + sd) / 1e3, color=col, alpha=.2)
        axes[1].plot(day, fmean, color=col, label=cond)
        axes[1].fill_between(day, fmean - fsd, fmean + fsd, color=col, alpha=.2)
        summ.append((cond, mean[0], mean[-1], fmean[-1], len(series)))

    axes[0].set(xlabel="Day", ylabel="Tumour size (thousands of cancer cells)",
                title="Absolute tumour size")
    axes[1].set(xlabel="Day", ylabel="Fold change", title="Cancer fold change")
    for ax in axes:
        ax.legend(frameon=False); ax.grid(alpha=.3)
    fig.tight_layout()
    out = os.path.join(a.resdir, "Fig_cancer_timecourse.png")
    fig.savefig(out, dpi=150)
    print("wrote", out)

    print("\nday-21 summary (cond: start -> end cancer, fold, nseed)")
    for cond, c0, c1, f1, n in summ:
        print("  %-12s %8.0f -> %8.0f   fold %.3f   (n=%d)" % (cond, c0, c1, f1, n))
    # male/female fold ratio (the thesis "1.29")
    d = {c: f for c, c0, c1, f, n in summ}
    if "male" in d and "female" in d:
        print("\n  male/female fold ratio = %.3f  (thesis ~1.29)"
              % (d["male"] / d["female"]))


if __name__ == "__main__":
    main()
