import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

BENCH = "MJ_cni_benchmark_202607.xls.xls"
PANELS = [
    ("Historical", "ci_replication/hist",
     [("GPR_index.csv", "GPRH"), ("GPT_index.csv", "GPRHT"), ("GPA_index.csv", "GPRHA")]),
    ("Recent", "ci_replication/recent",
     [("GPR_index.csv", "GPR"), ("GPT_index.csv", "GPRT"), ("GPA_index.csv", "GPRA")]),
]

def to_month(s):
    return pd.to_datetime(s, errors="coerce").dt.to_period("M").dt.to_timestamp()

def load_ours(path):
    df = pd.read_csv(path)
    date_col = next(c for c in df.columns if to_month(df[c]).notna().mean() > 0.9)
    num = [c for c in df.select_dtypes("number").columns if c != date_col]
    val = next((c for c in num if "index" in c.lower() or c.upper().startswith("GP")), num[-1])
    result = pd.Series(df[val].values, index=to_month(df[date_col]))
    return result.loc[~result.index.duplicated(keep='first')].dropna()

bench = pd.read_excel(BENCH)
bdate = next(c for c in bench.columns if to_month(bench[c]).notna().mean() > 0.9)
bench.index = to_month(bench[bdate])
bench = bench.loc[~bench.index.duplicated(keep='first')]

fig, axes = plt.subplots(2, 3, figsize=(18, 8))
for r, (label, folder, pairs) in enumerate(PANELS):
    for c, (fname, bcol) in enumerate(pairs):
        ax = axes[r][c]
        ours = load_ours(f"{folder}/{fname}")
        theirs = bench[bcol].dropna()
        # Align to overlapping dates
        both = pd.DataFrame({"ours": ours, "ci": theirs}).dropna()
        if len(both) > 0:
            corr = both["ours"].corr(both["ci"])
            ax.plot(both.index, both["ci"], color="black", lw=1.2, label=f"C&I {bcol}")
            ax.plot(both.index, both["ours"], color="tab:red", lw=0.9, alpha=0.8, label=f"Ours {fname[:3]}")
            ax.set_title(f"{label}: {fname[:3]} vs {bcol}  (r = {corr:.3f}, n = {len(both)})", fontsize=10)
            ax.legend(fontsize=8, loc="upper left")
            ax.grid(alpha=0.3)

fig.suptitle("Our replication vs Caldara & Iacoviello published indices", fontsize=13)
fig.tight_layout()
fig.savefig("ours_vs_ci.png", dpi=150)
print("saved ~/Downloads/ours_vs_ci.png")
plt.show()
