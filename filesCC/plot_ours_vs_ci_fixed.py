import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime

BENCH = "MJ_cni_benchmark_202607.xls.xls"
PANELS = [
    ("Historical", "ci_replication/hist",
     [("GPR_index.csv", "GPRH"), ("GPT_index.csv", "GPRHT"), ("GPA_index.csv", "GPRHA")]),
    ("Recent", "ci_replication/recent",
     [("GPR_index.csv", "GPR"), ("GPT_index.csv", "GPRT"), ("GPA_index.csv", "GPRA")]),
]

def parse_date(val):
    """Try multiple date formats"""
    if pd.isna(val):
        return None
    if isinstance(val, (int, float)):
        # Excel serial date
        try:
            return pd.Timestamp("1899-12-30") + pd.Timedelta(days=val)
        except:
            return None
    try:
        return pd.to_datetime(str(val))
    except:
        return None

def load_ours(path):
    """Load our CSV, detect date column, return series indexed by YM"""
    df = pd.read_csv(path)
    # Find date column (has lots of dates)
    date_col = None
    for col in df.columns:
        parsed = df[col].apply(parse_date)
        if parsed.notna().sum() > len(df) * 0.9:
            date_col = col
            break
    
    if date_col is None:
        raise ValueError(f"No date column found in {path}")
    
    df['date'] = df[date_col].apply(parse_date)
    df = df.dropna(subset=['date'])
    
    # Get index column (first numeric column that's not date)
    num_cols = [c for c in df.select_dtypes('number').columns if c != date_col]
    idx_col = num_cols[0]
    
    df['YM'] = df['date'].dt.to_period('M')
    result = df.set_index('YM')[idx_col]
    return result.loc[~result.index.duplicated(keep='first')]

# Load benchmark
bench = pd.read_excel(BENCH)
print(f"Benchmark shape: {bench.shape}")
print(f"Benchmark columns: {bench.columns.tolist()}")

# Identify date column in benchmark
date_col = None
for col in bench.columns:
    parsed = bench[col].apply(parse_date)
    if parsed.notna().sum() > len(bench) * 0.9:
        date_col = col
        break

if date_col is None:
    raise ValueError("Could not find date column in benchmark")

bench['date'] = bench[date_col].apply(parse_date)
bench['YM'] = bench['date'].dt.to_period('M')
bench = bench.set_index('YM')

# Plot
fig, axes = plt.subplots(2, 3, figsize=(18, 10))
fig.suptitle("Our replication vs Caldara & Iacoviello published indices", fontsize=14, fontweight='bold')

for r, (label, folder, pairs) in enumerate(PANELS):
    for c, (fname, bcol) in enumerate(pairs):
        ax = axes[r][c]
        try:
            ours = load_ours(f"{folder}/{fname}")
            theirs = bench[bcol].dropna()
            
            # Merge on index
            both = pd.DataFrame({'ours': ours, 'ci': theirs}).dropna()
            
            if len(both) == 0:
                ax.text(0.5, 0.5, 'No overlapping dates', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f"{label}: {fname[:3]} vs {bcol}")
                continue
            
            corr = both['ours'].corr(both['ci'])
            
            ax.plot(both.index.astype(str), both['ci'], 'ko-', lw=1.5, markersize=3, label=f'C&I {bcol}', alpha=0.7)
            ax.plot(both.index.astype(str), both['ours'], 'r-', lw=1, alpha=0.7, label=f'Ours {fname[:3]}')
            
            ax.set_title(f"{label}: {fname[:3]} vs {bcol}\n(r = {corr:.4f}, n = {len(both)})", fontsize=10)
            ax.legend(fontsize=8, loc='upper left')
            ax.grid(alpha=0.3)
            
            # Thin out x-axis labels for readability
            if len(both) > 50:
                step = len(both) // 10
                ax.set_xticks(range(0, len(both), step))
                ax.set_xticklabels([str(both.index[i]) for i in range(0, len(both), step)], rotation=45, fontsize=7)
            else:
                ax.set_xticklabels([str(i) for i in both.index], rotation=45, fontsize=7)
                
        except Exception as e:
            ax.text(0.5, 0.5, f'Error: {str(e)[:30]}', ha='center', va='center', transform=ax.transAxes, fontsize=8)
            ax.set_title(f"{label}: {fname[:3]} vs {bcol}")

fig.tight_layout()
fig.savefig("ours_vs_ci_fixed.png", dpi=150, bbox_inches='tight')
print("\n✓ Saved ours_vs_ci_fixed.png")
plt.show()
