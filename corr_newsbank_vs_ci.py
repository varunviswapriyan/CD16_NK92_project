#!/usr/bin/env python3
"""ours (NewsBank / Liu & Zhang) vs C&I GPRC_* country indices."""
import argparse, glob, shutil, subprocess, sys, tempfile
import numpy as np, pandas as pd

ISO = {"DEU": "Germany", "JPN": "Japan", "GBR": "United Kingdom",
       "CHN": "Chinese mainland"}

def read_bench(path):
    def finish(df):
        cols = {str(c).strip().upper(): c for c in df.columns}
        if "GPR" not in cols or "GPRH" not in cols:
            return None
        keep = [k for k in cols if k.startswith("GPRC_") or k == "GPR"]
        out = pd.DataFrame({k: pd.to_numeric(df[cols[k]], errors="coerce")
                            for k in keep})
        out = out[df[cols["GPRH"]].notna().values].reset_index(drop=True)
        out.index = pd.period_range("1900-01", periods=len(out),
                                    freq="M").to_timestamp()
        return out
    try:
        x = pd.ExcelFile(path)
        for sh in x.sheet_names:
            g = finish(x.parse(sh))
            if g is not None:
                return g
    except ImportError:
        so = shutil.which("soffice") or shutil.which("libreoffice")
        if not so:
            sys.exit("install xlrd (pip install xlrd) or LibreOffice")
        td = tempfile.mkdtemp()
        subprocess.run([so, "--headless", "--convert-to", "csv",
                        "--outdir", td, path], capture_output=True, timeout=600)
        for f in glob.glob(f"{td}/*.csv"):
            g = finish(pd.read_csv(f, low_memory=False))
            if g is not None:
                return g
    sys.exit(f"could not read {path}")

def read_ours(path):
    d = pd.read_csv(path).rename(columns=lambda c: str(c).strip())
    d = d.rename(columns={d.columns[0]: "month"})
    d["month"] = pd.to_datetime(d["month"], errors="coerce")
    d = d.dropna(subset=["month"]).set_index("month").sort_index()
    d.index = d.index.to_period("M").to_timestamp()
    return d.apply(pd.to_numeric, errors="coerce")

def nw_p(x, y, lags=12):
    try:
        import statsmodels.api as sm
    except ImportError:
        return np.nan
    x = (x - x.mean()) / x.std(); y = (y - y.mean()) / y.std()
    m = sm.OLS(y.values, sm.add_constant(x.values)).fit(
        cov_type="HAC", cov_kwds={"maxlags": lags})
    return float(m.pvalues[1])

def scores(o, t):
    j = pd.concat([o.rename("o"), t.rename("t")], axis=1).dropna()
    if len(j) < 24:
        return None
    d = j.diff(12).dropna()
    return {"n": len(j), "r_level": j.o.corr(j.t),
            "r_log": np.log1p(j.o.clip(lower=0)).corr(np.log1p(j.t.clip(lower=0))),
            "r_d12": d.o.corr(d.t),
            "rho": j.o.corr(j.t, method="spearman"),
            "p_nw": nw_p(j.o, j.t)}

ap = argparse.ArgumentParser()
ap.add_argument("--ours", required=True)
ap.add_argument("--bench", default="MJ_cni_benchmark_202607.xls.xls")
ap.add_argument("--from", dest="y0", type=int, default=1985)
ap.add_argument("--to", dest="y1", type=int, default=2019)
ap.add_argument("--out", default="newsbank_vs_ci.csv")
a = ap.parse_args()

B, O = read_bench(a.bench), read_ours(a.ours)
win = lambda s: s[(s.index.year >= a.y0) & (s.index.year <= a.y1)].dropna()

rows, used = [], []
for code in O.columns:
    c = str(code).strip().upper()
    col = f"GPRC_{c}"
    if col not in B.columns:
        print(f"  {c}: no {col} in the workbook, skipped"); continue
    s = scores(win(O[code]), win(B[col]))
    if s:
        rows.append({"country": ISO.get(c, c), "iso3": c, **s}); used.append(code)

t = pd.DataFrame(rows)
print(f"\nPer country, {a.y0}-{a.y1}, ours vs C&I GPRC_*\n")
print(t.to_string(index=False, float_format=lambda x: f"{x:.3f}"))
t.to_csv(a.out, index=False)

z = lambda s: (s - s.mean()) / s.std()
if used and "GPR" in B.columns:
    avg = pd.concat([z(win(O[c])) for c in used], axis=1).mean(axis=1)
    s = scores(avg, win(B["GPR"]))
    print("\nEqual-weighted average of your countries vs C&I global GPR:")
    print(f"  r_level {s['r_level']:.3f}   r_d12 {s['r_d12']:.3f}   n = {s['n']}")
    print("  Liu & Zhang report 0.58 equal-weighted / 0.78 GDP-weighted "
          "across all 53.")
print(f"\nwrote {a.out}")
