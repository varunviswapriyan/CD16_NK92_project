"""
fit_ca_variants_c10_c20.py
==========================
Expanded Ca-stage sweep, continuing the c0x series. 11 variants, one command.

PER-CONDITION REPLICATION variants (same ODE, one parameter replicated per
condition -- these are the ones that can actually close the FceRIg gap):

  c10  per-condition C1   (input-coupling gain)
  c11  per-condition g    (downstream gain/decay)
  c12  per-condition Ca0  (baseline offset control)
  c19  per-condition C2   (input sensitivity / EC50)
  c20  per-condition k3   (inactivation drive)

SHARED-PARAMETER variants (all 3 conditions share one param set -- objective
and identifiability tweaks; controls + best possible shared-param compromise):

  c13  k4 fixed (default 0.0618 from c08), raw SSR
  c14  condition-balanced SSR (per-condition variance normalization)
  c15  peak-window weighting (3x weight, t = 40-130 s)
  c17  log-residual SSR (relative-error fit, tames amplitude dominance)
  c18  combo: k4 fixed + balanced + peak-window (best shared-param shot)

AUTO-MAIN:

  c16  winner among {c10, c11, c19, c20} + k4 fixed + balanced SSR
       + big multistart. Built automatically after the others finish.

Every variant uses the SAME ODE / model rules. Each is wrapped in try/except
so one failure never kills the sweep: you always end with a summary table
(raw SSR reported for ALL variants -- apples-to-apples with c02 = 8.402e4),
plots plot_ca_c10.png ... plot_ca_c20.png in the c0x style, and results.json.

-------------------------------------------------------------------------
!! TWO THINGS TO WIRE UP BEFORE RUNNING (search for "ADAPTER") !!
  1) DATA PATHS in CONFIG (Ca traces + pZAP input curves per condition)
  2) ca_rhs(): paste your EXACT Ca ODE right-hand side from the c02 code.
     A plausible placeholder is provided so the script runs end-to-end,
     but the placeholder is NOT your model -- replace it.
-------------------------------------------------------------------------

Usage:
    python fit_ca_variants_c10_c20.py            # run everything
    python fit_ca_variants_c10_c20.py c10 c13    # run a subset
"""

import json
import sys
import time as _time
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import interp1d
from scipy.optimize import minimize

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ========================================================================
# CONFIG  -- ADAPTER (1/2): point these at your files
# ========================================================================
CONDITIONS = ["CD3z", "Hetero", "FceRIg"]
COLORS = {"CD3z": "tab:blue", "Hetero": "tab:orange", "FceRIg": "tab:green"}

CONFIG = {
    # Ca experimental data: CSV per condition with columns: time, ca  (+ optional se)
    "ca_data_files": {
        "CD3z":   "data/ca_CD3z.csv",
        "Hetero": "data/ca_Hetero.csv",
        "FceRIg": "data/ca_FceRIg.csv",
    },
    # pZAP input curves from the WINNING pZAP stage (c02 / v58):
    # CSV per condition with columns: time, pzap
    "pzap_files": {
        "CD3z":   "data/pzap_v58_CD3z.csv",
        "Hetero": "data/pzap_v58_Hetero.csv",
        "FceRIg": "data/pzap_v58_FceRIg.csv",
    },
    "out_dir": "out_ca_variants",

    # Shared-parameter bounds (wide, matching your c0x philosophy).
    # Order everywhere: C1, C2, g, k3, k4
    "bounds": {
        "C1": (1e-2, 1e2),
        "C2": (1e-4, 1e1),
        "g":  (1e-5, 1e0),
        "k3": (1e-4, 1e1),
        "k4": (1e-3, 1e2),
    },
    "log_scale_params": ["C1", "C2", "g", "k3", "k4"],  # sample/search in log10

    "k4_fixed_value": 0.0618,     # from c08; change if you prefer another pin
    "peak_window": (40.0, 130.0), # seconds, for c15 / c18
    "peak_weight": 3.0,

    # PSO settings
    "pso": {"n_particles": 40, "n_iter": 120, "w": 0.72, "c1": 1.5, "c2": 1.5},
    "n_multistart_default": 4,
    "n_multistart_main": 12,      # for c16
    "seed": 58,
}

# ========================================================================
# DATA LOADING
# ========================================================================
def _load_xy_csv(path, ycol_candidates):
    arr = np.genfromtxt(path, delimiter=",", names=True)
    names = arr.dtype.names
    tcol = names[0]
    ycol = next((c for c in ycol_candidates if c in names), names[1])
    out = {"t": np.asarray(arr[tcol], float), "y": np.asarray(arr[ycol], float)}
    if "se" in names:
        out["se"] = np.asarray(arr["se"], float)
    return out


def load_all_data(cfg):
    ca, pzap = {}, {}
    for c in CONDITIONS:
        d = _load_xy_csv(cfg["ca_data_files"][c], ("ca", "signal", "y"))
        ca[c] = d
        p = _load_xy_csv(cfg["pzap_files"][c], ("pzap", "y"))
        pzap[c] = interp1d(p["t"], p["y"], kind="linear",
                           bounds_error=False,
                           fill_value=(p["y"][0], p["y"][-1]))
    return ca, pzap

# ========================================================================
# MODEL  -- ADAPTER (2/2): paste your real Ca ODE here
# ========================================================================
PARAM_NAMES = ["C1", "C2", "g", "k3", "k4"]

def ca_rhs(t, y, p, pzap_func):
    """
    !! PLACEHOLDER ODE -- REPLACE WITH YOUR EXACT c02 Ca MODEL !!
    State y = [Ca, h].  Params p = dict with C1, C2, g, k3, k4 (+ Ca0 for c12).
    The placeholder below is only there so the pipeline runs end-to-end.
    """
    Ca, h = y
    z = pzap_func(t)
    act = p["C1"] * z / (p["C2"] + z) * h
    dCa = act - p["g"] * Ca
    dh = -p["k3"] * z * h + p["k4"] * (1.0 - h)
    return [dCa, dh]


def initial_state(p, ca_data_first):
    """Initial condition. c12 overrides Ca0 per condition; default = first data pt."""
    Ca0 = p.get("Ca0", ca_data_first)
    return [Ca0, 1.0]


def simulate(p, pzap_func, t_eval, ca_first):
    sol = solve_ivp(ca_rhs, (t_eval[0], t_eval[-1]),
                    initial_state(p, ca_first),
                    t_eval=t_eval, args=(p, pzap_func),
                    method="LSODA", rtol=1e-6, atol=1e-8)
    if not sol.success:
        return None
    return sol.y[0]

# ========================================================================
# VARIANT SPECS
# ========================================================================
# per_condition: params replicated per condition; fixed: {name: value};
# objective: "raw" | "balanced" | "peakweight" | "logres" | "balanced+peak";
# extra_per_condition: e.g. Ca0
VARIANTS = {
    # ---- per-condition replication family --------------------------------
    "c10": dict(label="per-condition C1 (coupling gain)",
                per_condition=["C1"], fixed={}, objective="raw"),
    "c11": dict(label="per-condition g (downstream gain)",
                per_condition=["g"], fixed={}, objective="raw"),
    "c12": dict(label="per-condition Ca0 (offset control)",
                per_condition=[], fixed={}, objective="raw",
                extra_per_condition=["Ca0"]),
    "c19": dict(label="per-condition C2 (input sensitivity)",
                per_condition=["C2"], fixed={}, objective="raw"),
    "c20": dict(label="per-condition k3 (inactivation drive)",
                per_condition=["k3"], fixed={}, objective="raw"),
    # ---- shared-parameter family -----------------------------------------
    "c13": dict(label=f"shared, k4 fixed @ {CONFIG['k4_fixed_value']}",
                per_condition=[], fixed={"k4": CONFIG["k4_fixed_value"]},
                objective="raw"),
    "c14": dict(label="shared, condition-balanced SSR",
                per_condition=[], fixed={}, objective="balanced"),
    "c15": dict(label="shared, peak-window weighted (3x, 40-130s)",
                per_condition=[], fixed={}, objective="peakweight"),
    "c17": dict(label="shared, log-residual SSR (relative error)",
                per_condition=[], fixed={}, objective="logres"),
    "c18": dict(label="shared combo: k4 fixed + balanced + peak-window",
                per_condition=[], fixed={"k4": CONFIG["k4_fixed_value"]},
                objective="balanced+peak"),
    # c16 (auto-MAIN) is built dynamically after the per-condition family runs
}

PER_COND_FAMILY = ["c10", "c11", "c12", "c19", "c20"]
SHARED_FAMILY = ["c13", "c14", "c15", "c17", "c18"]
# c16 winner is chosen among true parameter-replication variants (not the
# Ca0 offset control, which is expected to be insufficient alone):
TIER1_POOL = ["c10", "c11", "c19", "c20"]
TIER1_PARAM = {"c10": "C1", "c11": "g", "c19": "C2", "c20": "k3"}

# ========================================================================
# PARAMETER VECTOR MACHINERY
# ========================================================================
class ParamSpec:
    """Flattens shared + per-condition params into one optimization vector."""

    def __init__(self, variant, cfg, ca_data):
        self.entries = []  # (vector_name, base_name, condition_or_None, lo, hi, log?)
        b = cfg["bounds"]
        logset = set(cfg["log_scale_params"])
        percond = set(variant.get("per_condition", []))
        fixed = variant.get("fixed", {})
        self.fixed = dict(fixed)

        for name in PARAM_NAMES:
            if name in fixed:
                continue
            lo, hi = b[name]
            if name in percond:
                for c in CONDITIONS:
                    self.entries.append((f"{name}_{c}", name, c, lo, hi, name in logset))
            else:
                self.entries.append((name, name, None, lo, hi, name in logset))

        for name in variant.get("extra_per_condition", []):
            if name == "Ca0":
                for c in CONDITIONS:
                    y0 = ca_data[c]["y"][0]
                    self.entries.append((f"Ca0_{c}", "Ca0", c, 0.3 * y0, 3.0 * y0, False))

        self.n = len(self.entries)

    def bounds_vec(self):
        lo, hi = [], []
        for (_, _, _, l, h, is_log) in self.entries:
            lo.append(np.log10(l) if is_log else l)
            hi.append(np.log10(h) if is_log else h)
        return np.array(lo), np.array(hi)

    def to_params(self, x):
        """vector -> per-condition param dicts"""
        per = {c: dict(self.fixed) for c in CONDITIONS}
        for xi, (_, base, cond, _, _, is_log) in zip(x, self.entries):
            val = 10.0 ** xi if is_log else xi
            if cond is None:
                for c in CONDITIONS:
                    per[c][base] = val
            else:
                per[cond][base] = val
        return per

    def describe(self, x):
        out = {}
        for xi, (vname, _, _, _, _, is_log) in zip(x, self.entries):
            out[vname] = float(10.0 ** xi if is_log else xi)
        out.update({k: float(v) for k, v in self.fixed.items()})
        return out

# ========================================================================
# OBJECTIVES
# ========================================================================
def make_objective(spec, variant, ca_data, pzap, cfg):
    mode = variant["objective"]
    t0, t1 = cfg["peak_window"]
    pw = cfg["peak_weight"]
    var_per_cond = {c: np.var(ca_data[c]["y"]) for c in CONDITIONS}

    def sims_per_cond(x):
        per = spec.to_params(x)
        out = {}
        for c in CONDITIONS:
            t = ca_data[c]["t"]; y = ca_data[c]["y"]
            sim = simulate(per[c], pzap[c], t, y[0])
            if sim is None or not np.all(np.isfinite(sim)):
                return None
            out[c] = (sim, y, t)
        return out

    def objective(x):
        sims = sims_per_cond(x)
        if sims is None:
            return 1e12
        total = 0.0
        for c, (sim, y, t) in sims.items():
            if mode == "logres":
                # relative-error fit; data floor keeps logs safe
                eps = 1e-9
                res2 = (np.log(np.clip(sim, eps, None))
                        - np.log(np.clip(y, eps, None))) ** 2
                total += res2.sum()
                continue
            res2 = (sim - y) ** 2
            if mode == "balanced":
                total += res2.sum() / var_per_cond[c]
            elif mode == "peakweight":
                w = np.where((t >= t0) & (t <= t1), pw, 1.0)
                total += (w * res2).sum()
            elif mode == "balanced+peak":
                w = np.where((t >= t0) & (t <= t1), pw, 1.0)
                total += (w * res2).sum() / var_per_cond[c]
            else:  # raw
                total += res2.sum()
        return total

    def raw_total(x):
        sims = sims_per_cond(x)
        if sims is None:
            return np.inf
        return float(sum(((sim - y) ** 2).sum() for sim, y, _ in sims.values()))

    return objective, raw_total

# ========================================================================
# PSO + POLISH
# ========================================================================
def pso(objective, lo, hi, cfg, rng):
    ps = cfg["pso"]
    n, d = ps["n_particles"], len(lo)
    X = rng.uniform(lo, hi, (n, d))
    V = np.zeros_like(X)
    F = np.array([objective(x) for x in X])
    Pb, Pf = X.copy(), F.copy()
    gi = int(np.argmin(F)); gx, gf = X[gi].copy(), F[gi]
    for _ in range(ps["n_iter"]):
        r1 = rng.uniform(size=(n, d)); r2 = rng.uniform(size=(n, d))
        V = ps["w"] * V + ps["c1"] * r1 * (Pb - X) + ps["c2"] * r2 * (gx - X)
        X = np.clip(X + V, lo, hi)
        F = np.array([objective(x) for x in X])
        imp = F < Pf
        Pb[imp], Pf[imp] = X[imp], F[imp]
        gi = int(np.argmin(Pf))
        if Pf[gi] < gf:
            gx, gf = Pb[gi].copy(), Pf[gi]
    return gx, gf


def fit_variant(spec, variant, ca_data, pzap, cfg, n_starts, rng):
    objective, raw_total = make_objective(spec, variant, ca_data, pzap, cfg)
    lo, hi = spec.bounds_vec()
    best_x, best_f = None, np.inf
    for s in range(n_starts):
        x, f = pso(objective, lo, hi, cfg, rng)
        res = minimize(objective, x, method="Nelder-Mead",
                       options={"maxiter": 4000, "xatol": 1e-8, "fatol": 1e-8})
        xf, ff = (res.x, res.fun) if res.fun < f else (x, f)
        xf = np.clip(xf, lo, hi)
        ff = objective(xf)
        if ff < best_f:
            best_x, best_f = xf, ff
        print(f"    start {s+1}/{n_starts}: obj={ff:.4g}  (best so far {best_f:.4g})")
    return best_x, best_f, raw_total(best_x)

# ========================================================================
# PLOTTING (c0x style: pZAP inputs left, Ca fit right)
# ========================================================================
def plot_variant(cid, variant, spec, x, raw_ssr, ca_data, pzap, cfg, outdir):
    per = spec.to_params(x)
    desc = spec.describe(x)
    fig, (axL, axR) = plt.subplots(1, 2, figsize=(19.2, 7.2))

    tt = np.linspace(0, 300, 600)
    for c in CONDITIONS:
        axL.plot(tt, pzap[c](tt), color=COLORS[c], label=c)
    axL.set_xlabel("Time (s)"); axL.set_ylabel("pZAP input (uM)")
    axL.set_title("pZAP input to Ca ODE"); axL.legend()

    for c in CONDITIONS:
        t, y = ca_data[c]["t"], ca_data[c]["y"]
        axR.plot(t, y, ".", ms=3, color=COLORS[c], alpha=0.6, label=f"{c} exp")
        sim = simulate(per[c], pzap[c], t, y[0])
        if sim is not None:
            axR.plot(t, sim, "-", lw=2.5, color=COLORS[c], label=f"{c} model")
    axR.set_xlabel("Time (s)"); axR.set_ylabel("Ca signal")
    axR.set_title(f"Ca fit   total SSR={raw_ssr:.3e}"); axR.legend(ncol=2, fontsize=8)

    ptxt = "  ".join(f"{k}={v:.3g}" for k, v in desc.items())
    fig.suptitle(f"{cid}: {variant['label']}   {ptxt}", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    path = Path(outdir) / f"plot_ca_{cid}.png"
    fig.savefig(path, dpi=100); plt.close(fig)
    return str(path)

# ========================================================================
# MAIN SWEEP
# ========================================================================
def main():
    cfg = CONFIG
    outdir = Path(cfg["out_dir"]); outdir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(cfg["seed"])

    requested = [a for a in sys.argv[1:] if a.startswith("c")] or None

    print("Loading data...")
    ca_data, pzap = load_all_data(cfg)

    results = {}

    def run_one(cid, variant, n_starts):
        print(f"\n=== {cid}: {variant['label']} ===")
        t0 = _time.time()
        try:
            spec = ParamSpec(variant, cfg, ca_data)
            x, obj, raw_ssr = fit_variant(spec, variant, ca_data, pzap, cfg,
                                          n_starts, rng)
            plot = plot_variant(cid, variant, spec, x, raw_ssr,
                                ca_data, pzap, cfg, outdir)
            results[cid] = dict(label=variant["label"], ok=True,
                                objective_value=float(obj),
                                raw_ssr=float(raw_ssr),
                                params=spec.describe(x),
                                plot=plot,
                                seconds=round(_time.time() - t0, 1))
            print(f"  DONE  raw SSR = {raw_ssr:.4e}  ({results[cid]['seconds']} s)")
        except Exception as e:  # keep the sweep alive no matter what
            results[cid] = dict(label=variant["label"], ok=False, error=repr(e))
            print(f"  FAILED: {e!r}  -- continuing with remaining variants")

    # Cheap/informative shared runs first, then the replication family
    order = ["c13", "c14", "c15", "c17", "c18",
             "c10", "c11", "c19", "c20", "c12"]
    for cid in order:
        if requested and cid not in requested:
            continue
        run_one(cid, VARIANTS[cid], cfg["n_multistart_default"])

    # ---- c16: auto-MAIN, combine winners ---------------------------------
    if (requested is None) or ("c16" in requested):
        tier1 = {k: results[k]["raw_ssr"] for k in TIER1_POOL
                 if k in results and results[k].get("ok")}
        winner_param = "C1"  # sensible default if tier1 missing
        if tier1:
            winner = min(tier1, key=tier1.get)
            winner_param = TIER1_PARAM[winner]
            print(f"\nTier-1 winner: {winner} (per-condition {winner_param})")
        c16 = dict(label=f"MAIN: per-cond {winner_param} + k4 fixed + balanced SSR",
                   per_condition=[winner_param],
                   fixed={"k4": cfg["k4_fixed_value"]},
                   objective="balanced")
        run_one("c16", c16, cfg["n_multistart_main"])

    # ---- summary ----------------------------------------------------------
    print("\n" + "=" * 76)
    print(f"{'id':<5} {'family':<10} {'raw SSR':>12}   status / label")
    print("-" * 76)
    baseline = 8.402e4  # c02 for reference

    def fam(cid):
        if cid in PER_COND_FAMILY:
            return "per-cond"
        if cid in SHARED_FAMILY:
            return "shared"
        return "MAIN"

    for cid in sorted(results):
        r = results[cid]
        if r.get("ok"):
            print(f"{cid:<5} {fam(cid):<10} {r['raw_ssr']:>12.4e}   OK  {r['label']}")
        else:
            print(f"{cid:<5} {fam(cid):<10} {'--':>12}   FAILED  {r['label']}: {r['error']}")
    print("-" * 76)
    print(f"{'c02':<5} {'(ref)':<10} {baseline:>12.4e}   previous MAIN, for reference")

    with open(outdir / "results.json", "w") as f:
        json.dump(results, f, indent=2)
    print(f"\nAll outputs in: {outdir}/  (plots + results.json)")


if __name__ == "__main__":
    main()
