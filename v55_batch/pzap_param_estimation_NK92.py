"""
V55 - pZAP70 (Tyr493) fitting in NK92 cells with PSO  (N-parameter version)
Based on Indrani Nayak's pzap_param_estimation_NK92.py (Aug 31, 2026) and the V35 script.

What's new vs. V35:
  * Any number of fitted parameters. The first four are always the originals
    [lig0, kd10, ZAP0, SYK0]; any further names are .bngl parameters (multipliers) that
    are set on the model before each simulation via the 'extra' hook in calculate_SSR_*.py.
  * A parameter name may be prefixed with a model kind to apply to one model only,
    e.g. "gamma:KDL_MULT" (per-condition degradation rate) or "mixed:adaptor0" (hetero density).
  * FIXED: extra .bngl parameters held constant for the whole run (same naming rules).
  * All fitted values are log10-transformed; bounds are log10.

Config keys (v_config.json): PARAMS, LB, UB, FIXED, N_REPS, HETERO_WEIGHT, ORDER_WEIGHT,
                             RATIO_WEIGHT, SHAPE_MODE, PSO_OPTIONS, FIT_TMAX
Usage:  python pzap_param_estimation_NK92.py <n_particles> <n_iterations>
Output: analysis_param_residue.dat  (check_fit.py reads the PARAMS header line)
"""

import sys
import math
import time
import multiprocessing
from multiprocessing import Pool

import numpy as np
import pandas as pd

from pyswarms.single.global_best import GlobalBestPSO
from calculate_SSR_gamma import gamma_PZAP, calculate_residue
from calculate_SSR_zeta import zeta_PZAP
from calculate_SSR_mixed import mixed_PZAP


# ============================== CONFIG (defaults; v_config.json overrides) ==============================
PARAMS        = ['lig0', 'kd10', 'ZAP0', 'SYK0']     # first four are fixed in meaning
LB            = [+1.4, -4.5, +2.0, +0.5]              # log10
UB            = [+2.4, -2.3, +3.2, +2.5]              # log10
FIXED         = {}                                    # e.g. {"KZBG_FRAC": 0.5}
N_REPS        = 3
HETERO_WEIGHT = 1.0
ORDER_WEIGHT  = 5.0
RATIO_WEIGHT  = 30.0
SHAPE_MODE    = "own"
PSO_OPTIONS   = {'c1': 1.5, 'c2': 1.5, 'w': 0.5}
FIT_TMAX      = 300.0

import os, json
_cfg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'v_config.json')
if os.path.exists(_cfg_path):
    with open(_cfg_path) as _f:
        _cfg = json.load(_f)
    PARAMS        = _cfg.get('PARAMS', PARAMS)
    LB            = _cfg.get('LB', LB)
    UB            = _cfg.get('UB', UB)
    FIXED         = _cfg.get('FIXED', FIXED)
    N_REPS        = _cfg.get('N_REPS', N_REPS)
    HETERO_WEIGHT = _cfg.get('HETERO_WEIGHT', HETERO_WEIGHT)
    ORDER_WEIGHT  = _cfg.get('ORDER_WEIGHT', ORDER_WEIGHT)
    RATIO_WEIGHT  = _cfg.get('RATIO_WEIGHT', RATIO_WEIGHT)
    SHAPE_MODE    = _cfg.get('SHAPE_MODE', SHAPE_MODE)
    PSO_OPTIONS   = _cfg.get('PSO_OPTIONS', PSO_OPTIONS)
    FIT_TMAX      = _cfg.get('FIT_TMAX', FIT_TMAX)
    print(f'Loaded config from {_cfg_path}')

BASE4 = ['lig0', 'kd10', 'ZAP0', 'SYK0']
assert PARAMS[:4] == BASE4, "PARAMS must start with lig0, kd10, ZAP0, SYK0"
assert len(PARAMS) == len(LB) == len(UB), "PARAMS/LB/UB length mismatch"
# ========================================================================================================


def extras_for(kind, values):
    """Build the extra-parameter dict for one model kind from FIXED + fitted values.
    Keys 'gamma:X' apply only to the gamma model; unprefixed keys apply to all models."""
    out = {}
    for src in (FIXED, values):
        for k, v in src.items():
            if k in BASE4:
                continue
            if ':' in k:
                mk, name = k.split(':', 1)
                if mk == kind:
                    out[name] = v
            else:
                out[k] = v
    return out


def _interp_sorted(tgrid, t, y):
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    order = np.argsort(t)
    return np.interp(tgrid, t[order], y[order])


def process_particle(args):
    i, param, outdir, exp = args
    vals = {name: float(param[i, j]) for j, name in enumerate(PARAMS)}
    lig0, kd10, ZAP0, SYK0 = (vals[k] for k in BASE4)
    dir_name = str(outdir) + str(i)

    try:
        print(f'directory_name={dir_name}')
        tg, pg = gamma_PZAP(N_REPS, lig0, kd10, ZAP0, SYK0, dir_name, extra=extras_for('gamma', vals))
        tz, pz = zeta_PZAP(N_REPS, lig0, kd10, ZAP0, SYK0, dir_name, extra=extras_for('zeta', vals))
        tm, pm = mixed_PZAP(N_REPS, lig0, kd10, ZAP0, SYK0, dir_name, extra=extras_for('mixed', vals))

        pg = np.asarray(pg, dtype=float)
        pz = np.asarray(pz, dtype=float)
        pm = np.asarray(pm, dtype=float)
        g_max = max(float(np.max(pg)), 1e-12)
        z_max = max(float(np.max(pz)), 1e-12)
        m_max = max(float(np.max(pm)), 1e-12)

        if SHAPE_MODE == "own":
            r_g = calculate_residue(tg, pg / g_max, exp['t'], exp['g_own'])
            r_z = calculate_residue(tz, pz / z_max, exp['t'], exp['z_own'])
            r_m = calculate_residue(tm, pm / m_max, exp['t'], exp['m_own'])
        else:
            r_g = calculate_residue(tg, pg / g_max, exp['t'], exp['g_gam'])
            r_z = calculate_residue(tz, pz / g_max, exp['t'], exp['z_gam'])
            r_m = calculate_residue(tm, pm / g_max, exp['t'], exp['m_gam'])
        residue_total = r_g + r_z + HETERO_WEIGHT * r_m

        if ORDER_WEIGHT > 0:
            tgrid = np.linspace(20.0, FIT_TMAX, 100)
            gi = _interp_sorted(tgrid, tg, pg / g_max)
            zi = _interp_sorted(tgrid, tz, pz / g_max)
            mi = _interp_sorted(tgrid, tm, pm / g_max)
            viol = np.maximum(0.0, gi - mi) + np.maximum(0.0, zi - gi)
            residue_total += ORDER_WEIGHT * float(np.mean(viol ** 2))

        if RATIO_WEIGHT > 0:
            model_hg = m_max / g_max
            model_zg = z_max / g_max
            residue_total += RATIO_WEIGHT * ((model_hg - exp['hg']) ** 2 + (model_zg - exp['zg']) ** 2)

        extras_str = ' '.join(f'{k}={v:.3g}' for k, v in vals.items() if k not in BASE4)
        print(f'particle {i}: r_g={r_g:.4f} r_z={r_z:.4f} r_m={r_m:.4f} '
              f'hg={m_max / g_max:.3f}(exp {exp["hg"]:.3f}) zg={z_max / g_max:.3f}(exp {exp["zg"]:.3f}) '
              f'total={residue_total:.4f}  [{extras_str}]')

    except Exception as e:
        print(f"Error processing particle {i}: {e}")
        residue_total = float('inf')

    return residue_total


def cost_func(param):
    param = 10 ** param
    args = [(i, param, outdir, EXP) for i in range(param.shape[0])]
    n_cpus = multiprocessing.cpu_count()
    print(f"number of CPUs available: {n_cpus}")
    with Pool(processes=n_cpus) as pool:
        ssr_total = pool.map(process_particle, args)
    return np.asarray(ssr_total, dtype=float)


def main():
    global outdir, EXP
    outdir = 'analysis'
    number_of_particles = int(sys.argv[1])
    iteration_number = int(sys.argv[2])

    print("=== V55 CONFIG ===")
    print(f"PARAMS={PARAMS}")
    print(f"LB={LB}  UB={UB}  FIXED={FIXED}")
    print(f"N_REPS={N_REPS}  HETERO_WEIGHT={HETERO_WEIGHT}  ORDER_WEIGHT={ORDER_WEIGHT}  "
          f"RATIO_WEIGHT={RATIO_WEIGHT}  SHAPE_MODE={SHAPE_MODE}  PSO_OPTIONS={PSO_OPTIONS}")

    data = pd.read_csv('data/pZAP70_Tyr493_mean.csv')
    exp_time = data['time'].to_numpy(dtype=float)
    mask = exp_time <= FIT_TMAX
    t = exp_time[mask]
    z = data['mean_zeta'].to_numpy(dtype=float)[mask]
    g = data['mean_gamma'].to_numpy(dtype=float)[mask]
    m = data['mean_hetero'].to_numpy(dtype=float)[mask]
    g_max = max(float(np.max(g)), 1e-12)
    z_max = max(float(np.max(z)), 1e-12)
    m_max = max(float(np.max(m)), 1e-12)
    EXP = {'t': t,
           'g_gam': g / g_max, 'z_gam': z / g_max, 'm_gam': m / g_max,
           'g_own': g / g_max, 'z_own': z / z_max, 'm_own': m / m_max,
           'hg': m_max / g_max, 'zg': z_max / g_max}
    print(f"Experimental time points: {t}")
    print(f"Experimental peak ratios: hetero/gamma={EXP['hg']:.4f}  zeta/gamma={EXP['zg']:.4f}")

    start_time = time.time()
    optimizer = GlobalBestPSO(n_particles=number_of_particles, dimensions=len(PARAMS),
                              options=PSO_OPTIONS, bounds=(np.array(LB), np.array(UB)))
    residue, optimized_param = optimizer.optimize(cost_func, iters=iteration_number)

    print(optimized_param)
    print(residue)
    xx = math.sqrt(residue)
    yy = math.sqrt(xx)
    print(f' residue ** (1/4) = {yy}')

    optimized_param = np.asarray(optimized_param)
    with open(str(outdir) + "_param_residue.dat", 'w') as f:
        f.write("The optimized parameter is " + str(optimized_param) + "\n")
        f.write("Residue = " + str(residue) + "\n")
        f.write("Square root of residue = " + str(xx) + "\n")
        f.write("fourth root of residue = " + str(yy) + "\n")
        f.write("outdir\t" + str(outdir) + "\n")
        f.write("\t".join(PARAMS) + "\tresidue_1_4th\n")
        f.write("\t".join(str(v) for v in optimized_param) + "\t" + str(yy) + "\n")
        f.write(f"config\tPARAMS={PARAMS} FIXED={FIXED} N_REPS={N_REPS} HETERO_WEIGHT={HETERO_WEIGHT} "
                f"ORDER_WEIGHT={ORDER_WEIGHT} RATIO_WEIGHT={RATIO_WEIGHT} SHAPE_MODE={SHAPE_MODE}\n")
        f.write("linear\t" + "  ".join(f"{n}={10 ** v:.4g}" for n, v in zip(PARAMS, optimized_param)) + "\n")

    print(f'The time taken for the process is {(time.time() - start_time) / 60.0} min')
    print('The process ends now')


if __name__ == "__main__":
    main()
