"""
V35 - pZAP70 (Tyr493) fitting in NK92 cells with PSO
Based on Indrani Nayak's pzap_param_estimation_NK92.py (Aug 31, 2026).

Changes vs. original:
  1. SHAPE_MODE="own": each model curve is normalized by ITS OWN max and fit to the
     experimental curve normalized by ITS OWN max (pure shape fit). This removes the
     artificial ceiling where hetero could never exceed gamma because everything was
     divided by gamma's max.
  2. Absolute peak-ratio penalty: model peak(hetero)/peak(gamma) and peak(zeta)/peak(gamma)
     must match the experimental peak ratios (computed from the data file).
  3. Ordering penalty (hetero > gamma > zeta on 20-300 s), on gamma-normalized curves.
  4. HETERO_WEIGHT on the hetero shape residual.
  5. All knobs live in the CONFIG block at the top - edit those, nothing else.

Usage:  python pzap_param_estimation_NK92.py <n_particles> <n_iterations>
Output: analysis_param_residue.dat  (same format as before, so check_fit.py can read it)
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


# ============================== CONFIG (edit here only) ==============================
# log10 bounds for [lig0, kd10, ZAP0, SYK0]
LB = [+0.5, -3.3, +0.5, +0.5]
UB = [+1.8, -2.6, +3.2, +2.5]

N_REPS        = 3      # BioNetGen replicates averaged per particle
HETERO_WEIGHT = 1.0    # multiplier on hetero shape residual
ORDER_WEIGHT  = 5.0    # ordering penalty weight (hetero > gamma > zeta); 0 disables
RATIO_WEIGHT  = 30.0   # absolute peak-ratio penalty weight; 0 disables
SHAPE_MODE    = "own"  # "own" = each curve / its own max (new) ; "gamma" = original behaviour
PSO_OPTIONS   = {'c1': 1.5, 'c2': 1.5, 'w': 0.5}
FIT_TMAX      = 300.0  # use experimental points with t <= FIT_TMAX

# If a v_config.json exists next to this script, its keys override the defaults above.
import os, json
_cfg_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'v_config.json')
if os.path.exists(_cfg_path):
    with open(_cfg_path) as _f:
        _cfg = json.load(_f)
    LB            = _cfg.get('LB', LB)
    UB            = _cfg.get('UB', UB)
    N_REPS        = _cfg.get('N_REPS', N_REPS)
    HETERO_WEIGHT = _cfg.get('HETERO_WEIGHT', HETERO_WEIGHT)
    ORDER_WEIGHT  = _cfg.get('ORDER_WEIGHT', ORDER_WEIGHT)
    RATIO_WEIGHT  = _cfg.get('RATIO_WEIGHT', RATIO_WEIGHT)
    SHAPE_MODE    = _cfg.get('SHAPE_MODE', SHAPE_MODE)
    PSO_OPTIONS   = _cfg.get('PSO_OPTIONS', PSO_OPTIONS)
    FIT_TMAX      = _cfg.get('FIT_TMAX', FIT_TMAX)
    print(f'Loaded config from {_cfg_path}')
# =====================================================================================


def _interp_sorted(tgrid, t, y):
    """Interpolate y(t) onto tgrid, sorting t first (np.interp needs increasing x)."""
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    order = np.argsort(t)
    return np.interp(tgrid, t[order], y[order])


def process_particle(args):
    """Evaluate the cost for one particle."""
    i, param, outdir, exp = args

    lig0 = param[i, 0]
    kd10 = param[i, 1]
    ZAP0 = param[i, 2]
    SYK0 = param[i, 3]

    dir_name = str(outdir) + str(i)

    try:
        print(f'directory_name={dir_name}')

        tg, pg = gamma_PZAP(N_REPS, lig0, kd10, ZAP0, SYK0, dir_name)
        tz, pz = zeta_PZAP(N_REPS, lig0, kd10, ZAP0, SYK0, dir_name)
        tm, pm = mixed_PZAP(N_REPS, lig0, kd10, ZAP0, SYK0, dir_name)

        pg = np.asarray(pg, dtype=float)
        pz = np.asarray(pz, dtype=float)
        pm = np.asarray(pm, dtype=float)

        g_max = max(float(np.max(pg)), 1e-12)
        z_max = max(float(np.max(pz)), 1e-12)
        m_max = max(float(np.max(pm)), 1e-12)

        # ---------- 1. Shape residuals ----------
        if SHAPE_MODE == "own":
            r_g = calculate_residue(tg, pg / g_max, exp['t'], exp['g_own'])
            r_z = calculate_residue(tz, pz / z_max, exp['t'], exp['z_own'])
            r_m = calculate_residue(tm, pm / m_max, exp['t'], exp['m_own'])
        else:  # original: everything divided by gamma max
            r_g = calculate_residue(tg, pg / g_max, exp['t'], exp['g_gam'])
            r_z = calculate_residue(tz, pz / g_max, exp['t'], exp['z_gam'])
            r_m = calculate_residue(tm, pm / g_max, exp['t'], exp['m_gam'])

        residue_total = r_g + r_z + HETERO_WEIGHT * r_m

        # ---------- 2. Ordering penalty (on gamma-normalized = absolute relative) ----------
        if ORDER_WEIGHT > 0:
            tgrid = np.linspace(20.0, FIT_TMAX, 100)
            gi = _interp_sorted(tgrid, tg, pg / g_max)
            zi = _interp_sorted(tgrid, tz, pz / g_max)
            mi = _interp_sorted(tgrid, tm, pm / g_max)
            viol = np.maximum(0.0, gi - mi) + np.maximum(0.0, zi - gi)
            residue_total += ORDER_WEIGHT * float(np.mean(viol ** 2))

        # ---------- 3. Absolute peak-ratio penalty ----------
        if RATIO_WEIGHT > 0:
            model_hg = m_max / g_max
            model_zg = z_max / g_max
            residue_total += RATIO_WEIGHT * ((model_hg - exp['hg']) ** 2 + (model_zg - exp['zg']) ** 2)

        print(f'particle {i}: r_g={r_g:.4f} r_z={r_z:.4f} r_m={r_m:.4f} '
              f'hg={m_max / g_max:.3f}(exp {exp["hg"]:.3f}) zg={z_max / g_max:.3f}(exp {exp["zg"]:.3f}) '
              f'total={residue_total:.4f}')

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

    print("=== V35 CONFIG ===")
    print(f"LB={LB}  UB={UB}")
    print(f"N_REPS={N_REPS}  HETERO_WEIGHT={HETERO_WEIGHT}  ORDER_WEIGHT={ORDER_WEIGHT}  "
          f"RATIO_WEIGHT={RATIO_WEIGHT}  SHAPE_MODE={SHAPE_MODE}  PSO_OPTIONS={PSO_OPTIONS}")

    # ---------- experimental data ----------
    data = pd.read_csv('data/pZAP70_Tyr493_mean.csv')
    exp_time = data['time'].to_numpy(dtype=float)
    mean_zeta = data['mean_zeta'].to_numpy(dtype=float)
    mean_gamma = data['mean_gamma'].to_numpy(dtype=float)
    mean_hetero = data['mean_hetero'].to_numpy(dtype=float)

    mask = exp_time <= FIT_TMAX
    t = exp_time[mask]
    z = mean_zeta[mask]
    g = mean_gamma[mask]
    m = mean_hetero[mask]

    g_max = max(float(np.max(g)), 1e-12)
    z_max = max(float(np.max(z)), 1e-12)
    m_max = max(float(np.max(m)), 1e-12)

    EXP = {
        't': t,
        # original gamma-normalized targets
        'g_gam': g / g_max,
        'z_gam': z / g_max,
        'm_gam': m / g_max,
        # own-normalized (shape) targets
        'g_own': g / g_max,
        'z_own': z / z_max,
        'm_own': m / m_max,
        # experimental absolute peak ratios
        'hg': m_max / g_max,
        'zg': z_max / g_max,
    }
    print(f"Experimental time points: {t}")
    print(f"Experimental peak ratios: hetero/gamma={EXP['hg']:.4f}  zeta/gamma={EXP['zg']:.4f}")

    start_time = time.time()

    optimizer = GlobalBestPSO(n_particles=number_of_particles, dimensions=4,
                              options=PSO_OPTIONS, bounds=(LB, UB))
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
        f.write("lig0\tkd10\tZAP0\tSYK0\tresidue_1_4th\n")
        f.write(str(optimized_param[0]) + "\t" + str(optimized_param[1]) + "\t" +
                str(optimized_param[2]) + "\t" + str(optimized_param[3]) + "\t" + str(yy) + "\n")
        f.write(f"config\tLB={LB} UB={UB} N_REPS={N_REPS} HETERO_WEIGHT={HETERO_WEIGHT} "
                f"ORDER_WEIGHT={ORDER_WEIGHT} RATIO_WEIGHT={RATIO_WEIGHT} SHAPE_MODE={SHAPE_MODE}\n")

    print(f'The time taken for the process is {(time.time() - start_time) / 60.0} min')
    print('The process ends now')


if __name__ == "__main__":
    main()
