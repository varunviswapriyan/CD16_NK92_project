"""ci_report.py -- the final pZAP CI table, with containment enforced.

The problem this exists to solve: the reference fit fell OUTSIDE its own
percentile interval.  That is not something a different formula papers over --
it means the bootstrap refits landed somewhere the reference fit did not, and
the first job is to say WHY.

So this prints, in order:

  1. CONVERGENCE DIAGNOSTICS.  If the bootstrap refits are under-converged
     relative to the reference fit, every interval below is measuring
     optimizer scatter as though it were sampling variability.  The tell is
     median(bootstrap SSR) / reference SSR much above 1, and a large spread of
     SSR across samples.  Fix that before trusting any interval.

  2. THREE INTERVALS per parameter, with containment marked:
       percentile  exact quantiles of the weighted bootstrap distribution.
                   Standard, tightest, but NOT guaranteed to contain the point
                   estimate.
       BC          bias-corrected percentile.  Shifts the quantile levels by
                   the weighted fraction of samples below the point estimate.
                   UNDEFINED when the point estimate is outside the whole
                   bootstrap support -- which is exactly the failure case, so
                   BC is not a rescue, only a refinement.
       normal      point +/- 1.96 * weighted SD.  Always contains the point
                   estimate by construction.  Wider, and assumes a normal
                   shape that a 10-atom distribution does not have.

  3. A REPORTABLE choice per parameter under a stated rule, plus the
     alternative of quoting the bootstrap median as the point estimate, which
     is contained by construction.

Because all 10 distinct day-resamples are enumerated and weighted by their
exact multinomial probabilities, the bootstrap distribution is known exactly:
a discrete distribution on 10 atoms.  Its quantiles are exact, not estimated,
so the usual "you need B >= 1000" caveat does not apply.

    BOOT_TAG=_conv python ~/CD16_NK92_project/filesCC/ci_report.py
    BOOT_TMAX=600  python ~/CD16_NK92_project/filesCC/ci_report.py
"""
import os, sys, json
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bootstrap_pzap as BP           # reuse its root, weights and parser

SUB   = BP.SUB
ROOT  = BP.BOOT_ROOT
FIT6  = BP.FIT6
# multipliers -> physical units
PHYS  = {'KZP_MULT': ('kzp', 0.03, '(uM s)^-1'),
         'KPR_MULT': ('KPR', 0.01, 's^-1')}

# convergence thresholds; a run that fails these makes the intervals unsafe
MAX_SSR_RATIO  = float(os.environ.get('CI_MAX_SSR_RATIO', '1.5'))
MAX_SSR_SPREAD = float(os.environ.get('CI_MAX_SSR_SPREAD', '10'))


def ncdf(x):
    from math import erf, sqrt
    return 0.5 * (1.0 + erf(x / sqrt(2.0)))


def nppf(p):
    """Inverse normal CDF (Acklam); good to ~1e-9, plenty here."""
    if p <= 0.0 or p >= 1.0:
        return float('-inf') if p <= 0 else float('inf')
    a = [-3.969683028665376e+01, 2.209460984245205e+02, -2.759285104469687e+02,
         1.383577518672690e+02, -3.066479806614716e+01, 2.506628277459239e+00]
    b = [-5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
         6.680131188771972e+01, -1.328068155288572e+01]
    c = [-7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
         -2.549732539343734e+00, 4.374664141464968e+00, 2.938163982698783e+00]
    d = [7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
         3.754408661907416e+00]
    pl, ph = 0.02425, 1 - 0.02425
    if p < pl:
        q = np.sqrt(-2 * np.log(p))
        return (((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
               ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    if p > ph:
        q = np.sqrt(-2 * np.log(1 - p))
        return -(((((c[0]*q+c[1])*q+c[2])*q+c[3])*q+c[4])*q+c[5]) / \
                ((((d[0]*q+d[1])*q+d[2])*q+d[3])*q+1)
    q = p - 0.5; r = q * q
    return (((((a[0]*r+a[1])*r+a[2])*r+a[3])*r+a[4])*r+a[5])*q / \
           (((((b[0]*r+b[1])*r+b[2])*r+b[3])*r+b[4])*r+1)


def wquant(x, w, q):
    """Quantile of a discrete distribution with atoms x and weights w."""
    x = np.asarray(x, float); w = np.asarray(w, float)
    o = np.argsort(x); x, w = x[o], w[o]
    c = np.cumsum(w) / w.sum()
    i = int(np.searchsorted(c, q, side='left'))
    return float(x[min(i, len(x) - 1)])


def load():
    """-> (reference dict, [(weight, ssr, values)], coverage)"""
    ref, rows, cov = None, [], 0.0
    r0 = BP.parse_result(os.path.join(ROOT, 'emp000'))
    if r0:
        ref = dict(ssr=r0[0], vals=r0[1])
    for i, ms in enumerate(BP.EMP_SETS, start=1):
        r = BP.parse_result(os.path.join(ROOT, 'emp%03d' % i))
        if not r:
            continue
        w = BP.emp_weight(ms)
        rows.append((w, r[0], r[1])); cov += w
    return ref, rows, cov


def bounds_of():
    """log10 search box from the reference sample's config, for !bound flags."""
    try:
        c = json.load(open(os.path.join(ROOT, 'emp000', SUB, 'v_config.json')))
    except Exception:
        return {}
    return {BP.canon(n): (c['LB'][j], c['UB'][j]) for j, n in enumerate(c['PARAMS'])}


def main():
    print('root: %s\n' % ROOT)
    ref, rows, cov = load()
    if ref is None or len(rows) < 3:
        raise SystemExit('not enough finished samples yet (reference=%s, samples=%d)'
                         % (ref is not None, len(rows)))

    W  = np.array([r[0] for r in rows], float)
    SS = np.array([r[1] if r[1] is not None else np.nan for r in rows], float)
    Wn = W / W.sum()

    # ---------------- 1. convergence ----------------------------------------
    print('=' * 74)
    print('1. CONVERGENCE  -- are the refits as well optimised as the reference?')
    print('=' * 74)
    print('  day-sets fitted      : %d of 10   (%.0f%% of bootstrap probability)'
          % (len(rows), 100 * cov))
    ok_conv = True
    if ref['ssr'] and np.isfinite(SS).any():
        med = wquant(SS[np.isfinite(SS)], Wn[np.isfinite(SS)], 0.5)
        ratio = med / ref['ssr']
        spread = np.nanmax(SS) / np.nanmin(SS)
        print('  reference SSR        : %.4e' % ref['ssr'])
        print('  median bootstrap SSR : %.4e   ratio %.2f   (target <= %.1f)'
              % (med, ratio, MAX_SSR_RATIO))
        print('  SSR spread max/min   : %.1fx                (target <= %.0fx)'
              % (spread, MAX_SSR_SPREAD))
        if ratio > MAX_SSR_RATIO:
            ok_conv = False
            print('  --> FAIL: refits are landing at much worse optima than the')
            print('            reference.  The intervals below then measure optimizer')
            print('            scatter, not sampling variability.  Raise the PSO')
            print('            budget (BOOT_PARTICLES / BOOT_ITERS) and rerun.')
        if spread > MAX_SSR_SPREAD:
            ok_conv = False
            print('  --> FAIL: SSR varies too much across resamples of the same 3 days.')
        if ok_conv:
            print('  --> PASS: refits are converged comparably to the reference.')
    print()

    # ---------------- 2. intervals ------------------------------------------
    BND = bounds_of()
    print('=' * 74)
    print('2. INTERVALS   (c = contains the point estimate)')
    print('=' * 74)
    print('%-9s %10s %10s  %-22s %-22s %-22s' %
          ('param', 'point', 'boot med', 'percentile', 'BC', 'normal +/-1.96SD'))
    print('-' * 118)

    chosen, allrows = {}, {}
    for n in FIT6:
        x = np.array([r[2].get(n, np.nan) for r in rows], float)
        m = np.isfinite(x)
        if m.sum() < 3 or n not in ref['vals']:
            print('%-9s %s' % (n, '(not parsed)')); continue
        xs, ws = x[m], Wn[m] / Wn[m].sum()
        th = float(ref['vals'][n])
        mu = float(np.sum(ws * xs))
        sd = float(np.sqrt(np.sum(ws * (xs - mu) ** 2)))
        med = wquant(xs, ws, 0.5)

        pc = (wquant(xs, ws, 0.025), wquant(xs, ws, 0.975))
        nm = (th - 1.96 * sd, th + 1.96 * sd)
        frac = float(np.sum(ws[xs < th]))
        if 0.0 < frac < 1.0:
            z0 = nppf(frac)
            bc = (wquant(xs, ws, ncdf(2 * z0 - 1.959964)),
                  wquant(xs, ws, ncdf(2 * z0 + 1.959964)))
        else:
            bc = None                       # point estimate outside the support

        def fmt(ci):
            if ci is None:
                return '%-22s' % '  undefined'
            c = 'c' if ci[0] <= th <= ci[1] else ' '
            return ('%s[%.4g, %.4g]' % (c, ci[0], ci[1])).ljust(22)

        print('%-9s %10.4g %10.4g  %-22s %-22s %-22s'
              % (n, th, med, fmt(pc), fmt(bc), fmt(nm)))
        allrows[n] = dict(th=th, med=med, sd=sd, pc=pc, bc=bc, nm=nm, frac=frac,
                          xs=xs, ws=ws)

        # selection rule, stated once below the table
        if pc[0] <= th <= pc[1]:
            chosen[n] = ('percentile', pc)
        elif bc is not None and bc[0] <= th <= bc[1]:
            chosen[n] = ('BC', bc)
        else:
            chosen[n] = ('normal', nm)
    print('-' * 118)
    print('rule: percentile if it contains the point estimate, else BC, else normal.')
    print('BC is undefined when the point estimate lies outside every bootstrap')
    print('sample -- it corrects a shifted distribution, it cannot rescue an')
    print('estimate outside the support.')
    print()

    # ---------------- 3. reportable -----------------------------------------
    print('=' * 74)
    print('3. REPORTABLE TABLE')
    print('=' * 74)
    print('%-9s %12s %26s %7s  %s' % ('param', 'estimate', '95% CI', 'w', 'form'))
    print('-' * 74)
    nbad = 0
    for n in FIT6:
        if n not in chosen:
            continue
        form, ci = chosen[n]
        a = allrows[n]; th = a['th']
        w = (ci[1] - ci[0]) / abs(th) if th else float('inf')
        flag = ''
        lo, hi = BND.get(n, (None, None))
        if lo is not None:
            lg = np.log10(np.abs(a['xs']))
            if np.any(np.abs(lg - lo) < 5e-3) or np.any(np.abs(lg - hi) < 5e-3):
                flag = ' !bound'
        if form != 'percentile':
            nbad += 1
        print('%-9s %12.4g %26s %7.2f  %s%s'
              % (n, th, '[%.4g, %.4g]' % ci, w, form, flag))
    print('-' * 74)
    for n, (nm_, mult, unit) in PHYS.items():
        if n in chosen:
            _, ci = chosen[n]
            print('  %-4s = %.4g  95%% CI [%.4g, %.4g]  %s'
                  % (nm_, allrows[n]['th'] * mult, ci[0] * mult, ci[1] * mult, unit))
    print()

    # ---------------- verdict ------------------------------------------------
    print('=' * 74)
    if not ok_conv:
        print('VERDICT: convergence FAILED.  Do not quote these intervals yet --')
        print('rerun with a larger PSO budget, then re-read this table.')
    elif nbad == 0:
        print('VERDICT: every parameter uses the standard percentile interval and')
        print('contains its point estimate.  This table is reportable.')
    else:
        print('VERDICT: %d parameter(s) needed a fallback interval because the' % nbad)
        print('percentile interval did not contain the point estimate.  Those are')
        print('honest but wider.  The alternative is to quote the BOOTSTRAP MEDIAN')
        print('as the estimate, which the percentile interval contains by')
        print('construction:')
        print()
        print('  %-9s %12s %26s %7s' % ('param', 'boot median', '95% CI (percentile)', 'w'))
        for n in FIT6:
            if n not in allrows:
                continue
            a = allrows[n]; pc = a['pc']; md = a['med']
            print('  %-9s %12.4g %26s %7.2f'
                  % (n, md, '[%.4g, %.4g]' % pc, (pc[1] - pc[0]) / abs(md)))
        print()
        print('  Quoting the median is defensible -- with a stochastic simulator a')
        print('  single optimisation run is itself a noisy estimator -- but the')
        print('  fitted CURVES in the figures come from the reference fit, so the')
        print('  figures would no longer correspond to the quoted numbers.')
    print('=' * 74)


if __name__ == '__main__':
    main()
