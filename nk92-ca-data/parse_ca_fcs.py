#!/usr/bin/env python3
"""
parse_ca_fcs.py -- Indo-1 ratio vs. time -> binned mean + SE, from Oscar's
exported (already FlowJo-gated) FCS 3.0 files.

Command-line equivalent of Indrani's Ca_flow_analysis_mean_se.ipynb, with no
flowkit dependency (numpy + pandas only).

Usage
-----
  python3 parse_ca_fcs.py <fcs_dir> -o <out.csv> [options]

  --grid indrani   250 bins spanning 0 .. max time over ALL tubes  (HER default, default here)
  --grid fixed     fixed-width bins, set with --bin (seconds)
  --nbins 250      number of bins in indrani mode
  --bin 1.0        bin width in fixed mode
  --schema indrani time_seconds, mean_{zeta,gamma,hetero}, SE_{...}   (drop-in Ca_NK92.csv)
  --schema full    adds n_* and the parental tube
  --trim           drop bins where any of zeta/gamma/hetero has no data
  --minn 1         drop bins with fewer than this many cells (1 = keep everything, like her)

Tube -> model condition (ITAM count in parentheses):
  NK92 / Parental            -> parental   (no CD16 chimera)
  CD3z / KI-1 CD3z_CD16      -> zeta       (6 ITAMs)
  FcRg / KI-2 Fcg_CD16       -> gamma      (2 ITAMs)
  CD3z+FcRg / KI-6 ...       -> hetero     (4 ITAMs)

Time scale
----------
seconds = raw * $PnG * $TIMESTEP  (0.01 * 0.01 = 1e-4 in Oscar's files).

The gain factor is NOT optional. In the notebook it is applied by FlowKit inside
get_events(source="raw"); here it is applied explicitly. The check is $P7R, the
declared range of the Time parameter: Oscar's files say $P7R = 331, and
raw_max * 1e-4 = 331.3. Multiplying by $TIMESTEP alone gives 33133 s (9.2 h),
which contradicts $P7R and the ~5 min acquisition. This script asserts the two
agree and refuses to continue if they do not.
"""
import argparse
import os
import signal
import sys

import numpy as np
import pandas as pd

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

RATIO_CH = 'Indo-1 high_low Ratio'
TIME_CH = 'Time'

LABEL = {
    'parental': 'parental', 'nk92': 'parental',
    'cd3z': 'zeta', 'ki-1 cd3z_cd16': 'zeta', 'ki-1': 'zeta', 'zeta': 'zeta',
    'fcrg': 'gamma', 'ki-2 fcg_cd16': 'gamma', 'ki-2': 'gamma', 'gamma': 'gamma',
    'cd3z+fcrg': 'hetero', 'ki-6 fcg_cd3z_cd16': 'hetero', 'ki-6': 'hetero',
    'hetero': 'hetero',
}
FIT_CONDS = ['zeta', 'gamma', 'hetero']
ALL_CONDS = FIT_CONDS + ['parental']


# ----------------------------------------------------------------- FCS reader
def read_fcs(path):
    """Minimal FCS 3.0 list-mode reader. Returns (keyword dict, DataFrame)."""
    with open(path, 'rb') as fh:
        b = fh.read()
    head = b[:58].decode('ascii', 'replace')
    ts, te = int(head[10:18]), int(head[18:26])
    txt = b[ts:te + 1].decode('latin-1')
    d = txt[0]
    parts = txt.split(d)[1:]
    kv = dict(zip(parts[0::2], parts[1::2]))

    if kv.get('$MODE', 'L') != 'L':
        raise ValueError(f'{path}: $MODE={kv.get("$MODE")}, only list mode supported')
    if kv.get('$DATATYPE') != 'F':
        raise ValueError(f'{path}: $DATATYPE={kv.get("$DATATYPE")}, only float supported')

    n = int(kv['$PAR'])
    tot = int(kv['$TOT'])
    ds = int(kv.get('$BEGINDATA') or head[26:34])
    de = int(kv.get('$ENDDATA') or head[34:42])
    dt = '>f4' if kv['$BYTEORD'].startswith('4,3') else '<f4'
    a = np.frombuffer(b[ds:de + 1], dtype=dt, count=tot * n).reshape(tot, n)
    cols = [kv[f'$P{i}N'] for i in range(1, n + 1)]
    return kv, pd.DataFrame(a, columns=cols)


def time_in_seconds(kv, df, strict=True):
    """raw Time channel -> seconds, with the $PnR consistency check."""
    tidx = list(df.columns).index(TIME_CH) + 1
    gain = float(kv.get(f'$P{tidx}G', 1.0) or 1.0)
    step = float(kv.get('$TIMESTEP', 1.0) or 1.0)
    raw = df[TIME_CH].to_numpy(float)
    t = raw * gain * step

    rng = kv.get(f'$P{tidx}R')
    note = f'gain={gain} timestep={step} -> raw_max {np.nanmax(raw):.0f} = {np.nanmax(t):.1f} s'
    if rng is not None:
        r = float(rng)
        note += f'  ($P{tidx}R={r:.0f})'
        if r > 0 and not (0.5 * r <= np.nanmax(t) <= 2.0 * r):
            msg = (f'time scale fails the $P{tidx}R check: converted max '
                   f'{np.nanmax(t):.1f} s vs declared range {r:.0f}. '
                   f'gain={gain}, timestep={step}.')
            if strict:
                raise ValueError(msg)
            print('WARNING: ' + msg, file=sys.stderr)
    return t, note


# ------------------------------------------------------------------- plumbing
def tube_label(path):
    base = os.path.basename(path)
    if '(' in base and ')' in base:
        tag = base[base.rfind('(') + 1:base.rfind(')')]
    else:
        tag = os.path.splitext(base)[0]
    key = tag.strip().lower()
    if key in LABEL:
        return LABEL[key], tag
    for k, v in sorted(LABEL.items(), key=lambda kv: -len(kv[0])):
        if k in key:
            return v, tag
    return None, tag


def find_fcs(root):
    out = []
    for dirpath, dirnames, filenames in os.walk(root):
        dirnames[:] = [x for x in dirnames if x != '__MACOSX']
        for f in filenames:
            if f.lower().endswith('.fcs') and not f.startswith('._'):
                out.append(os.path.join(dirpath, f))
    return sorted(out)


def load_events(path, strict=True):
    kv, df = read_fcs(path)
    for ch in (TIME_CH, RATIO_CH):
        if ch not in df.columns:
            raise ValueError(f'{path}: channel {ch!r} missing. Have: {list(df.columns)}')
    t, note = time_in_seconds(kv, df, strict=strict)
    r = df[RATIO_CH].to_numpy(float)
    ok = np.isfinite(t) & np.isfinite(r)
    return t[ok], r[ok], note


def bin_events(t, r, edges, minn):
    """mean / SE / n per bin, on the given edges. Empty bins -> NaN."""
    nb = len(edges) - 1
    k = np.digitize(t, edges) - 1
    ok = (k >= 0) & (k < nb)
    k, v = k[ok], r[ok]

    cnt = np.bincount(k, minlength=nb).astype(float)
    s1 = np.bincount(k, weights=v, minlength=nb)
    s2 = np.bincount(k, weights=v * v, minlength=nb)

    with np.errstate(invalid='ignore', divide='ignore'):
        mean = np.where(cnt > 0, s1 / cnt, np.nan)
        var = np.where(cnt > 1, (s2 - cnt * mean ** 2) / (cnt - 1), np.nan)
        se = np.sqrt(np.maximum(var, 0.0) / cnt)

    thin = cnt < max(minn, 1)
    mean[thin] = np.nan
    se[thin] = np.nan
    return mean, se, cnt


# ----------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('fcs_dir')
    ap.add_argument('-o', '--out', required=True)
    ap.add_argument('--grid', choices=['indrani', 'fixed'], default='indrani')
    ap.add_argument('--nbins', type=int, default=250)
    ap.add_argument('--bin', type=float, default=1.0)
    ap.add_argument('--schema', choices=['indrani', 'full'], default='indrani')
    ap.add_argument('--minn', type=int, default=1)
    ap.add_argument('--trim', action='store_true',
                    help='drop bins where any of zeta/gamma/hetero is empty')
    ap.add_argument('--long', default=None)
    ap.add_argument('--no-strict-time', action='store_true',
                    help='warn instead of failing if the $PnR time check does not pass')
    a = ap.parse_args()

    files = find_fcs(a.fcs_dir)
    if not files:
        sys.exit(f'no .fcs files under {a.fcs_dir}')

    ev = {}
    print(f'{"condition":10s} {"tube":14s} {"cells":>8s} {"t_max_s":>9s}  time scale')
    for p in files:
        cond, tag = tube_label(p)
        if cond is None:
            print(f'{"SKIP":10s} {tag:14s}')
            continue
        if cond in ev:
            sys.exit(f'two files map to condition {cond!r}; clean up {a.fcs_dir}')
        t, r, note = load_events(p, strict=not a.no_strict_time)
        ev[cond] = (t, r)
        print(f'{cond:10s} {tag:14s} {len(r):8d} {t.max():9.1f}  {note}')

    missing = [c for c in FIT_CONDS if c not in ev]
    if missing:
        sys.exit(f'missing model conditions: {missing} -- cannot build a fit-ready CSV')

    # --- bin edges. indrani mode: 250 bins over 0..max across ALL tubes loaded.
    gmax = max(t.max() for t, _ in ev.values())
    if a.grid == 'indrani':
        edges = np.linspace(0.0, gmax, a.nbins + 1)
        print(f'\ngrid: indrani -- {a.nbins} bins over 0..{gmax:.2f} s '
              f'(width {gmax / a.nbins:.4f} s)')
    else:
        nb = int(np.ceil(gmax / a.bin))
        edges = np.arange(nb + 1) * a.bin
        print(f'\ngrid: fixed -- {nb} bins of {a.bin} s over 0..{edges[-1]:.2f} s')
    centers = (edges[:-1] + edges[1:]) / 2

    out = pd.DataFrame({'time_seconds': centers})
    counts = {}
    for c in ALL_CONDS:
        if c not in ev:
            continue
        m, s, n = bin_events(ev[c][0], ev[c][1], edges, a.minn)
        out[f'mean_{c}'] = m
        out[f'SE_{c}'] = s
        counts[c] = n

    # --- report empties, which is where Exp 2 bites
    print()
    for c in FIT_CONDS:
        bad = int(np.isnan(out[f'mean_{c}']).sum())
        if bad:
            t_bad = out.time_seconds[np.isnan(out[f'mean_{c}'])]
            print(f'  {c:8s}: {bad} empty bins, {t_bad.min():.1f}-{t_bad.max():.1f} s '
                  f'(tube ends at {ev[c][0].max():.1f} s)')

    if a.trim:
        before = len(out)
        keep = ~out[[f'mean_{c}' for c in FIT_CONDS]].isna().any(axis=1)
        out = out[keep].reset_index(drop=True)
        for c in counts:
            counts[c] = counts[c][keep.to_numpy()]
        print(f'  --trim: dropped {before - len(out)} bins, '
              f'kept {out.time_seconds.min():.1f}-{out.time_seconds.max():.1f} s')

    # --- column schema
    if a.schema == 'indrani':
        cols = ['time_seconds'] + [f'mean_{c}' for c in FIT_CONDS] + [f'SE_{c}' for c in FIT_CONDS]
        out = out[cols]
    else:
        for c in ALL_CONDS:
            if c in counts:
                out[f'n_{c}'] = counts[c].astype(int)
        cols = ['time_seconds']
        for c in ALL_CONDS:
            if c in counts:
                cols += [f'n_{c}', f'mean_{c}', f'SE_{c}']
        out = out[cols]

    os.makedirs(os.path.dirname(os.path.abspath(a.out)) or '.', exist_ok=True)
    out.to_csv(a.out, index=False)
    print(f'\nwrote {a.out}   {len(out)} rows x {len(out.columns)} cols  (schema={a.schema})')
    print(f'   time range: {out.time_seconds.min():.2f} - {out.time_seconds.max():.2f} s')
    for c in FIT_CONDS:
        m = out[f'mean_{c}']
        base = m[out.time_seconds < 20].mean()
        print(f'   {c:8s} baseline(<20s)={base:8.3f}  peak={m.max():8.3f} '
              f'@ {out.time_seconds[m.idxmax()]:7.2f}s')

    if a.long:
        rows = []
        for c, (t, r) in ev.items():
            rows.append(pd.DataFrame({'condition': c, 'time_seconds': t, 'indo_ratio': r}))
        pd.concat(rows).to_csv(a.long, index=False)
        print(f'wrote {a.long} (event-level long format)')


if __name__ == '__main__':
    main()
