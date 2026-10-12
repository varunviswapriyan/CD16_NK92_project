#!/usr/bin/env python3
"""
compare_ca_csv.py -- is a newly parsed Ca CSV the same data the fits were run on?

  python3 compare_ca_csv.py <old.csv> <new.csv>

Run this BEFORE you drop anything into ca_data/Ca_NK92.csv. If the means differ,
every SSR already quoted (orig = 5.206e4 etc.) is stale and the fits have to be
re-run before those numbers go in an email.

Exit status: 0 = identical within tolerance, 1 = differs.
"""
import signal
import sys

import numpy as np
import pandas as pd

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

TOL = 1e-6


def main():
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    old, new = sys.argv[1], sys.argv[2]
    a = pd.read_csv(old)
    b = pd.read_csv(new)

    print(f'old: {old}   {a.shape[0]} rows, {a.shape[1]} cols')
    print(f'new: {new}   {b.shape[0]} rows, {b.shape[1]} cols')

    ca, cb = set(a.columns), set(b.columns)
    if ca != cb:
        print(f'\ncolumns differ:')
        print(f'  only in old: {sorted(ca - cb)}')
        print(f'  only in new: {sorted(cb - ca)}')

    tcol = 'time_seconds' if 'time_seconds' in ca & cb else list(ca & cb)[0]
    m = a.merge(b, on=tcol, how='outer', suffixes=('_old', '_new'), indicator=True)
    only_old = (m['_merge'] == 'left_only').sum()
    only_new = (m['_merge'] == 'right_only').sum()
    if only_old or only_new:
        print(f'\ntime grid differs: {only_old} bins only in old, {only_new} only in new')

    shared = m[m['_merge'] == 'both']
    worst = 0.0
    print()
    for c in sorted(ca & cb):
        if c == tcol:
            continue
        x = shared[f'{c}_old'].to_numpy(float)
        y = shared[f'{c}_new'].to_numpy(float)
        ok = np.isfinite(x) & np.isfinite(y)
        if not ok.any():
            continue
        d = np.abs(x[ok] - y[ok])
        rel = d / np.maximum(np.abs(x[ok]), 1e-12)
        worst = max(worst, d.max())
        flag = '' if d.max() <= TOL else '   <-- DIFFERS'
        print(f'  {c:16s} max|diff|={d.max():12.6g}  max rel={rel.max():10.4g}{flag}')

    same = worst <= TOL and not only_old and not only_new and ca == cb
    print()
    if same:
        print('IDENTICAL -- safe to swap in; existing SSRs still valid.')
        return 0
    print('DIFFERENT -- do NOT quote the old SSRs. Re-run the fits after swapping.')
    return 1


if __name__ == '__main__':
    sys.exit(main())
