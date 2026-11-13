#!/usr/bin/env python
"""Plot the pZAP70 best-fit model output against the experimental data. Read-only.

    python plot_pzap_fit.py
    python plot_pzap_fit.py --model optimized_model_pzap/model_output_pzap_v77.csv \
                            --data  estimate_params_pzap/data/<file>.xlsx

Column names are auto-detected, and whatever is found is printed. If the data file
is not given, the script looks for one and otherwise plots the model alone.
"""
import argparse, glob, os, re, sys
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

CONDS = ['zeta', 'gamma', 'hetero']
COLOR = {'zeta': '#0000CC', 'gamma': '#C00000', 'hetero': '#2BB8C4'}
DATA_GLOBS = ['estimate_params_pzap/data/*.xlsx', 'estimate_params_pzap/data/*.csv',
              'data/*.xlsx', 'data/*.csv']


def read_any(path):
    if path.lower().endswith(('.xlsx', '.xlsm', '.xls')):
        try:
            import openpyxl  # noqa: F401
        except ImportError:
            sys.exit('openpyxl needed for xlsx:  pip install --user openpyxl')
        return pd.read_excel(path)
    return pd.read_csv(path)


def find_time(df):
    for c in df.columns:
        if re.fullmatch(r'\s*time(_s(econds)?)?\s*', str(c), re.I):
            return c
    for c in df.columns:
        if 'time' in str(c).lower():
            return c
    return df.columns[0]


def find_cond(df, cond, kind):
    """kind: 'mean' or 'se'. Score columns containing the condition name."""
    want = {'mean': ('mean', 'avg', 'pzap'), 'se': ('se', 'sem', 'std', 'err')}[kind]
    best, score = None, 0
    for c in df.columns:
        s = str(c).lower()
        if cond not in s:
            continue
        v = 1 + sum(w in s for w in want)
        if kind == 'mean' and any(w in s for w in ('se', 'sem', 'std', 'err')):
            continue
        if v > score:
            best, score = c, v
    return best


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--model', default='optimized_model_pzap/model_output_pzap_v77.csv')
    ap.add_argument('--data', default=None)
    ap.add_argument('--dir', default='.')
    ap.add_argument('--out', default='fit_pzap.png')
    a = ap.parse_args()
    os.chdir(a.dir)

    if not os.path.exists(a.model):
        sys.exit(f'model output not found: {os.path.abspath(a.model)}\n'
                 f'run from the dir that holds optimized_model_pzap/, or pass --model')
    m = read_any(a.model)
    tm = find_time(m)
    print(f'model : {a.model}\n        time column -> {tm!r}')

    d = None
    if a.data is None:
        for g in DATA_GLOBS:
            hit = sorted(glob.glob(g))
            if hit:
                a.data = hit[0]
                print(f'        (auto-selected data file: {a.data}; {len(hit)} matched)')
                break
    if a.data and os.path.exists(a.data):
        d = read_any(a.data)
        td = find_time(d)
        print(f'data  : {a.data}\n        time column -> {td!r}\n        columns: {list(d.columns)}')
    else:
        print('data  : none found -- plotting model only')

    fig, axes = plt.subplots(1, 3, figsize=(11, 3.4), dpi=130, sharex=True)
    for ax, k in zip(axes, CONDS):
        cm = find_cond(m, k, 'mean')
        if cm is not None:
            ax.plot(m[tm], m[cm], color=COLOR[k], lw=1.9, label=f'model ({cm})')
        else:
            ax.text(.5, .5, f'no model column for {k}', ha='center', transform=ax.transAxes)
        if d is not None:
            dm = find_cond(d, k, 'mean'); ds = find_cond(d, k, 'se')
            if dm is not None:
                yerr = d[ds] if ds is not None else None
                ax.errorbar(d[td], d[dm], yerr=yerr, fmt='o', ms=4, lw=0.9,
                            color='0.35', capsize=2, label=f'data ({dm})')
        ax.set_title(k, fontsize=10); ax.set_xlabel('Time')
        ax.legend(fontsize=6.5, frameon=False)
        for sp in ('top', 'right'): ax.spines[sp].set_visible(False)
    axes[0].set_ylabel('pZAP70')
    fig.suptitle('pZAP70 best fit (v77) vs experiment', fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(a.out)
    print(f'wrote {os.path.abspath(a.out)}')


if __name__ == '__main__':
    main()
