#!/usr/bin/env python3
"""
plot_ca_data.py -- mean +/- SE Indo-1 ratio vs. time, so you can eyeball the parsed
CSV against Oscar's slide panels before anything gets fitted.

  python3 plot_ca_data.py <ca.csv> [<ca2.csv> ...] -o ca_data.png
"""
import argparse
import signal

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

CONDS = [('zeta', 'CD3z (zeta-zeta, 6 ITAM)', 'tab:blue'),
         ('gamma', 'FcRg (gamma-gamma, 2 ITAM)', 'tab:red'),
         ('hetero', 'CD3z+FcRg (hetero, 4 ITAM)', 'tab:green'),
         ('parental', 'NK92 parental', '0.5')]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('csvs', nargs='+')
    ap.add_argument('-o', '--out', default='ca_data.png')
    a = ap.parse_args()

    n = len(a.csvs)
    fig, axes = plt.subplots(1, n, figsize=(6.0 * n, 4.4), squeeze=False)
    for ax, path in zip(axes[0], a.csvs):
        df = pd.read_csv(path)
        t = df['time_seconds']
        for key, lab, col in CONDS:
            if f'mean_{key}' not in df.columns:
                continue
            m, s = df[f'mean_{key}'], df[f'SE_{key}']
            ax.plot(t, m, color=col, lw=1.4, label=lab)
            ax.fill_between(t, m - s, m + s, color=col, alpha=0.22, lw=0)
        ax.set_xlabel('time (s)')
        ax.set_ylabel('Indo-1 high/low ratio')
        ax.set_title(path.split('/')[-1], fontsize=9)
        ax.legend(fontsize=7, frameon=False)
    fig.tight_layout()
    fig.savefig(a.out, dpi=150)
    print(f'wrote {a.out}')


if __name__ == '__main__':
    main()
