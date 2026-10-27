#!/usr/bin/env python3
"""
clean_fig.py -- the Ca fit figure, with a plain title and no shaded band.

    python3 clean_fig.py "Experiment 1" -o Experiment_1.png
    python3 clean_fig.py "Experiment 2" -o Experiment_2.png --variant orig

Run it from inside the fit directory (the one holding ca_data/ and out_boot_ca2/).
It imports bootstrap_ca2.py and reuses that file's own load_base/solve/rhs, so the
curve drawn here is the same curve `band` draws -- only the title and the shading
differ. Nothing is refitted.
"""
import argparse
import json
import os
import runpy
import signal
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt  # noqa: E402

signal.signal(signal.SIGPIPE, signal.SIG_DFL)

STYLE = [('zeta',   'CD3z',   'tab:blue'),
         ('gamma',  'FceRIg', 'tab:green'),
         ('hetero', 'Hetero', 'tab:orange')]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('title', help='plot title, e.g. "Experiment 1"')
    ap.add_argument('-o', '--out', default=None)
    ap.add_argument('--variant', default='orig')
    ap.add_argument('--boot', default=None,
                    help='path to bootstrap_ca2.py (default: ~/CD16_NK92_project/filesCC)')
    ap.add_argument('--sample', type=int, default=0,
                    help='which boot_NNN.json to draw (0 = the reference fit)')
    a = ap.parse_args()

    boot = a.boot or os.path.expanduser('~/CD16_NK92_project/filesCC/bootstrap_ca2.py')
    if not os.path.exists(boot):
        sys.exit(f'bootstrap_ca2.py not found at {boot} (pass --boot)')

    # import without running its CLI
    g = runpy.run_path(boot, run_name='_clean_fig')
    for need in ('VARIANTS', 'load_base', 'resample', 'solve', 'BASEFIX', 'N01'):
        if need not in g:
            sys.exit(f'{boot} has no {need} -- is this the right script?')

    if a.variant not in g['VARIANTS']:
        sys.exit(f'unknown variant {a.variant!r}; known: {" ".join(g["VARIANTS"])}')
    cfg = g['VARIANTS'][a.variant]

    jf = os.path.join(g.get('ROOT', 'out_boot_ca2'), a.variant, f'boot_{a.sample:03d}.json')
    if not os.path.exists(jf):
        sys.exit(f'no {jf} -- run this from the fit directory (the one with ca_data/)')
    d = json.load(open(jf))
    if not d.get('ok'):
        sys.exit(f'{jf} records a failed fit')

    # exactly how bootstrap_ca2 rebuilds a parameter set
    p = dict(g['BASEFIX']); p.update(g['N01']); p.update(cfg['fix'])
    p.update(d.get('fixed', {})); p.update(d.get('params', {}))

    data = g['resample'](g['load_base'](), seed=1, jitter=False)

    # solve() reads CFG/PN/DATA from the module namespace
    for fn in ('solve', 'rhs'):
        g[fn].__globals__.update({'CFG': cfg, 'PN': cfg['fit'], 'DATA': data})

    fig, ax = plt.subplots(figsize=(7.2, 4.8), constrained_layout=True)
    for key, label, col in STYLE:
        if key not in data:
            continue
        dd = data[key]
        ax.errorbar(dd['t'], dd['mean'], yerr=dd['se'], fmt='.', ms=2.5,
                    lw=0.6, color=col, alpha=0.55, label=f'{label} exp', zorder=1)
        ax.plot(dd['t'], g['solve'](dd, p, cfg), color=col, lw=2.4,
                label=f'{label} model', zorder=3)

    ax.set_title(a.title, fontsize=13)
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('Ca signal')
    h, l = ax.get_legend_handles_labels()
    order = [i for i, s in enumerate(l) if s.endswith('model')] + \
            [i for i, s in enumerate(l) if s.endswith('exp')]
    ax.legend([h[i] for i in order], [l[i] for i in order],
              ncol=2, fontsize=8, frameon=True)

    out = a.out or (a.title.replace(' ', '_') + '.png')
    fig.savefig(out, dpi=200)
    ssr = sum(float(np.sum((data[k]['mean'] - g['solve'](data[k], p, cfg)) ** 2))
              for k, _, _ in STYLE if k in data)
    print(f'wrote {out}')
    print(f'  variant={a.variant} sample={a.sample}  SSR={ssr:.4e} '
          f'(json says {d.get("ssr", float("nan")):.4e})')


if __name__ == '__main__':
    main()
