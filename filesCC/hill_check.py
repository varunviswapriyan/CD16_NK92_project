"""hill_check.py -- does this experiment actually probe the Hill threshold?

The calcium drive in Indrani's ODE is

    F = z**nn / (z**nn + k3**nn)   +   k4 * z
        \_______ Hill _______/        \_ linear _/

with z = pZAP per ITAM.  The Hill term is a steep switch (nn = 9): it is ~0
well below k3, ~1 well above, and only carries information where z crosses k3.

Removing it costs no SSR.  That has two very different explanations:

  (a) there is no threshold in the biology, or
  (b) z never crosses k3 in this data, so the term is a CONSTANT here and
      the switch may be entirely real outside the range probed.

Only (b) is checkable without new experiments, and it is checkable with no
fitting at all -- just compare the range of z to k3.  This script does that.

    cd ~/Ca_fit_c02
    python ~/CD16_NK92_project/filesCC/hill_check.py
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import bootstrap_ca2 as B                      # reuse its loader, never re-guess the CSV format

K3 = float(os.environ.get('HILL_K3', B.N01['k3']))
NN = float(B.BASEFIX['nn'])


def hill(z):
    return z**NN / (z**NN + K3**NN)


def main():
    for f in (B.CA_PATH, B.PZAP_PATH):
        if os.path.isdir(f):
            raise SystemExit('ERROR: %s is a directory, not a file.' % f)
        if not os.path.isfile(f):
            raise SystemExit('ERROR: %s not found (run this from ~/Ca_fit_c02)' % f)

    base = B.load_base()
    print('Hill threshold k3 = %.4g   steepness nn = %g' % (K3, NN))
    print('z = pZAP per ITAM, over the fitted time window\n')
    print('%-8s %6s %11s %11s %11s %11s %9s' %
          ('adaptor', 'ITAMs', 'min z', 'max z', 'Hill(min)', 'Hill(max)', 'crosses?'))
    print('-' * 76)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    spans, any_cross = [], False
    for k in B.CONDS:
        d = base[k]
        zc = d['in_raw'] / d['itam']
        z = np.interp(d['t'], d['in_t'], zc, left=zc[0], right=zc[-1])
        h = hill(z)
        lo, hi = float(z.min()), float(z.max())
        # "crosses" = the threshold sits inside the range the data actually visits
        cross = (lo < K3 < hi)
        any_cross = any_cross or cross
        spans.append((k, lo, hi, float(h.min()), float(h.max()), cross))
        print('%-8s %6g %11.4g %11.4g %11.4f %11.4f %9s' %
              (B.DISPLAY[k], d['itam'], lo, hi, h.min(), h.max(),
               'YES' if cross else 'no'))

        ax1.plot(d['t'], z, color=B.COLORS[k], lw=2, label=B.DISPLAY[k])
        ax2.plot(d['t'], h, color=B.COLORS[k], lw=2, label=B.DISPLAY[k])

    ax1.axhline(K3, color='k', ls='--', lw=1.5, label='k3 threshold')
    ax1.set_ylabel('z = pZAP per ITAM'); ax1.set_yscale('log')
    ax1.legend(fontsize=9); ax1.set_title('Does z cross the Hill threshold?')
    ax2.set_ylabel('Hill term'); ax2.set_xlabel('Time (s)')
    ax2.set_ylim(-0.05, 1.05); ax2.legend(fontsize=9)
    ax2.set_title('Hill term over time -- a flat line means it acts as a constant')
    fig.tight_layout()
    out = os.path.join(B.ROOT if os.path.isdir(B.ROOT) else '.', 'hill_check.png')
    fig.savefig(out, dpi=110); plt.close(fig)

    print('-' * 76)
    hmin = min(s[3] for s in spans); hmax = max(s[4] for s in spans)
    swing = hmax - hmin
    print('Hill term ranges over %.4f - %.4f across all adaptors (swing %.4f)\n'
          % (hmin, hmax, swing))

    if not any_cross and swing < 0.05:
        where = 'above' if spans[0][1] > K3 else 'below'
        print('VERDICT: the data never crosses k3 -- z stays %s it throughout,' % where)
        print('so the Hill term is effectively a CONSTANT over this experiment.')
        print('It contributes no shape, which is why removing it costs no SSR.')
        print()
        print('This is a statement about the EXPERIMENT, not about the biology.')
        print('Say to Indrani: the threshold is not probed by these three')
        print('conditions, so k3 cannot be estimated and keeping the term makes')
        print('the remaining parameters unidentifiable.  Do NOT say the')
        print('switch does not exist -- this data cannot address that.')
        print()
        print('To probe it you would need conditions spanning z ~ %.3g,' % K3)
        print('i.e. a weaker or stronger stimulus than any of the three here.')
    elif any_cross and swing > 0.2:
        print('VERDICT: z DOES cross k3 and the Hill term varies substantially.')
        print('The term is being exercised by this data, so dropping it is a')
        print('real change to the model -- and the fact that SSR does not move')
        print('deserves a closer look before removing it.  Report full1 or')
        print('k3fix rather than lin1, and tell me this came out.')
    else:
        print('VERDICT: borderline -- partial crossing, or a small swing.')
        print('Treat the Hill removal as a modelling choice that needs stating')
        print('explicitly, not as something the data settles.')

    print('\nfigure: %s' % os.path.abspath(out))


if __name__ == '__main__':
    main()
