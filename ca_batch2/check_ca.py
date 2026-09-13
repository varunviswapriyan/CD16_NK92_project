#!/usr/bin/env python
"""
check_ca.py  -  plot + summarize a finished Ca fit. Run inside a Ca fit folder (needs Ca_model.csv;
uses ca_result.json and the pZAP csv if present). No ODE re-solve: Ca_model.csv already holds
model vs exp at every experimental time point.
    python check_ca.py [out.png]
"""
import sys, os, json
import numpy as np, pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

out = sys.argv[1] if len(sys.argv) > 1 else 'ca_plot.png'
m = pd.read_csv('Ca_model.csv')
res = json.load(open('ca_result.json')) if os.path.exists('ca_result.json') else {}
pz_path = res.get('pzap', 'optimized_model_pzap/model_output_pzap.csv')
C = {'zeta': '#1f77b4', 'hetero': '#ff7f0e', 'gamma': '#2ca02c'}
lab = {'zeta': 'CD3z', 'hetero': 'Hetero', 'gamma': 'FceRIg'}

print(f"== {res.get('note', '')}")
if res: print('params: ' + '  '.join(f'{k}={v:.4g}' for k, v in res['linear'].items()) + '  fixed=' + str({k: v for k, v in res['fixed'].items() if k not in res['params']}))
tot = 0.0
for k in ['zeta', 'gamma', 'hetero']:
    d = m[m.curve == k]
    ssr = float(np.sum(d.residual ** 2)); tot += ssr
    rms = np.sqrt(ssr / len(d)); rel = rms / d.exp_ca.mean()
    print(f'  {k:6s} SSR={ssr:.4e}  RMS={rms:.2f} ({rel * 100:.1f}% of mean)  peak exp={d.exp_ca.max():.1f} @ {d.time_abs_s[d.exp_ca.idxmax()]:.0f}s  model={d.model_ca.max():.1f} @ {d.time_abs_s[d.model_ca.idxmax()]:.0f}s  end exp={d.exp_ca.iloc[-1]:.1f} model={d.model_ca.iloc[-1]:.1f}')
print(f'  total SSR={tot:.4e}')

fig, (a1, a2) = plt.subplots(1, 2, figsize=(16, 6))
if os.path.exists(pz_path):
    pz = pd.read_csv(pz_path); tc = 'time' if 'time' in pz.columns else 'time_seconds'
    for k in C:
        col = f'mean_pZAP_{k}' if f'mean_pZAP_{k}' in pz.columns else f'mean_{k}'
        a1.plot(pz[tc], pz[col] / (25 * 602), color=C[k], lw=2, label=lab[k])
    a1.set_xlabel('Time (s)'); a1.set_ylabel('pZAP input (uM)'); a1.set_title('pZAP input to Ca ODE'); a1.legend(); a1.grid(alpha=.3)
for k in C:
    d = m[m.curve == k]
    a2.plot(d.time_abs_s, d.exp_ca, color=C[k], lw=1, alpha=.5)
    a2.scatter(d.time_abs_s, d.exp_ca, color=C[k], s=8, alpha=.5, label=f'{lab[k]} exp')
    a2.plot(d.time_abs_s, d.model_ca, color=C[k], lw=2.5, label=f'{lab[k]} model')
a2.set_xlabel('Time (s)'); a2.set_ylabel('Ca signal'); a2.set_title(f'Ca fit   total SSR={tot:.3e}'); a2.legend(ncol=2, fontsize=9); a2.grid(alpha=.3)
fig.suptitle(res.get('note', 'Ca fit') + '   ' + '  '.join(f'{k}={v:.3g}' for k, v in res.get('linear', {}).items()), fontsize=11)
plt.tight_layout(); plt.savefig(out, dpi=120); print('Saved', out)
