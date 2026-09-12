#!/usr/bin/env python3
"""
patch_ssr_extra.py <calculate_SSR_x.py>
Adds an optional 'extra' dict argument to the *_PZAP function and applies it to model.parameters
right after SS0 is set.  Idempotent.
"""
import sys, re
p = sys.argv[1]; s = open(p).read()
if 'extra=None' in s:
    print(f'    {p.split("/")[-1]}: already patched'); sys.exit(0)
s2, n = re.subn(r'(def\s+\w+_PZAP\s*\([^)]*?dir_name)\s*\)\s*:', r'\1, extra=None):', s, count=1)
if n != 1: print(f'    ERROR: def line not found in {p}'); sys.exit(1)
lines = s2.split('\n'); out = []; done = False
for l in lines:
    out.append(l)
    if not done and re.match(r'^\s*model\.parameters\.SS0\s*=', l):
        ind = re.match(r'^(\s*)', l).group(1)
        out.append(f'{ind}for _k, _v in (extra or {{}}).items():')
        out.append(f'{ind}    setattr(model.parameters, _k, _v)')
        done = True
if not done: print(f'    ERROR: "model.parameters.SS0 =" line not found in {p}'); sys.exit(1)
open(p, 'w').write('\n'.join(out))
print(f'    {p.split("/")[-1]}: extra-params hook added')
