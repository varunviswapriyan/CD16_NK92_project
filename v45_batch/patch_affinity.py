#!/usr/bin/env python3
"""
patch_affinity.py  <file.bngl>  <fraction>  <rule_name> [<rule_name> ...]

Adds parameter   kzb_g kzb * <fraction>   right after the 'kzb' line, and
switches the listed ZAP-binding rules from 'kzb, kzu' to 'kzb_g, kzu'.
Exits non-zero if the parameter line could not be inserted or any rule was not patched.
"""
import sys, re
path, frac, rules = sys.argv[1], sys.argv[2], sys.argv[3:]
lines = open(path).read().split('\n')
inserted, patched = False, []
i = 0
while i < len(lines):
    ln = lines[i]
    if not inserted and re.match(r'^\s*kzb\s', ln):
        lines.insert(i + 1, f'kzb_g kzb * {frac}')
        inserted = True
        i += 2
        continue
    stripped = ln.lstrip()
    hit = next((r for r in rules if stripped.startswith(r + ':')), None)
    if hit:
        for j in range(i, min(i + 4, len(lines))):
            if re.search(r'\bkzb\s*,', lines[j]):
                lines[j] = re.sub(r'\bkzb\s*,', 'kzb_g,', lines[j], count=1)
                patched.append(hit)
                break
    i += 1
open(path, 'w').write('\n'.join(lines))
missing = [r for r in rules if r not in patched]
print(f'    affinity patch {path.split("/")[-1]}: kzb_g={frac}, rules patched {len(patched)}/{len(rules)}')
if not inserted or missing:
    print(f'    ERROR: inserted={inserted} missing={missing}')
    sys.exit(1)
