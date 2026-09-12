#!/usr/bin/env python3
"""
patch_model.py <file.bngl> <kind: zeta|gamma|mixed> [--frac] [--recdeg]

ALWAYS: inserts multiplier constants (all 1.0 = original behaviour) after the 'Vc' line and
        multiplies the matching rate lines by them, so the fitter can scale any of them:
          KZP_MULT -> kzp   (ZAP70 phosphorylation by Lck; controls rise)
          KPR_MULT -> KPR   (ligand off-rate / proofreading lifetime)
          KZU_MULT -> kzu   (ZAP70 unbinding from ITAM)
          KP_MULT  -> kp    (ITAM phosphorylation by Lck)
          KZD_MULT -> kzd   (SHP dephosphorylation of pZAP)
          KDL_MULT -> kdl   (CBL degradation; per-file so per-condition rates are possible)
          KZBG_FRAC        (gamma-type ITAM ZAP70 on-rate fraction; only used with --frac)
--frac   : add kzb_g = kzb*KZBG_FRAC and switch gamma-type ITAM ZAP70-binding rules to it
--recdeg : replace the 4 per-bound-kinase CBL rules with per-ligated-receptor rule(s)
Exits non-zero if any expected line is not found.
"""
import sys, re
path, kind = sys.argv[1], sys.argv[2]
flags = set(sys.argv[3:])
L = open(path).read().split('\n')
log = []
def die(msg):
    print(f'    ERROR ({path.split("/")[-1]}): {msg}'); sys.exit(1)

idx = next((i for i, l in enumerate(L) if re.match(r'^\s*Vc\s', l)), None)
if idx is None: die("no 'Vc' parameter line")
L[idx+1:idx+1] = ['KZBG_FRAC 1.0', 'KZP_MULT 1.0', 'KPR_MULT 1.0', 'KZU_MULT 1.0',
                  'KP_MULT 1.0', 'KZD_MULT 1.0', 'KDL_MULT 1.0']
for name, mult in [('kzp', 'KZP_MULT'), ('KPR', 'KPR_MULT'), ('kzu', 'KZU_MULT'),
                   ('kp', 'KP_MULT'), ('kzd', 'KZD_MULT'), ('kdl', 'KDL_MULT')]:
    hits = [i for i, l in enumerate(L) if re.match(rf'^\s*{name}\s', l)]
    if len(hits) != 1: die(f'expected 1 "{name}" param line, found {len(hits)}')
    L[hits[0]] = L[hits[0]].rstrip() + f' * {mult}'
log.append('multipliers on kzp,KPR,kzu,kp,kzd,kdl')

if '--frac' in flags:
    hits = [i for i, l in enumerate(L) if re.match(r'^\s*kzb\s', l)]
    if len(hits) != 1: die(f'expected 1 kzb param line, found {len(hits)}')
    L.insert(hits[0] + 1, 'kzb_g kzb * KZBG_FRAC')
    rules = {'zeta': [],
             'gamma': ['r9_PZAP_binding_ITAM1', 'r10_PZAP_binding_ITAM2',
                       'r7_UZAP_binding_ITAM1', 'r8_UZAP_binding_ITAM2'],
             'mixed': ['rh14_PZAP_binding_ITAM4', 'rh18_UZAP_binding_ITAM4',
                       'rg9_PZAP_binding_ITAM1', 'rg10_PZAP_binding_ITAM2',
                       'rg7_UZAP_binding_ITAM1', 'rg8_UZAP_binding_ITAM2']}[kind]
    done = []
    for i, l in enumerate(L):
        hit = next((r for r in rules if l.lstrip().startswith(r + ':')), None)
        if not hit: continue
        for j in range(i, min(i + 4, len(L))):
            if re.search(r'\bkzb\s*,', L[j]):
                L[j] = re.sub(r'\bkzb\s*,', 'kzb_g,', L[j], count=1); done.append(hit); break
    miss = [r for r in rules if r not in done]
    if miss: die(f'affinity rules not patched: {miss}')
    log.append(f'kzb_g on {len(done)} rules')

if '--recdeg' in flags:
    keep, removed = [], 0
    for l in L:
        if '_bound_CBL_degrade_' in l: removed += 1
        else: keep.append(l)
    if removed != 4: die(f'expected 4 CBL degrade rules, found {removed}')
    L = keep
    new = {'zeta':  ['rdeg_zeta: Zeta(receptor!+) -> ligand(receptor) kdl'],
           'gamma': ['rdeg_gamma: Gamma(receptor!+) -> ligand(receptor) kdl'],
           'mixed': ['rdeg_zeta: Zeta(receptor!+) -> ligand(receptor) kdl',
                     'rdeg_gamma: Gamma(receptor!+) -> ligand(receptor) kdl',
                     'rdeg_hetero: Hetero(receptor!+) -> ligand(receptor) kdl']}[kind]
    end = next((i for i, l in enumerate(L) if l.strip() == 'end reaction rules'), None)
    if end is None: die("no 'end reaction rules'")
    L[end:end] = new
    log.append(f'CBL: 4 per-kinase rules -> {len(new)} per-receptor rule(s)')

open(path, 'w').write('\n'.join(L))
print(f'    patched {path.split("/")[-1]}: ' + '; '.join(log))
