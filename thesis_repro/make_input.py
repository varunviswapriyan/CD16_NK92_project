#!/usr/bin/env python3
"""
make_input.py  -- build an in.simulation from Giuseppe's in.full.

It reproduces EXACTLY the transformation in his run_samples.py
(first_changes_to_input_file + the seed edit), so the resulting input
file is byte-for-byte the model/params he ran, with only:
  - the random seed set,
  - read_sites pointed at the chosen day-14 slide,
  - (female only) his female rate edits applied,
  - (optional) AR+ IFN-g uptake set to 0 for the no-sequestration case.

Usage:
  python3 make_input.py --infull in.full --out in.simulation \
      --seed 8377 --sex male   --slide d14_M_4_SPPARKS_init
  python3 make_input.py --infull in.full --out in.simulation \
      --seed 8377 --sex female --slide d14_M_4_madeF_SPPARKS_init
  add --no-seq to zero the AR+ IFN-g uptake (sequestration) rate.

Female factor is Giuseppe's: make_female = 0.9*0.7 = 0.63, so
AR- recruitment is multiplied by 1/0.63, AR+ recruitment -> 0,
Tp->Tm differentiation x2, Tp->AR+Tm differentiation -> 0.
"""
import argparse

MAKE_FEMALE = 0.9 * 0.7          # 0.63, exactly as in run_samples.py
FEMALE_DIFF_FACTOR = 2


def female_line(ls):
    """Port of first_changes_to_input_file's make_female branch, verbatim logic."""
    recruit_factor = 1.0 / MAKE_FEMALE
    new_line = ls
    if ls.find('recruitment') >= 0 and ls.find('app_style') < 0:
        if ls.find('AR+') >= 0:                        # AR+ recruitment -> 0
            value = float(ls.split()[5])
            new_line = ls.replace(str(value), str(0))
        elif ls.find('Tp recruitment') >= 0 or ls.find('Tm recruitment') >= 0:
            value = ls.split()[5]                      # male recruitment value
            new_value = recruit_factor * float(value)
            new_line = ls.replace(str(value), str(new_value))
    if ls.find('Tp differentiation') >= 0 and ls.find('app_style') < 0:
        if ls.find('AR+ Tm') >= 0:                     # Tp->AR+Tm diff -> 0
            value = float(ls.split()[5])
            new_line = ls.replace(str(value), str(0))
        else:                                          # Tp->Tm diff x2
            value = ls.split()[5]
            new_value = FEMALE_DIFF_FACTOR * float(value)
            new_line = ls.replace(str(value), str(new_value))
    return new_line


def zero_ar_uptake(ls):
    """On the active app_style line, set the AR+ IFN-g uptake rate (token after
    app_style's 13th numeric arg) to 0. In in.full that is the 0.0016667 field."""
    # strip trailing comment
    code, sep, comment = ls.partition('#')
    toks = code.split()
    # toks[0]=app_style toks[1]=rxn/diff/custom; uptake is toks[15]
    if len(toks) > 15:
        toks[15] = '0.0'
    return '    '.join(['app_style', toks[1]]) + ' ' + ' '.join(toks[2:]) + \
           (('  #' + comment) if sep else '') + ('\n' if ls.endswith('\n') else '')


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--infull', required=True)
    ap.add_argument('--out', required=True)
    ap.add_argument('--seed', type=int, required=True)
    ap.add_argument('--sex', choices=['male', 'female'], required=True)
    ap.add_argument('--slide', required=True)
    ap.add_argument('--no-seq', action='store_true')
    a = ap.parse_args()

    out = []
    with open(a.infull) as f:
        for i, ls in enumerate(f):
            new = ls
            # seed line (in.full: line index 2 is "seed    8377")
            if ls.split() and ls.split()[0] == 'seed':
                old = ls.split()[1]
                new = ls.replace(old, str(a.seed))
            # read_sites -> chosen slide
            elif ls.find('read_sites') >= 0:
                new = "read_sites       sex_bias_initialization_files/%s\n" % a.slide
            # active app_style: optionally zero AR+ uptake
            elif a.no_seq and ls.lstrip().startswith('app_style'):
                new = zero_ar_uptake(ls)
            # female rate edits
            if a.sex == 'female':
                new = female_line(new)
            out.append(new)

    with open(a.out, 'w') as f:
        f.writelines(out)
    print("wrote %s  (seed=%d sex=%s slide=%s no_seq=%s)"
          % (a.out, a.seed, a.sex, a.slide, a.no_seq))


if __name__ == '__main__':
    main()
