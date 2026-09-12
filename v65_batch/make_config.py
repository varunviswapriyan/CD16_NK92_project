#!/usr/bin/env python3
"""make_config.py <dest_json> <extra_params|-> <extra_lb|-> <extra_ub|-> <fixed_json> <shape> <rw> <ow> <hw> <zaplb> <pso_json> <note>"""
import sys, json
(dest, ep, elb, eub, fixed, shape, rw, ow, hw, zaplb, pso, note) = sys.argv[1:13]
extra = [] if ep == '-' else ep.split(',')
elb = [] if elb == '-' else [float(x) for x in elb.split(',')]
eub = [] if eub == '-' else [float(x) for x in eub.split(',')]
assert len(extra) == len(elb) == len(eub), 'extra params/bounds mismatch'
cfg = {
    'PARAMS': ['lig0', 'kd10', 'ZAP0', 'SYK0'] + extra,
    'LB': [1.4, -4.5, float(zaplb), 0.5] + elb,
    'UB': [2.4, -2.3, 3.2, 2.5] + eub,
    'FIXED': json.loads(fixed),
    'N_REPS': 3,
    'HETERO_WEIGHT': float(hw), 'ORDER_WEIGHT': float(ow), 'RATIO_WEIGHT': float(rw),
    'SHAPE_MODE': shape, 'PSO_OPTIONS': json.loads(pso), 'FIT_TMAX': 300.0,
    'NOTE': note,
}
json.dump(cfg, open(dest, 'w'), indent=2)
print(f'    config: {len(cfg["PARAMS"])} params {cfg["PARAMS"][4:]} FIXED={cfg["FIXED"]} {shape} RW={rw} OW={ow} HW={hw}')
