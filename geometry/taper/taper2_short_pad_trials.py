"""Native verification of measured-response proposals; no topology changes."""
from dataclasses import replace
from geometry.taper.taper2_idc import SETTINGS as baseline
from geometry.taper.parallel_trials import run

base = replace(baseline, cap=replace(baseline.cap, t_right=6),
    pad=replace(baseline.pad, l=100),
    coupling=replace(baseline.coupling, d=52, d_t1=36, d_regular=82))
CASES = [
    ('idc_rail8_right6_dl256_pad148_w148_p100_i52n70_j36_r82tg2', base,
     'Measured phase-aligned response proposal: shorten pads to lower T2 shunt C, '
     'lengthen input fingers and narrow T2/T1 gap to increase both input couplings; '
     'widen centered T1/unit junction to reduce its coupling. Predicted RMS 0.02378 '
     'is not acceptance; native fixed-target response and independent explicit-CAP L required.'),
    ('idc_rail8_right6_dl256_pad148_w148_p100_i52n70_j40_r82tg2',
     replace(base,coupling=replace(base.coupling,d_t1=40)),
     'Controlled companion to j36: increase only T2/T1 gap by 4 um. '
     'Check whether predicted coupling sensitivity survives native EM.')]

if __name__ == '__main__':
    run(CASES, 'shortpad100_input52_j36_vs40_r82')
