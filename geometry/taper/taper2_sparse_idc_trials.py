"""Tune couplings on the body with native complex-fit RMS 0.00155.

Its explicit-CAP inductance is cached: 4.330243614 nH. No body changes.
"""
from dataclasses import replace
from geometry.taper.taper2_ground_pad_trials import base as measured
from geometry.taper.parallel_trials import run

base=replace(measured,coupling=replace(measured.coupling,d=100,d_t1=32,d_regular=74))
CASES=[
    ('idc_n24fg20_right6_dl256_dr42_pad148_w152_p594_i100_j32_r74tg2',base,
     'Parent sparse-IDC fit RMS 0.00155: fitted input Cg 117.45 fF is low, '
     'T2/T1 Cg 36.92 fF is high. Increase finger length 50 to 98 um '
     'without increasing rail width or finger count; increase T2/T1 gap '
     '8 to 32 um. Restore centered regular junction gap74 used by the '
     'best previous target response. Body and independent L are identical.'),
    ('idc_n24fg20_right6_dl256_dr42_pad148_w152_p594_i100_j40_r74tg2',
     replace(base,coupling=replace(base.coupling,d_t1=40)),
     'Controlled companion: T2/T1 gap32 to40 um only, to measure its '
     'effect on the outer coupling, effective T1 shunt C and actual peak balance.')]

if __name__=='__main__':
    run(CASES,'sparse_IDC_input100_gap32_vs40')
