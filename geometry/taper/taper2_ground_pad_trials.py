"""Parameter-only alternative suggested by the fixed-background circuit check.

The specified target circuit and acceptance window are unchanged. Effective
fitted capacitances are guidance only. Both IDCs occupy the same 556 um height
as the approved 70-finger IDC: 24*4 + 23*20 = 556 um.
"""
from dataclasses import replace
from geometry.taper.taper2_idc import SETTINGS as baseline
from geometry.taper.parallel_trials import run

base=replace(baseline,cap=replace(baseline.cap,t_right=6),
    meander=replace(baseline.meander,w=152),
    clearance=replace(baseline.clearance,d_right=42),
    pad=replace(baseline.pad,l=594),
    coupling=replace(baseline.coupling,d=52,n=24,d_finger=20,
                     d_t1=8,d_regular=58))
CASES=[
    ('idc_n24fg20_right6_dl256_dr42_pad148_w152_p594_i52_j8_r58tg2',base,
     'Circuit-only fixed-T1-C0 diagnostic suggests an alternative equivalent response '
     'with higher T2 C0 and lower input Cg. Increase shunt C using the existing '
     'right ground pad 14 um from the coil (formerly 94), reduce IDC finger count '
     '70 to 24 at identical occupied height, and strengthen both downstream '
     'junctions. This is a heuristic; nominal targets are unchanged and native '
     'S21 RMS<=0.025 plus explicit independent L must still pass.'),
    ('idc_n24fg20_right6_dl256_dr58_pad148_w152_p594_i52_j8_r58tg2',
     replace(base,clearance=replace(base.clearance,d_right=58)),
     'Controlled companion: right coil-to-cap clearance 42 to 58 um, moving '
     'the existing right ground pad from 14 to 30 um from the coil. '
     'Measure its shunt-C and inductance effect before choosing further changes.')]

if __name__=='__main__':
    run(CASES,'sparse_IDC_close_right_pad42_vs58')
