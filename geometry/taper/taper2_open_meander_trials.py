"""Input-coupling control plus an isolated, more open meander L check."""
from dataclasses import replace
from geometry.taper.taper2_sparse_idc_trials import base as sparse
from geometry.taper.taper2_idc import SETTINGS as original
from geometry.taper.parallel_trials import run

input_case='idc_n24fg20_right6_dl256_dr42_pad148_w152_p594_i160_j40_r74tg2'
open_case='idc_n70fg4_right30_dl256_dr122_pad148_w164_n58fg6_p100_i100_j60_r82tg2'
input_settings=replace(sparse,coupling=replace(sparse.coupling,d=160,d_t1=40))
open_settings=replace(original,
    meander=replace(original.meander,w=164,n=58,d=6),
    pad=replace(original.pad,l=100),
    coupling=replace(original.coupling,d=100,d_t1=60,d_regular=82))
CASES=[(input_case,input_settings,
    'Only input gap/finger length changes from i100,j40,r74: extracted input '
    'Cg99.27 fF is below153 fF. Test i160 while retaining all other geometry '
    'and reusing the explicit4.33024nH body calibration.'),
    (open_case,open_settings,
    'Inductance-only screen of a more open ordinary meander. Change n72/gap4 '
    'to n58/gap6 at the same596um total length,4um trace and centered leads. '
    'First width164um is an estimate, not a validated inductance. '
    'No full response is authorized by this job flag until L is checked.')]

if __name__=='__main__':
    run(CASES,'input160_and_open_meander164_Lonly',inductance_only_cases=[open_case])
