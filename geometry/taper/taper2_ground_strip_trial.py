"""Controlled ground-strip thickness test on the actual best native body."""
import json
from dataclasses import replace,asdict
from geometry.taper.taper2_idc import SETTINGS as baseline
from geometry.taper import taper2_plain as w

SETTINGS=replace(baseline,coupling=replace(baseline.coupling,t_regular_ground=6))
CASE='idc_rail8_right30_dl256_pad148_w148_p200_i42n70_j42_r74tg6'

if __name__=='__main__':
    project,out,_=w.prepare(CASE,SETTINGS)
    (out/'proposal.json').write_text(json.dumps(dict(case=CASE,settings=asdict(SETTINGS),
        reason='Only centered ground-strip thickness2 to6um at fixed74um '
        'T1/unit gap, on the best0.03947RMS body. Ambient ground can alter '
        'both mutual and shunt capacitances despite unchanged component bodies. '
        'Test actual peak balance; fixed extraction background is not proof '
        'that the physical unit ground capacitance remains constant.'),indent=2)+'\n')
    measured=w.measure(CASE,SETTINGS)
    if abs(measured['target_error_percent'])>2:raise ValueError('Independent L failed')
    score,fitted=w.response(CASE,SETTINGS,measured['L0_nH'])
    print(json.dumps(dict(target=score,fit=fitted),indent=2),flush=True)
    if score['pass'] and fitted['curve_pass']:w.finalize(CASE)
