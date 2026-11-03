"""Parameter-only wider-gap meander: explicit L gate before full response."""
import json
from dataclasses import replace,asdict
from geometry.taper.taper2_open_meander_trials import open_settings as earlier
from geometry.taper import taper2_plain as w

SETTINGS=replace(earlier,meander=replace(earlier.meander,w=200,n=37,d=12))
CASE='idc_n70fg4_right30_dl256_dr122_pad148_w200_n37fg12_p100_i100_j60_r82tg2'

if __name__=='__main__':
    project,out,_=w.prepare(CASE,SETTINGS)
    (out/'proposal.json').write_text(json.dumps(dict(case=CASE,settings=asdict(SETTINGS),
        reason='Inductance-gated wider-gap ordinary-meander test:37 runs with12um '
        'gaps,4um trace and596um total length. Reduce inter-run distributed '
        'capacitance while keeping the approved cap/IDC/ground architecture. '
        'Initial width168um measured3.69480nH, so increase width to200um '
        'without changing the trace thickness. No response runs unless independent L is '
        'within2percent. T1 and regular/output bodies remain unchanged.'),indent=2)+'\n')
    measured=w.measure(CASE,SETTINGS)
    print(json.dumps(dict(stage='inductance_check',**measured),indent=2),flush=True)
    if abs(measured['target_error_percent'])>2:raise ValueError('Independent L failed; response skipped')
    score,fitted=w.response(CASE,SETTINGS,measured['L0_nH'])
    print(json.dumps(dict(target=score,fit=fitted),indent=2),flush=True)
    if score['pass'] and fitted['curve_pass']:w.finalize(CASE)
