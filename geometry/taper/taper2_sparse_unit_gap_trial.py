"""Small, controlled centered-junction change; no body or L change."""
import json
from dataclasses import replace,asdict
from geometry.taper.taper2_open_meander_trials import input_settings
from geometry.taper import taper2_plain as w

SETTINGS=replace(input_settings,coupling=replace(input_settings.coupling,d_regular=78))
CASE='idc_n24fg20_right6_dl256_dr42_pad148_w152_p594_i160_j40_r78tg2'

if __name__=='__main__':
    project,out,_=w.prepare(CASE,SETTINGS)
    (out/'proposal.json').write_text(json.dumps(dict(case=CASE,settings=asdict(SETTINGS),
        reason='Only centered T1/unit gap74 to78um. Parent lower/upper peaks '
        '4.8870/5.1091GHz are below4.8938/5.1118GHz targets. Extracted '
        'unit coupling8.509fF exceeds8.2865fF. A slightly larger gap tests '
        'frequency correction and supplies an independent native coupling sensitivity.'),indent=2)+'\n')
    measured=w.measure(CASE,SETTINGS)
    if abs(measured['target_error_percent'])>2:raise ValueError('Independent L failed')
    score,fitted=w.response(CASE,SETTINGS,measured['L0_nH'])
    print(json.dumps(dict(target=score,fit=fitted),indent=2),flush=True)
    if score['pass'] and fitted['curve_pass']:w.finalize(CASE)
