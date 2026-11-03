"""Full native response only after the isolated open-meander L screen passed."""
import json
from geometry.taper.taper2_open_meander_trials import open_case,open_settings
from geometry.taper import taper2_plain as w

if __name__=='__main__':
    measured=w.measure(open_case,open_settings)
    if abs(measured['target_error_percent'])>2:
        raise ValueError('Open meander failed the independent L screen')
    score,fitted=w.response(open_case,open_settings,measured['L0_nH'])
    print(json.dumps(dict(target=score,fit=fitted),indent=2),flush=True)
    if score['pass'] and fitted['curve_pass']:w.finalize(open_case)
