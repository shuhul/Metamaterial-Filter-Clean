"""Edit these settings, then run with the local Python; no CLI arguments."""
from dataclasses import asdict
from hashlib import sha256
from pathlib import Path
import json
import sys
from contextlib import redirect_stdout

ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from geometry.meander import MeanderParams
from geometry.taper.layout import TaperParams,build,validate
from geometry.taper.experiment import run_case,measure_unit
from geometry.taper.fitting import fit

# Change CASE_NAME for each new simulated design. Completed trials are retained.
CASE_NAME='finger536_t20_right310_input6_j56'
SIMULATE_INDUCTANCE=False
SIMULATE_MIXED=False
FIT_RESPONSE=True

TAPER=TaperParams(
    w=328,t_top=14,t_left=24,t_right=6,
    l_left=584,l_right=310,
    meander=MeanderParams(w=240,l=596,l_top=10,n=58,t=4,d=6),
)
D_INPUT=6
D_JUNCTION=56
INPUT_STYLE='finger'  # Native taper-unit fixture uses a finger feed, not an end meander.
L_INPUT=536
T_INPUT=20
INDUCTANCE_RESULT_FILE=ROOT/'output_data/taper_exploration/finger536_t20_right310_input6_j56_isolation/results.json'


def main():
    reference=INDUCTANCE_RESULT_FILE
    if SIMULATE_INDUCTANCE:reference,_=measure_unit(CASE_NAME,TAPER)
    if SIMULATE_MIXED:
        result=run_case(CASE_NAME,TAPER,d_input=D_INPUT,d_junction=D_JUNCTION,input_style=INPUT_STYLE,
                        l_input=L_INPUT,t_input=T_INPUT,inductance_reference=reference,fit_response=FIT_RESPONSE)
        print(result.get('snp'))
        return result
    if not FIT_RESPONSE:
        print('No mixed simulation or fit requested.');return
    # Refit only the response matching these settings; never silently fit stale geometry.
    folder=ROOT/'output_data/taper_exploration'/CASE_NAME
    trial_path=folder/'trial.json'
    project_dir=Path(json.loads(trial_path.read_text())['project']).parent if trial_path.exists() else ROOT/'sonnet_files/taper_exploration'/CASE_NAME
    gds=ROOT/'geometry/taper/output'/f'{CASE_NAME}.gds'
    metadata=json.loads(gds.with_suffix('.json').read_text())
    _,expected=build(TAPER,d_input=D_INPUT,d_junction=D_JUNCTION,input_style=INPUT_STYLE,l_input=L_INPUT,t_input=T_INPUT)
    for key in ('w_top_slot','l_top_slot'):
        if key not in metadata:expected.pop(key)
    if expected!=metadata:raise ValueError('Changed geometry: choose a new CASE_NAME and enable simulations')
    validate(gds,metadata)
    validation=json.loads((project_dir/'validation.json').read_text())
    if validation['source_gds_sha256']!=sha256(gds.read_bytes()).hexdigest():raise ValueError('Source GDS changed since import')
    status=json.loads((folder/'run_status.json').read_text())
    record=json.loads(trial_path.read_text()) if trial_path.exists() else {}
    project=Path(record['project']) if record else project_dir/'generated_ncell_2.sonx'
    snp=Path(record['snp']) if record else project.with_suffix('.snp')
    if status['status']!='complete' or status['project_sha256']!=sha256(project.read_bytes()).hexdigest() or status['snp_sha256']!=sha256(snp.read_bytes()).hexdigest():
        raise ValueError('Project/response changed since the recorded completed simulation')
    iso_metadata=json.loads((ROOT/'geometry/taper/output'/reference.parent.name).with_suffix('.json').read_text())
    if iso_metadata['taper']!=asdict(TAPER):raise ValueError('Inductance reference does not match taper; enable SIMULATE_INDUCTANCE')
    l=json.loads(reference.read_text())['L0_nH']
    cached=folder/'cached_fit';cached.mkdir(parents=True,exist_ok=True)
    with (cached/'fit_details.log').open('w') as log,redirect_stdout(log):
        result=fit(snp,cached,l_taper=l,l_regular=4.922723945068658,c_regular=188.8781163861302,free_regular=True)
    p=result['complex_with_delay']
    print(f"L={l:.6f} nH; C0={p['C_taper_fF']:.6f} fF; input Cg={p['Cg_input_fF']:.6f} fF; junction Cg={p['Cg_junction_fF']:.6f} fF")
    print('Plot:',folder/'cached_fit/s21_fit_complex_with_delay.png')
    print('Linear plot:',folder/'cached_fit/s21_fit_complex_with_delay_linear.png')
    return result

if __name__=='__main__':main()
