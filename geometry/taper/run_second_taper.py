"""Editable taper-2 trial settings. Default: render the saved candidate, no EM.

This candidate is NOT an accepted final design. See TAPER_2_README.md.
Always use a new CASE_NAME before generating a changed simulated design.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
from hashlib import sha256
from dataclasses import asdict
import json
import sys

ROOT=Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:sys.path.insert(0,str(ROOT))
from geometry.meander import MeanderParams
from geometry.taper.layout import TaperParams
from geometry.taper.second_layout import SecondParams,from_dict,validate,require_same_isolation_geometry
from geometry.taper.second_experiment import measure,run,fit_response,fit_network,plot_final
from geometry.taper.experiment import measure_unit

CASE_NAME='t2_plate276_slot264_tooth96_l400'
SIMULATE_INDUCTANCE=False
SIMULATE_TAPER1_INDUCTANCE=False
SIMULATE_MIXED=False
REFIT_RESPONSE=False
PLOT_SAVED_RESPONSE=True

GEOMETRY=SecondParams(
    taper=TaperParams(w=328,t_top=12,t_left=24,t_right=20,l_left=400,l_right=214,
                     meander=MeanderParams(w=100,l=598,l_top=10,n=4,t=32,d=6)),
    taper1=TaperParams(w=280,t_top=10,t_left=24,t_right=6,l_left=200,l_right=300,
                      meander=MeanderParams(w=152,l=600,l_top=10,n=72,t=2,d=6)),
    idc_pairs=15,idc_t=4,idc_d=4,idc_span=332,d_idc_top=12,
    t_ground_tail=44,d_t2_t1=6,t_coupling_finger=32,l_coupling_finger=182,
    n_coupling_teeth=3,w_coupling_teeth=96,w_upper_plate=276,
    l_coupling_bridge=222,l_cal_output=200,w_t1_ground_slot=264,l_t1_ground_slot=100,
    d_t1_regular=44,t_junction_ground=4,
)
INDUCTANCE_RESULT_FILE=ROOT/'output_data/taper_exploration/t2_plate276_slot264_tooth96_l400_isolation/results.json'
TAPER1_INDUCTANCE_RESULT_FILE=ROOT/'output_data/taper_exploration/t1_bar10_r300_slot264_isolation/results.json'


def main():
    reference=INDUCTANCE_RESULT_FILE
    reference1=TAPER1_INDUCTANCE_RESULT_FILE
    if SIMULATE_TAPER1_INDUCTANCE:
        if GEOMETRY.taper1 is None:raise ValueError('Set an embedded taper1 geometry before calibrating it')
        reference1,_=measure_unit(CASE_NAME+'_taper1',GEOMETRY.taper1,
                                 sweep_range=(3.8,6.4),
                                 w_top_slot=GEOMETRY.w_t1_ground_slot,
                                 l_top_slot=GEOMETRY.l_t1_ground_slot)
    if GEOMETRY.taper1 is not None and reference1 is None:
        raise ValueError('Embedded t1 variant requires a matching L calibration')
    if SIMULATE_INDUCTANCE:reference,_=measure(CASE_NAME,GEOMETRY)
    if SIMULATE_MIXED:return run(CASE_NAME,GEOMETRY,reference,inductance_t1_reference=reference1)
    if not (REFIT_RESPONSE or PLOT_SAVED_RESPONSE):return
    output=ROOT/'output_data/taper_exploration'/CASE_NAME
    gds=ROOT/'geometry/taper/output'/f'{CASE_NAME}.gds'
    metadata=json.loads(gds.with_suffix('.json').read_text())
    if from_dict(metadata['second'])!=GEOMETRY:raise ValueError('Changed geometry: use a new case name and simulate')
    validate(gds,metadata)
    status=json.loads((output/'run_status.json').read_text())
    for filename,key in [('response.snp','snp_sha256'),('project.sonx','project_sha256')]:
        if sha256((output/filename).read_bytes()).hexdigest()!=status[key]:raise ValueError('Archived simulation changed')
    trial=json.loads((output/'trial.json').read_text())
    if status['status']!='complete':raise ValueError('Simulation was not completed')
    recorded=Path(trial['inductance_reference']).resolve()
    if recorded!=reference.resolve():raise ValueError('Changed L reference')
    expected_hash=trial.get('hashes',{}).get(trial['inductance_reference'])
    if expected_hash and sha256(reference.read_bytes()).hexdigest()!=expected_hash:
        raise ValueError('Inductance record changed since fitting')
    ref_meta=json.loads((ROOT/'geometry/taper/output'/f'{reference.parent.name}.json').read_text())
    p=from_dict(ref_meta['second'])
    for key in ('taper','idc_pairs','idc_t','idc_d','idc_span','d_idc_top','x_meander','t_ground_tail',
                't_coupling_finger','l_coupling_finger','n_coupling_teeth','w_coupling_teeth',
                'w_upper_plate','l_coupling_bridge','l_cal_output','w_t1_ground_slot','l_t1_ground_slot'):
        if asdict(p)[key]!=asdict(GEOMETRY)[key]:raise ValueError('L reference geometry differs')
    require_same_isolation_geometry(GEOMETRY,ROOT/'geometry/taper/output'/f'{reference.parent.name}.gds')
    l2=json.loads(reference.read_text())['L0_nH']
    l1=6.014105932449249
    if GEOMETRY.taper1 is not None:
        reference1=Path(reference1)
        m1=json.loads((ROOT/'geometry/taper/output'/f'{reference1.parent.name}.json').read_text())
        if m1['taper']!=asdict(GEOMETRY.taper1):raise ValueError('T1 calibration geometry differs')
        if m1.get('w_top_slot',0)!=GEOMETRY.w_t1_ground_slot or m1.get('l_top_slot',100)!=GEOMETRY.l_t1_ground_slot:
            raise ValueError('T1 calibration ground slot differs')
        if Path(trial['inductance_t1_reference']).resolve()!=reference1.resolve():raise ValueError('Changed t1 L reference')
        expected_hash=trial['hashes'][trial['inductance_t1_reference']]
        if sha256(reference1.read_bytes()).hexdigest()!=expected_hash:raise ValueError('Changed t1 L record')
        l1=json.loads(reference1.read_text())['L0_nH']
    if REFIT_RESPONSE:
        cached=output/'cached_fit';cached.mkdir(exist_ok=True)
        result=fit_network(output/'response.snp',cached,l2,l1=l1)
        fit_response(output/'response.snp',cached,l2,l1=l1,seed_result=result)
    else:result=json.loads((output/'fit_network.json').read_text())
    path=output/'cached_fit/s21_fit_linear.png' if REFIT_RESPONSE else output/'s21_fit_linear.png'
    plot_final(output/'response.snp',result,path,title='Taper 2 Candidate: Targets Not Met')
    print('Candidate project:',trial['project'])
    print('Candidate plot:',path)
    print('This is an exploratory candidate, not a final target pass.')


if __name__=='__main__':main()
