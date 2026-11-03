"""Audit taper-2 metal/settings, response identity and both L calibrations."""
from pathlib import Path
from dataclasses import asdict
from hashlib import sha256
import json
import numpy as np
from scipy.stats import linregress
from circuit_common import read_snp,single_sweep_resonance
from geometry.taper.audit import check_geometry
from geometry.taper.second_layout import validate,from_dict,require_same_isolation_geometry
from geometry.specifications import TAPER_2,TAPER_1,UNIT_CELL,specification

ROOT=Path(__file__).resolve().parents[2]


def calibration(path,*,second):
    path=Path(path).resolve();r=json.loads(path.read_text())
    gds=ROOT/'geometry/taper/output'/f'{path.parent.name}.gds'
    meta=json.loads(gds.with_suffix('.json').read_text())
    xy=[]
    for point in r['sweep']:
        stem=path.parent/f"cadd_{point['C_add_fF']:g}ff"
        for ext,key in [('.sonx','project_sha256'),('.snp','snp_sha256')]:
            if sha256(stem.with_suffix(ext).read_bytes()).hexdigest()!=point[key]:raise ValueError('Calibration archive hash changed')
        kwargs={'geometry_validator':validate} if second else {}
        check_geometry(stem.with_suffix('.sonx'),gds,meta,**kwargs)
        log=stem.with_suffix('.log').read_text()
        if 'Em errors detected: 0    Em warnings detected: 0.' not in log:raise ValueError('Calibration warnings/errors')
        f=single_sweep_resonance(read_snp(stem.with_suffix('.snp')))
        if abs(f-point['frequency_Hz'])>1:raise ValueError('Archived calibration peak mismatch')
        xy.append([point['C_add_fF']*1e-15,1/(2*np.pi*f)**2])
    slope=linregress(*np.array(xy).T).slope*1e9
    if not np.isclose(slope,r['L0_nH'],rtol=1e-12):raise ValueError('Calibration slope mismatch')
    return r,meta


def audit(name):
    output=ROOT/'output_data/taper_exploration'/name
    trial=json.loads((output/'trial.json').read_text())
    gds=Path(trial['gds']);meta=json.loads(gds.with_suffix('.json').read_text());p=from_dict(meta['second'])
    project=Path(trial['project']);snp=Path(trial['snp'])
    check_geometry(project,gds,meta,geometry_validator=validate)
    for filename,h in trial['hashes'].items():
        if sha256(Path(filename).read_bytes()).hexdigest()!=h:raise ValueError(f'Changed trial artifact: {filename}')
    status=json.loads((output/'run_status.json').read_text())
    if status['status']!='complete':raise ValueError('Mixed run incomplete')
    for path,key in [(project,'project_sha256'),(snp,'snp_sha256'),(output/'project.sonx','project_sha256'),(output/'response.snp','snp_sha256')]:
        if sha256(path.read_bytes()).hexdigest()!=status[key]:raise ValueError('Mixed run/archive identity mismatch')
    if 'Em errors detected: 0    Em warnings detected: 0.' not in (output/'em.log').read_text():raise ValueError('Mixed EM warnings/errors')
    lr,iso=calibration(trial['inductance_reference'],second=True)
    q=from_dict(iso['second'])
    keys=['taper','idc_pairs','idc_t','idc_d','idc_span','d_idc_top','x_meander','t_ground_tail',
          't_coupling_finger','l_coupling_finger','n_coupling_teeth','w_coupling_teeth',
          'w_upper_plate','l_coupling_bridge','l_cal_output','w_t1_ground_slot','l_t1_ground_slot']
    if p.t_coupling_finger:keys+=['d_t2_t1']
    for key in keys:
        if asdict(p)[key]!=asdict(q)[key]:raise ValueError(f'T2 calibration geometry mismatch: {key}')
    require_same_isolation_geometry(p,ROOT/'geometry/taper/output'/f'{Path(trial["inductance_reference"]).parent.name}.gds')
    primary=trial['network_fit']
    l1=6.014105932449249
    if p.taper1 is not None:
        r1,m1=calibration(trial['inductance_t1_reference'],second=False)
        if m1['taper']!=asdict(p.taper1):raise ValueError('T1 calibration differs from embedded metal')
        if m1.get('w_top_slot',0)!=p.w_t1_ground_slot or m1.get('l_top_slot',100)!=p.l_t1_ground_slot:
            raise ValueError('T1 ground slot calibration differs')
        l1=r1['L0_nH']
    if not np.isclose(primary.get('L_t1_nH',l1),l1,rtol=1e-12):raise ValueError('Fitter used wrong L1')
    if not np.isclose(primary['L_t2_nH'],lr['L0_nH'],rtol=1e-12):raise ValueError('Fitter used wrong L2')
    values={'L0_t2':(lr['L0_nH'],TAPER_2.L0_nH),'C0_t2':(primary['C0_t2_fF'],TAPER_2.C0_fF),
            'Cg_t2':(primary['Cg_t2_fF'],TAPER_2.Cg_input_fF),'Cg_t1':(primary['Cg_t1_fF'],TAPER_1.Cg_input_fF),
            'Cg':(primary['Cg_fF'],UNIT_CELL.Cg_input_fF),'L0_t1':(l1,TAPER_1.L0_nH),
            'C0_t1':(primary['C0_t1_fF'],TAPER_1.C0_fF),'C0_regular':(primary['C0_regular_fF'],UNIT_CELL.C0_fF)}
    errors={k:100*(v/t-1) for k,(v,t) in values.items()}
    passed=primary['success'] and not any(primary['active_bounds'][:7]) and all(abs(v)<=2 for v in errors.values())
    report={'case':name,'status':'target_pass' if passed else 'targets_not_met',
            'errors_percent':errors,'gds_project_xor_empty':True,'mesh_um':[2,2],
            'native_warnings_errors':0,'materials_and_stack_preserved':True,'ports':2,'reference_planes_um':[400,400],
            'mixed_hashes_and_archives_match':True,'calibrations_recomputed':2 if p.taper1 else 1,
            't2_calibration_metal_xor_empty':True,
            'specification':specification(),'geometry':asdict(p),'fit':primary}
    (output/'audit.json').write_text(json.dumps(report,indent=2)+'\n')
    return report
