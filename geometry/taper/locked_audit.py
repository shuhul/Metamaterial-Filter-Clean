"""Audit/promote only T2 trials made with verified T1 and two fitted parameters."""
from pathlib import Path
from dataclasses import asdict
from hashlib import sha256
import json
import os
import shutil
import numpy as np
from scipy.signal import find_peaks
from circuit_common import read_snp
from geometry.taper.locked_taper2 import ROOT,SOURCE,FIXED,locked_validate,require_native_locked,response,plot
from geometry.taper.second_layout import from_dict,require_same_isolation_geometry
from geometry.taper.audit import check_geometry
from geometry.taper.second_audit import calibration
from geometry.specifications import TAPER_2


def audit(name):
    out=ROOT/'output_data/taper_exploration'/name
    trial=json.loads((out/'trial_locked.json').read_text());fit=trial['fit']
    for file,digest in trial['hashes'].items():
        if sha256(Path(file).read_bytes()).hexdigest()!=digest:raise ValueError('Changed artifact: '+file)
    if fit['free_parameters']!=['C0_t2_fF','Cg_t2_fF'] or len(fit['x'])!=2:
        raise ValueError('Fit is not limited to two capacitances')
    if fit['fixed_model']!=FIXED or fit['fixed_reference_sha256']!=sha256(SOURCE.read_bytes()).hexdigest():
        raise ValueError('Fixed model/reference changed')
    gds=Path(trial['gds']);meta=json.loads(gds.with_suffix('.json').read_text());p=from_dict(meta['second'])
    check_geometry(trial['project'],gds,meta,geometry_validator=locked_validate)
    require_native_locked(trial['project'],meta)
    status=json.loads((out/'run_status.json').read_text())
    if status['status']!='complete':raise ValueError('Mixed simulation incomplete')
    for file,key in [(trial['project'],'project_sha256'),(trial['snp'],'snp_sha256'),
                     (out/'project.sonx','project_sha256'),(out/'response.snp','snp_sha256')]:
        if sha256(Path(file).read_bytes()).hexdigest()!=status[key]:raise ValueError('Mixed archive identity mismatch')
    log=(out/'em.log').read_text()
    if 'Em errors detected: 0    Em warnings detected: 0.' not in log:raise ValueError('EM warning/error')
    if 'using 4 thread(s)' not in log:raise ValueError('Four-thread policy not recorded')
    inductance,cal_meta=calibration(trial['inductance_reference'],second=True)
    ref_gds=ROOT/'geometry/taper/output'/f'{Path(trial["inductance_reference"]).parent.name}.gds'
    require_same_isolation_geometry(p,ref_gds)
    if not np.isclose(inductance['L0_nH'],fit['L_t2_nH'],rtol=1e-12):raise ValueError('L mismatch')
    net=read_snp(out/'response.snp');pred=response(net.f,fit['x'],fit['L_t2_nH'])
    error=abs(pred[:,1,0]-net.s[:,1,0]);rmse=float(np.sqrt(np.mean(error**2)))
    if not np.isclose(rmse,fit['s21_complex_rmse'],rtol=1e-10,atol=1e-12):raise ValueError('Fit response mismatch')
    peaks,_=find_peaks(abs(net.s[:,1,0]),prominence=.02)
    peaks=[i for i in peaks if net.f[i]<6e9]
    band=(net.f>=net.f[peaks[0]]-.35e9)&(net.f<=net.f[peaks[-1]]+.35e9) if len(peaks)==2 else np.ones(len(net.f),dtype=bool)
    quality={'s21_complex_rmse':rmse,'lower_passband_peak_count':len(peaks),
             'passband_complex_rmse':float(np.sqrt(np.mean(error[band]**2))),
             'passband_max_complex_error':float(error[band].max())}
    curve_pass=len(peaks)==2 and quality['passband_complex_rmse']<=.015 and quality['passband_max_complex_error']<=.06
    errors={'L0_t2':100*(fit['L_t2_nH']/TAPER_2.L0_nH-1),
            'C0_t2':100*(fit['C0_t2_fF']/TAPER_2.C0_fF-1),
            'Cg_t2':100*(fit['Cg_t2_fF']/TAPER_2.Cg_input_fF-1)}
    target_pass=fit['success'] and not any(fit['active_bounds']) and all(abs(v)<=2 for v in errors.values())
    report={'case':name,'status':'pass' if target_pass and curve_pass else 'not_accepted',
            'target_pass':bool(target_pass),'curve_pass':bool(curve_pass),'errors_percent':errors,
            'quality':quality,'T1_and_right_hand_metal_xor_empty':True,
            'native_Sonnet_T1_and_right_hand_metal_xor_empty':True,
            'source_imported_metal_and_template_settings_checked':True,'mesh_um':[2,2],
            'ports':2,'reference_planes_um':[400,400],'native_warnings_errors':0,
            'fixed_model_reference_sha256':fit['fixed_reference_sha256'],
            'inductance_regression_recomputed_from_archived_peaks':True,
            'inductance_calibration_metal_xor_empty':True,
            'free_parameters':fit['free_parameters'],'fixed_model':FIXED,'geometry':asdict(p)}
    (out/'audit_locked.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)
    return report


def finalize(name):
    report=audit(name)
    if report['status']!='pass':raise ValueError('Trial not accepted; target/curve audit failed')
    out=ROOT/'output_data/taper_exploration'/name
    record=json.loads((out/'trial_locked.json').read_text());source=Path(record['project']).parent.resolve()
    dest=(ROOT/'sonnet_files/taper_exploration/final_designs'/name).resolve()
    if not source.is_relative_to(ROOT) or not dest.is_relative_to(ROOT):raise ValueError('Move outside workspace')
    if dest.exists():raise ValueError('Final folder already exists')
    plot(out/'response.snp',record['fit'],out/'s21_fit_linear.png')
    dest.parent.mkdir(parents=True,exist_ok=True);shutil.move(str(source),str(dest))
    for path in list(out.glob('*.json'))+list(dest.glob('*.json')):
        text=path.read_text();updated=text.replace(str(source).replace('\\','\\\\'),str(dest).replace('\\','\\\\'))
        if updated!=text:path.write_text(updated)
    record=json.loads((out/'trial_locked.json').read_text())
    for filename in ('s21_fit_linear.png','geometry_review.png','geometry_check.json','fit_locked.json','trial_locked.json','audit_locked.json','run_status.json'):
        shutil.copyfile(out/filename,dest/filename)
    for ext in ('.gds','.json'):shutil.copyfile(Path(record['gds']).with_suffix(ext),dest/f'geometry{ext}')
    shutil.copyfile(out/'s21_fit_linear.png',out.parent/'second_taper_s21_fit_linear.png')
    selected={'stage':'second_taper','status':'selected_screening_pass','mode':'verified_T1_locked_two_capacitances',
              'case':name,'project':str(dest/'generated_ncell_3.sonx'),'snp':str(dest/'generated_ncell_3.snp'),
              'plot':str(dest/'s21_fit_linear.png'),'geometry_review':str(dest/'geometry_review.png'),
              'geometry':record['parameters'],'fit':record['fit'],'audit':str(out/'audit_locked.json'),
              'inductance_reference':record['inductance_reference'],'first_taper_reference':str(SOURCE),
              'limitations':['Fixed 2 um grid; not a mesh-convergence study.',
                             'T1/regular parameters and reference delays held fixed as requested.']}
    (ROOT/'geometry/taper/selected_second_design.json').write_text(json.dumps(selected,indent=2)+'\n')
    rows='\n'.join(f'| {k} | {v:+.3f}% |' for k,v in report['errors_percent'].items())
    link=os.path.relpath(ROOT/'sonnet_files/taper_exploration/LOCKED_TAPER2_README.md',dest).replace('\\','/')
    (dest/'README.md').write_text(f'''# Selected Taper 2 with verified Taper 1 unchanged

Open [generated_ncell_3.sonx](generated_ncell_3.sonx).
Share [s21_fit_linear.png](s21_fit_linear.png).
Inspect [geometry_review.png](geometry_review.png): blue is T2, orange is
verified T1, gray is ground. All metal is on GDS layer 1/datatype 0.

Only C0_t2 and input Cg_t2 were fitted, using complex S21 real/imaginary parts.
L2 was independently calibrated. All T1/regular values, couplings and phase
delays were fixed to the saved first-taper model. The complete metal region
from T1 to the right box wall matches the verified geometry by polygon XOR.

| Target | Signed error |
|---|---:|
{rows}

The project was promoted only after archived-response identity, native geometry,
materials, mesh, ports, inductance regression, target and curve-quality checks.
The full native project folder and EM cache are retained here.

[Workflow and assumptions]({link}).
See `fit_locked.json`, `audit_locked.json`, and `trial_locked.json` for machine
records. This is a fixed-grid equivalent-circuit screening result.
''')
    return selected
