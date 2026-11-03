"""Promote only an audited, target-passing taper-2 trial; retain native cache."""
from pathlib import Path
from dataclasses import asdict
import json
import shutil
import os
import numpy as np
from scipy.signal import find_peaks
from circuit_common import read_snp
from geometry.taper.second_audit import audit,ROOT
from geometry.taper.second_experiment import model_response,plot_final,L1


def finalize(name):
    output=ROOT/'output_data/taper_exploration'/name
    trial=json.loads((output/'trial.json').read_text())
    if not trial['passes_2_percent']:raise ValueError('Trial did not meet targets; no promotion')
    report=audit(name)
    if report['status']!='target_pass':raise ValueError('Independent audit did not confirm target pass')
    fit=report['fit'];net=read_snp(output/'response.snp')
    model=model_response(net.f,fit['x'],fit['L_t2_nH'],l1=fit.get('L_t1_nH',L1))
    peaks,_=find_peaks(abs(net.s[:,1,0]),prominence=.02)
    peaks=[i for i in peaks if net.f[i]<6e9]
    if len(peaks)!=2:raise ValueError('Unexpected lower passband modes')
    mask=(net.f>=net.f[peaks[0]]-.35e9)&(net.f<=net.f[peaks[-1]]+.35e9)
    error=abs(model[mask,1,0]-net.s[mask,1,0])
    quality={'passband_complex_rmse':float(np.sqrt(np.mean(error**2))),
             'passband_max_complex_error':float(error.max()),
             's11_complex_rmse':fit['s11_complex_rmse']}
    if quality['passband_complex_rmse']>.015 or quality['passband_max_complex_error']>.06 or fit['s11_complex_rmse']>.03:
        raise ValueError(f'Targets pass but fit quality needs review: {quality}')
    source=Path(trial['project']).parent.resolve()
    dest=(ROOT/'sonnet_files/taper_exploration/final_designs'/name).resolve()
    if not source.is_relative_to(ROOT) or not dest.is_relative_to(ROOT):raise ValueError('Project outside workspace')
    if dest.exists():raise ValueError('Final destination already exists')
    dest.parent.mkdir(parents=True,exist_ok=True)
    # Finish the reviewable plot before moving the complete project directory.
    plot_final(output/'response.snp',fit,output/'s21_fit_linear.png')
    shutil.move(str(source),str(dest))
    # Rewrite live JSON paths, including hash-map keys, without touching project bytes.
    records=list((ROOT/'output_data/taper_exploration').rglob('*.json'))+list(dest.glob('*.json'))
    for p in records:
        s=p.read_text(encoding='utf-8')
        updated=s.replace(str(source).replace('\\','\\\\'),str(dest).replace('\\','\\\\'))
        if updated!=s:p.write_text(updated,encoding='utf-8')
    trial=json.loads((output/'trial.json').read_text())
    shutil.copyfile(output/'s21_fit_linear.png',dest/'s21_fit_linear.png')
    shutil.copyfile(output/'s21_fit_linear.png',output.parent/'second_taper_s21_fit_linear.png')
    shutil.copyfile(trial['gds'],dest/'geometry.gds')
    shutil.copyfile(Path(trial['gds']).with_suffix('.json'),dest/'geometry.json')
    results=dest/'results';results.mkdir()
    for filename in ('trial.json','fit.json','fit_network.json','audit.json','run_status.json'):
        shutil.copyfile(output/filename,results/filename)
    selected={'stage':'second_taper','status':'selected_screening_pass','case':name,
              'geometry':trial['parameters'],'specification':report['specification'],
              'fit_method':trial['primary_fit'],'fit':fit,'errors_percent':report['errors_percent'],
              'fit_quality':quality,'project':str(dest/'generated_ncell_3.sonx'),
              'snp':str(dest/'generated_ncell_3.snp'),'plot':str(dest/'s21_fit_linear.png'),
              'inductance_t2':trial['inductance_reference'],'inductance_t1':trial.get('inductance_t1_reference'),
              'trial':str(output/'trial.json'),'audit':str(output/'audit.json'),
              'limitations':['Fixed 2 um mesh; not a mesh-convergence study.',
                             'Equivalent-circuit parameters depend on the fixture and model.',
                             'Isolated inductance calibrations transfer into the mixed environment.',
                             'The IDC has distributed parasitics; capacitances are fitted effective values.']}
    (ROOT/'geometry/taper/selected_second_design.json').write_text(json.dumps(selected,indent=2)+'\n')
    def link(path):return os.path.relpath(path,dest).replace('\\','/')
    rows='\n'.join(f'| {key} | {value:+.3f}% |' for key,value in report['errors_percent'].items())
    text=f'''# Selected second-taper assembly

Open [generated_ncell_3.sonx](generated_ncell_3.sonx).
Share [s21_fit_linear.png](s21_fit_linear.png).

The entire native project directory and EM cache were moved here after an
independent geometry/settings/response/calibration audit. All recorded target
errors are within 2% for the selected circuit extraction.

| Quantity | Signed target error |
|---|---:|
{rows}

The fit uses complex S21 and S11 with weights 3:1 and matched effective port
delays. All capacitances are free in the final fit. L values come from the
independent added-capacitance regressions. The annotated plot displays linear
|S21|; it is not a magnitude-only fit. Effective delay rotates phase only.

## Geometry

```json
{json.dumps(trial['parameters'],indent=2)}
```

If `taper1` is non-null above, this assembly uses an independently calibrated
embedded first-taper variant. It does not overwrite the earlier accepted
two-resonator first-taper project. The driven coupling finger belongs to the
taper-2 signal node; it is not a ground trace. The source GDS and metadata
are included as `geometry.gds` and `geometry.json`.

## Records and limits

- [Full trial record]({link(output/'trial.json')})
- [Independent audit]({link(output/'audit.json')})
- [Taper-2 inductance]({link(Path(trial['inductance_reference']))})
- [Complete exploration notes]({link(ROOT/'sonnet_files/taper_exploration/TAPER_2_README.md')})

Passband complex S21 RMSE: {quality['passband_complex_rmse']:.6g}.
Maximum passband complex S21 error: {quality['passband_max_complex_error']:.6g}.
Input-reflection complex RMSE: {quality['s11_complex_rmse']:.6g}.

This is an EM/equivalent-circuit screening pass on the requested 2 x 2 um grid,
not fabrication validation or a finer-mesh convergence result. The inherited
sheet material and stack were held fixed. IDC parasitics and isolated-to-mixed
L transfer remain modeling assumptions; see the full notes and alternative fits.
'''
    (dest/'README.md').write_text(text,encoding='utf-8')
    print(json.dumps(selected,indent=2))
    return selected
