"""Isolated compact-T2 L calibration from one four-port EM operator.

Four ideal capacitor values are applied in Python, without four EM sweeps.
An independent direct frequency is held out of the operator interpolation.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json,runpy
import xml.etree.ElementTree as ET
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.stats import linregress
from geometry.taper.compact_inductor_case import ROOT,write_compact,validate_compact,from_dict,MeanderParams
from geometry.taper.response_cases import W
from geometry.prepare_sonnet import main as prepare,tag,write_xml
from geometry.taper.sonnet import add_capacitor

AP=runpy.run_path(str(ROOT/'investigations/taper2_fit_2026-10-22/analyze_probes.py'))

def measure(name):
    original=json.loads((ROOT/'geometry/taper/output'/f'{name}.json').read_text())
    cfg=original['external_inductor'];p=from_dict(original['second']);iso=name+'_L_cal'
    gds,out=write_compact(iso,p,MeanderParams(**cfg['meander']),cfg['x_from_cap_left'],cfg['t_bridge'],isolation=True)
    project=ROOT/'sonnet_files/taper_exploration'/iso/'generated_ncell_1.sonx'
    if not project.exists():
        prepare(1,gds_path=gds,output_dir=project.parent,geometry_validator=validate_compact)
        add_capacitor(project,json.loads(gds.with_suffix('.json').read_text()))
        tree=ET.parse(project);comp=next(tree.iter(tag('Component')))
        comp.set('Type','PORTS');comp.attrib.pop('Units',None)
        for q in list(comp):
            if q.tag==tag('Value'):comp.remove(q)
        sweeps=tree.find(tag('Sweeps'))
        for q in list(sweeps):sweeps.remove(q)
        b=ET.SubElement(sweeps,tag('Set'),On='TRUE');fs=ET.SubElement(b,tag('Frequencies'))
        for f in (4,5,5.5,6,7):
            ET.SubElement(fs,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(f),Stop=str(f),Step='1')
        for q in tree.findall(f"{tag('OutputFiles')}//{tag('File')}"):
            if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','FALSE')
        write_xml(tree,project)
    return calibrate(project,out,name)

def calibrate(project,out,name,*,auto_loads=False,desired_GHz=None,target_L_nH=.74,
              held_out_GHz=5.5, frequency_range_GHz=(4,7)):
    if not np.isfinite(target_L_nH) or target_L_nH <= 0:
        raise ValueError('Calibration target must be positive and finite')
    marker=out/'operator_completed.json'
    if not marker.exists():
        net=W['em'](project,4)
        marker.write_text(json.dumps({'project_sha256':W['digest'](project),
            'snp_sha256':W['digest'](project.with_suffix('.snp'))},indent=2)+'\n')
    else:
        saved=json.loads(marker.read_text())
        if saved['project_sha256']!=W['digest'](project) or saved['snp_sha256']!=W['digest'](project.with_suffix('.snp')):
            raise ValueError('L calibration source changed')
        net=W['read_network'](project.with_suffix('.snp'),4)
    s,z0=AP['reduce_pairs'](net.s);y=AP['s_to_y'](s,z0)
    reconstruction=AP['add_capacitors'](y,net.f,[800])
    error=float(np.max(abs(reconstruction-W['terminate_capacitor'](net.s,net.f,.8))))
    if error>1e-7:raise ValueError('Capacitor termination check failed')
    held=np.isclose(net.f,held_out_GHz*1e9,rtol=0,atol=1);train=~held
    if held.sum()!=1 or train.sum()!=4:raise ValueError('Calibration needs four training points and one held-out point')
    basis=lambda f: (np.asarray(f)[:,None]/5e9)**np.array([-1,1,3,5])
    co=np.linalg.solve(basis(net.f[train]),y[train].reshape(4,-1))
    def operator(f):return (basis(f)@co).reshape(len(f),3,3)
    validation=float(np.max(abs(operator(net.f[held])-y[held])))
    if validation>1e-5:raise ValueError(f'Interpolation validation failed: {validation} S')
    band=np.asarray(frequency_range_GHz,dtype=float)*1e9
    if band.shape!=(2,) or band[0]>=band[1] or band[0]<net.f.min()-1 or band[1]>net.f.max()+1:
        raise ValueError('Calibration search must stay within measured frequencies')
    f=np.linspace(*band,30001);yy=operator(f);rows=[]
    loads=(.8,1.,1.2,1.4)
    if auto_loads:
        desired=np.array([6.4,5.8,5.2,4.6] if desired_GHz is None else desired_GHz)*1e9
        if desired.shape!=(4,) or not np.all(np.isfinite(desired)) or np.any(np.diff(desired)>=0) or desired.min()<band[0] or desired.max()>band[1]:
            raise ValueError('Choose four decreasing calibration frequencies within the measured band')
        yn=operator(desired)
        effective=yn[:,2,2]- (yn[:,2:3,:2]@np.linalg.solve(yn[:,:2,:2]+np.eye(2)[None,:,:]/50,yn[:,:2,2:3]))[:,0,0]
        loads=-effective.imag/(2*np.pi*desired)*1e12
        if np.any(loads<=0) or not np.all(np.diff(loads)>0):raise ValueError('No positive in-band capacitor-load set')
    for cap in loads:
        observed=AP['add_capacitors'](yy,f,[cap*1000])[:,1,0];i=int(np.argmax(abs(observed)))
        if i<2 or i>len(f)-3:raise ValueError('Calibration peak outside interval')
        def amplitude(freq):return -abs(AP['add_capacitors'](operator(np.array([freq])),np.array([freq]),[cap*1000])[0,1,0])**2
        q=minimize_scalar(amplitude,bounds=(f[i-2],f[i+2]),method='bounded',options={'xatol':1})
        rows.append({'C_add_fF':cap*1000,'frequency_Hz':float(q.x)})
    caps=np.array([v['C_add_fF'] for v in rows])*1e-15;freqs=np.array([v['frequency_Hz'] for v in rows])
    if not np.all(np.diff(freqs)<0):raise ValueError('Nonmonotonic capacitor sweep')
    fit=linregress(caps,1/(2*np.pi*freqs)**2)
    result={'case':name,'L0_nH':float(fit.slope*1e9),'target_L0_nH':target_L_nH,
        'target_error_percent':float(100*(fit.slope*1e9/target_L_nH-1)),'R_squared':float(fit.rvalue**2),
        'four_loads_zero_additional_EM':True,'held_out_max_Y_error_S':validation,
        'termination_identity_max_S_error':error,'sweep':rows,'scope':'Independent isolation calibration; not an S21 target acceptance'}
    result['load_selection']='Positive loads selected for in-band resonances' if auto_loads else 'Fixed 0.8/1.0/1.2/1.4 pF loads'
    if desired_GHz is not None:result['requested_calibration_frequencies_GHz']=list(desired_GHz)
    result['adjacent_pair_L0_nH']=(np.diff(1/(2*np.pi*freqs)**2)/np.diff(caps)*1e9).tolist()
    result['load_sweep_intercept_effective_C_fF']=float(fit.intercept/fit.slope*1e15)
    result['intercept_scope']='Includes calibration feed and output loading; not the assembled resonator ground capacitance C0'
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result

if __name__=='__main__':
    import sys
    measure(sys.argv[1])
