"""Nonperturbing series-port checks for the plain, directly fed T2 fixture.

T1 metal is untouched. Shorting the added zero-width ports must reproduce the
ordinary two-port simulation before any wider operator sweep is authorized.
This coupled EM model is a diagnostic; scalar circuit targets are not renamed.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json,runpy,shutil
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from scipy.optimize import least_squares
from circuit_common import read_snp
from geometry.taper.compact_stages import ROOT,GDS,NATIVE,OUT
from geometry.taper.run_bridged_stage import digest
from geometry.prepare_sonnet import tag,write_xml

G=runpy.run_path(str(ROOT/'investigations/taper2_fit_2026-10-22/gap_probes.py'))
A=runpy.run_path(str(ROOT/'investigations/taper2_fit_2026-10-22/analyze_series.py'))
W=G['W']


def control(name):
    source=NATIVE/name/'direct/generated_ncell_3.sonx'
    meta=json.loads((GDS/name/'direct/geometry.json').read_text())
    if any(meta['second'].get(k,0) for k in ('w_upper_plate','t_inner_ground','l_ground_rise','x_meander')):
        raise ValueError('This control is restricted to the plain centred-wire design')
    dest=NATIVE/name/'series_probes';out=OUT/name/'series_probes';out.mkdir(parents=True,exist_ok=True)
    p=G['prepare'](source,meta,dest,first_frequency=5.2)
    marker=out/'control.json'
    if marker.exists():
        report=json.loads(marker.read_text())
        if digest(dest/'control.sonx')!=report['project_sha256'] or digest(dest/'control.snp')!=report['snp_sha256']:
            raise ValueError('Control archive changed')
        return p,out,report
    n=W['em'](p,6);shorted=G['close_wires'](n.s)
    ref=read_snp(OUT/name/'direct/response.snp')
    ix=np.array([np.argmin(abs(ref.f-f)) for f in n.f])
    if np.max(abs(ref.f[ix]-n.f))>1:raise ValueError('Control frequency absent from reference')
    error=abs(shorted-ref.s[ix])
    for ext in ('.sonx','.snp'):shutil.copyfile(p.with_suffix(ext),dest/('control'+ext))
    shutil.copyfile(dest/'em.log',out/'control_em.log')
    report=dict(project_sha256=digest(p),snp_sha256=digest(p.with_suffix('.snp')),
        frequency_GHz=(n.f/1e9).tolist(),max_S21_error=float(error[:,1,0].max()),max_all_S_error=float(error.max()),
        physical_union_XOR_empty=True,passed=bool(error[:,1,0].max()<5e-4 and error.max()<.002))
    marker.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2),flush=True)
    return p,out,report


def run(name):
    p,out,report=control(name)
    if not report['passed']:raise ValueError('Series-port short reconstruction failed; no wider solve')
    marker=out/'completed.json'
    if marker.exists():
        r=json.loads(marker.read_text())
        if digest(p)!=r['project_sha256'] or digest(p.with_suffix('.snp'))!=r['snp_sha256']:
            raise ValueError('Series operator changed')
        return analyze(name)
    tree=ET.parse(p);fs=tree.find(f"{tag('Sweeps')}/{tag('Set')}/{tag('Frequencies')}")
    for q in list(fs):fs.remove(q)
    for f in (.5,1,2,4.65,4.9,5.15,5.2,5.4):
        ET.SubElement(fs,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(f),Stop=str(f),Step='1')
    write_xml(tree,p);n=W['em'](p,6)
    marker.write_text(json.dumps(dict(project_sha256=digest(p),snp_sha256=digest(p.with_suffix('.snp')),
        frequencies_GHz=(n.f/1e9).tolist(),held_out_GHz=5.2),indent=2)+'\n')
    return analyze(name)


def analyze(name):
    p=NATIVE/name/'series_probes/series_probe.sonx';out=OUT/name/'series_probes'
    n=W['read_network'](p.with_suffix('.snp'),6)
    polarity=np.diag([1,1,-1,-1,-1,-1]);s=polarity@n.s@polarity;z=A['s_to_z'](s)
    low=n.f<3e9;idx=[0,3,4,1]
    f=n.f[low];basis=np.c_[-1/(f/1e9),f/1e9,(f/1e9)**3]
    coef=np.linalg.solve(basis,z[low][:,idx,:][:,:,idx].imag.reshape(3,-1))
    cap=np.linalg.inv(coef[0].reshape(4,4))*1e6/(2*np.pi)
    low_full_coef=np.linalg.solve(basis,z[low].imag.reshape(3,-1))
    loop=np.array([-1.,0,1,0,0,0])
    loop_low_L=low_full_coef[1].reshape(6,6)/(2*np.pi)
    held=np.isclose(n.f,5.2e9,atol=1,rtol=0);train=(n.f>4e9)&~held
    powers=np.arange(-1,2*train.sum()-2,2)
    co=np.linalg.solve((n.f[train,None]/5e9)**powers,z[train].reshape(sum(train),-1))
    def predict_z(f):return ((np.asarray(f)[:,None]/5e9)**powers@co).reshape(len(f),6,6)
    ref=read_snp(OUT/name/'direct/response.snp');prediction=A['short_internal_z'](predict_z(ref.f))
    held_error=abs(A['short_internal_z'](predict_z(n.f[held]))-A['short_internal_z'](z[held]))
    error=abs(prediction[:,1,0]-ref.s[:,1,0])
    report=dict(diagnostic_only=True,accepted_T2_design=False,source=str(p),
        conductor_order=['Directly fed T2','T1','Regular','Output end'],
        ground_C_fF=cap.sum(axis=1).tolist(),mutual_T2_T1_C_fF=float(-cap[0,1]),
        static_C_matrix_fF=cap.tolist(),C_matrix_positive_definite=bool(np.linalg.eigvalsh((cap+cap.T)/2).min()>0),
        held_out_max_all_S_error=float(held_error.max()),held_out_S21_error=float(held_error[:,1,0].max()),
        dense_S21_complex_RMS=float(np.sqrt(np.mean(error**2))),dense_S21_max_error=float(error.max()),
        model_verified=bool(held_error[:,1,0].max()<.003 and error.max()<.01),
        interpretation='Whole-conductor static C and independently calibrated coupled response; not automatically the modal LC targets')
    report['defined_input_to_T2_wire_loop_L_nH']=float(loop@loop_low_L@loop)
    report['defined_loop_mutual_coefficients_nH']=(loop@loop_low_L).tolist()
    report['loop_coefficient_scope']='Laurent impedance coefficients of a specified closed current-loop excitation, not unique ideal-ladder element assignments'
    np.savez(out/'coupled_model.npz',coefficients=co,powers=powers,static_C_fF=cap,
             f=ref.f,actual=ref.s,predicted=prediction)
    (out/'analysis.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)
    if report['model_verified'] and report['C_matrix_positive_definite']:
        fit_caps(name,ref,predict_z,cap)
    return report


def change_caps(z,f,reference,x):
    """Two physical capacitance controls; all other EM coefficients are fixed."""
    c=reference.copy();c[0,0]+=x[0]-reference[0].sum()
    delta=x[1]+reference[0,1]
    c[0,0]+=delta;c[1,1]+=delta;c[0,1]-=delta;c[1,0]-=delta
    if np.linalg.eigvalsh((c+c.T)/2).min()<=0:raise ValueError('Nonpositive capacitance matrix')
    b=np.zeros((4,6))
    for port,node in enumerate([0,3,0,1,2,3]):b[node,port]=1
    dz=b.T@(np.linalg.inv(c*1e-15)-np.linalg.inv(reference*1e-15))@b
    return z+dz[None]/(2j*np.pi*np.asarray(f)[:,None,None])


def fit_caps(name,ref,predict_z,cap):
    f=ref.f;zz=predict_z(f);initial=np.array([cap[0].sum(),-cap[0,1]])
    def response(x):return A['short_internal_z'](change_caps(zz,f,cap,x))
    def residual(x):
        d=response(x)[:,1,0]-ref.s[:,1,0];return np.r_[d.real,d.imag]
    fits=[least_squares(residual,initial*q,bounds=([.1,.1],[2000,2000]),x_scale=initial,max_nfev=150)
          for q in ([1,1],[.98,1],[1.02,1],[1,.98],[1,1.02])]
    best=min(fits,key=lambda q:q.fun@q.fun);pred=response(best.x)
    result=dict(free_electrical_parameters=['whole_T2_source_island_ground_C_fF','whole_T2_T1_mutual_C_fF'],
        independent_reference_fF=initial.tolist(),fitted_fF=best.x.tolist(),
        calibrated_to_fitted_difference_percent=(100*(best.x/initial-1)).tolist(),
        complex_S21_RMS=float(np.sqrt(np.mean(abs(pred[:,1,0]-ref.s[:,1,0])**2))),
        max_S21_error=float(np.max(abs(pred[:,1,0]-ref.s[:,1,0]))),
        unfitted_S11_RMS=float(np.sqrt(np.mean(abs(pred[:,0,0]-ref.s[:,0,0])**2))),
        T1_parameters_fitted=False,accepted_T2_design=False,
        scope='Independently calibrated coupled model; static conductor capacitances retain their explicit definition, not the historical modal target definition')
    out=OUT/name/'series_probes';(out/'two_capacitor_fit.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez(out/'two_capacitor_fit.npz',f=f,actual=ref.s,predicted=pred)
    print(json.dumps(result,indent=2),flush=True)
    return result


if __name__=='__main__':run('w500_top10_wire52_gap6')
