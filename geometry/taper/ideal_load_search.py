"""Reuse ideal-component field data; verify native loading before load search.

Every source file is restored byte for byte after cache export. Results are
diagnostic directions, never physical geometry or independent C measurements.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from scipy.optimize import least_squares
from circuit_common import read_snp
from geometry.prepare_sonnet import tag,write_xml,BIN
from geometry.taper.solver_guard import run as guarded
from geometry.taper.solver_resources import solver_environment
from geometry.taper.run_bridged_stage import digest
from geometry.taper.response_cases import W
from geometry.taper.staged_fitting import circuit,LT1,CT1,LR,CR,GOUT
from geometry.taper.compact_fit import GJ
from geometry.taper.compact_stages import ROOT


def terminate(s,f,caps_fF):
    k=2*len(caps_fF);load=np.zeros((len(f),k,k),complex)
    for i,c in enumerate(caps_fF):
        v=2j*np.pi*np.asarray(f)*100*c*1e-15
        load[:,2*i,2*i]=load[:,2*i+1,2*i+1]=1/(1+v)
        load[:,2*i,2*i+1]=load[:,2*i+1,2*i]=v/(1+v)
    return s[:,:2,:2]+s[:,:2,2:]@load@np.linalg.solve(np.eye(k)[None]-s[:,2:,2:]@load,s[:,2:,:2])


def export(project,out,known_caps):
    project=Path(project);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    marker=out/'operator.json';saved=out/'operator.snp'
    if marker.exists():
        r=json.loads(marker.read_text())
        if digest(project)!=r['source_project_sha256'] or digest(project.with_suffix('.snp'))!=r['source_snp_sha256'] or digest(saved)!=r['operator_sha256']:
            raise ValueError('Cached operator/source changed')
        return W['read_network'](saved,2+2*len(known_caps))
    before=project.read_bytes();before_s=project.with_suffix('.snp').read_bytes()
    try:
        tree=ET.parse(project);components=list(tree.iter(tag('Component')))
        if len(components)!=len(known_caps):raise ValueError('Wrong ideal-component count')
        numbers=[]
        for q in components:
            numbers.extend(int(pin.get('PortNum')) for pin in q.findall(tag('Pin')))
            q.set('Type','PORTS');q.attrib.pop('Units',None)
            for e in list(q):
                if e.tag==tag('Value'):q.remove(e)
        if sorted(numbers)!=list(range(3,3+2*len(known_caps))):raise ValueError('Component port numbering is not unique')
        write_xml(tree,project);(out/'operator.sonx').write_bytes(project.read_bytes())
        with (out/'cache_export.log').open('w') as log:
            result=guarded([str(BIN/'em.exe'),'-v','-AbsNoDiscrete',str(project)],cwd=project.parent,stdout=log,env=solver_environment())
        result.check_returncode()
        if 'Total time for 0 frequencies:' not in (out/'cache_export.log').read_text():raise ValueError('Expected zero new field frequencies')
        saved.write_bytes(project.with_suffix('.snp').read_bytes());n=W['read_network'](saved,2+2*len(known_caps))
    finally:
        project.write_bytes(before);project.with_suffix('.snp').write_bytes(before_s)
    ref=read_snp(project.with_suffix('.snp'))
    if not np.array_equal(n.f,ref.f):raise ValueError('Operator grid mismatch')
    error=float(np.max(abs(terminate(n.s,n.f,known_caps)-ref.s)))
    if error>1e-7:raise ValueError(f'Native load identity failed: {error}')
    marker.write_text(json.dumps(dict(source_project_sha256=digest(project),source_snp_sha256=digest(project.with_suffix('.snp')),
        operator_sha256=digest(saved),known_capacitors_fF=known_caps,native_loading_max_S_error=error,
        zero_new_field_frequencies=True,source_restored_byte_for_byte=True),indent=2)+'\n')
    print('Native loading identity',error,flush=True);return n


def target(f):
    return circuit(f,[.7332215397274585,LT1,LR,LR],[100,CT1,CR,CR],
                   [459,56.7,GJ,GOUT],left_ps=5)


def optimize(n,out,*,upper=None,lower=None):
    out=Path(out);out.mkdir(exist_ok=True);ix=np.unique(np.linspace(0,len(n.f)-1,1501).astype(int))
    f=n.f[ix];s=n.s[ix];observed=target(f)[:,1,0];count=(s.shape[1]-2)//2
    hi=([1800]+[40,500][:count-1]) if upper is None else upper
    lo=list([0]*count if lower is None else lower)+[-20]
    hi=list(hi)+[20];scale=[500]+[10,100][:count-1]+[5]
    def pred(x):return terminate(s,f,x[:count])[:,1,0]*np.exp(-2j*np.pi*f*x[-1]*1e-12)
    def residual(x,smooth=False):
        v=pred(x);d=(1/v-1/observed)*np.maximum(abs(observed),.02)**2 if smooth else v-observed
        return np.r_[d.real,d.imag]
    seeds=[]
    for c in (0,200,420,700,1000,1400):
        for side in ((0,5,15,30) if count>=2 else (0,)):
            seed=[c]+([side] if count>=2 else [])+([20] if count>=3 else [])+[0]
            seeds.append(seed)
    smooth=[least_squares(lambda x:residual(x,True),x,bounds=(lo,hi),x_scale=scale,max_nfev=150) for x in seeds]
    fits=[least_squares(residual,q.x,bounds=(lo,hi),x_scale=scale,max_nfev=250)
          for q in sorted(smooth,key=lambda q:np.linalg.norm(residual(q.x)))[:4]]
    best=min(fits,key=lambda q:q.fun@q.fun)
    full=terminate(n.s,n.f,best.x[:count])
    phase=np.exp(-2j*np.pi*n.f*best.x[-1]*1e-12)
    full[:,1,0]*=phase;full[:,0,1]*=phase;full[:,0,0]*=phase**2
    expected=target(n.f);error=abs(full[:,1,0]-expected[:,1,0])
    result=dict(diagnostic_only=True,physical_geometry_accepted=False,
        ideal_extra_capacitors_fF=best.x[:count].tolist(),extra_fixture_delay_ps=float(best.x[-1]),
        complex_S21_RMS=float(np.sqrt(np.mean(error**2))),magnitude_RMS=float(np.sqrt(np.mean((abs(full[:,1,0])-abs(expected[:,1,0]))**2))),
        max_S21_error=float(error.max()),active_bounds=best.active_mask.tolist(),
        scope='Ideal-load algebra on cached EM fields; separate port nonperturbation control and ordinary physical EM validation are required')
    (out/'load_search.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez(out/'load_search.npz',f=n.f,target=expected,predicted=full)
    print(json.dumps(result,indent=2),flush=True);return result


if __name__=='__main__':
    p=ROOT/'sonnet_files/taper_exploration/ideal_shunt_diagnostics/plain_w500_wire52_idc30_span172_l504_c420/diagnostic_ideal_shunt.sonx'
    out=ROOT/'output_data/taper_exploration/ideal_shunt_diagnostics/plain_w500_wire52_idc30_span172_l504_c420/load_search'
    optimize(export(p,out,[420]),out)
