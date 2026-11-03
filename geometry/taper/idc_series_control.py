"""Series-gap control of the full uniform-IDC layout; no shunt node probes.

Five conductor islands are identified before closing four existing wires.
T1/regular fields stay fixed. Static conductor C values are explicitly not
renamed as the historical modal LC targets.
"""
import json
import xml.etree.ElementTree as ET
import numpy as np
from scipy.optimize import least_squares
from circuit_common import read_snp
from geometry.prepare_sonnet import tag,write_xml
from geometry.taper.compact_stages import ROOT,NATIVE,GDS,OUT
from geometry.taper.plain_series import G,A,W
from geometry.taper.run_bridged_stage import digest
from geometry.taper.ideal_load_search import target

NAME='plain_w500_wire52_idc26_span174_l440'


def change(z,f,cap,x):
    # x: T2 whole-island ground C, input mutual C, optional side mutual C.
    c=cap.copy();c[1,1]+=x[0]-cap[1].sum()
    for i,j,value in [(0,1,x[1])]+([(1,2,x[2])] if len(x)>2 else []):
        delta=value+cap[i,j];c[i,i]+=delta;c[j,j]+=delta;c[i,j]-=delta;c[j,i]-=delta
    if np.linalg.eigvalsh((c+c.T)/2).min()<=0:raise ValueError('Nonpositive conductor C matrix')
    b=np.zeros((5,6))
    for port,node in enumerate([0,4,1,2,3,4]):b[node,port]=1
    return z+(b.T@(np.linalg.inv(c*1e-15)-np.linalg.inv(cap*1e-15))@b)[None]/(2j*np.pi*np.asarray(f)[:,None,None])


def control():
    source=NATIVE/NAME/'idc/generated_ncell_3.sonx'
    meta=json.loads((GDS/NAME/'idc/geometry.json').read_text())
    folder=NATIVE/NAME/'idc_series';out=OUT/NAME/'idc_series';out.mkdir(exist_ok=True)
    p=G['prepare'](source,meta,folder,first_frequency=5.2)
    marker=out/'control.json'
    if marker.exists():
        saved=json.loads(marker.read_text())
        if digest(out/'control.sonx')!=saved['project_sha256'] or digest(out/'control.snp')!=saved['snp_sha256']:
            raise ValueError('Independent series control archive changed')
        return p,out,saved
    n=W['em'](p,6);s=G['close_wires'](n.s);ref=read_snp(source.with_suffix('.snp'))
    ix=np.array([np.argmin(abs(ref.f-f)) for f in n.f])
    if max(abs(ref.f[ix]-n.f))>1:raise ValueError('Control frequency absent')
    e=abs(s-ref.s[ix])
    r=dict(project_sha256=digest(p),snp_sha256=digest(p.with_suffix('.snp')),
           frequency_GHz=(n.f/1e9).tolist(),max_S21_error=float(e[:,1,0].max()),max_all_S_error=float(e.max()),
           physical_union_XOR_empty=True,passed=bool(e[:,1,0].max()<.0005 and e.max()<.002))
    for ext in ('.sonx','.snp'):(out/('control'+ext)).write_bytes(p.with_suffix(ext).read_bytes())
    marker.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2),flush=True)
    return p,out,r


def run():
    p,out,r=control()
    if not r['passed']:raise ValueError('Series short reconstruction failed; no expanded solve')
    marker=out/'operator_completed.json'
    if not marker.exists():
        tree=ET.parse(p);fs=tree.find(f"{tag('Sweeps')}/{tag('Set')}/{tag('Frequencies')}")
        for e in list(fs):fs.remove(e)
        for f in (.5,1,2,4.65,4.9,5.15,5.2,5.4):
            ET.SubElement(fs,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(f),Stop=str(f),Step='1')
        write_xml(tree,p);W['em'](p,6)
        marker.write_text(json.dumps(dict(project_sha256=digest(p),snp_sha256=digest(p.with_suffix('.snp'))),indent=2)+'\n')
    saved=json.loads(marker.read_text())
    if saved['project_sha256']!=digest(p) or saved['snp_sha256']!=digest(p.with_suffix('.snp')):raise ValueError('Operator changed')
    return analyze(p,out)


def analyze(p,out):
    n=W['read_network'](p.with_suffix('.snp'),6);sign=np.diag([1,1,-1,-1,-1,-1]);z=A['s_to_z'](sign@n.s@sign)
    low=n.f<3e9;f=n.f[low];idx=[0,2,3,4,1]
    co=np.linalg.solve(np.c_[-1/(f/1e9),f/1e9,(f/1e9)**3],z[low][:,idx,:][:,:,idx].imag.reshape(3,-1))
    cap=np.linalg.inv(co[0].reshape(5,5))*1e6/(2*np.pi)
    held=np.isclose(n.f,5.2e9,rtol=0,atol=1);train=(n.f>4e9)&~held;powers=np.array([-1,1,3,5])
    coeff=np.linalg.solve((n.f[train,None]/5e9)**powers,z[train].reshape(4,-1))
    predict_z=lambda f:((np.asarray(f)[:,None]/5e9)**powers@coeff).reshape(len(f),6,6)
    ref=read_snp(OUT/NAME/'idc/response.snp');zz=predict_z(ref.f)
    raw=A['short_internal_z'](zz);held_error=abs(A['short_internal_z'](predict_z(n.f[held]))-A['short_internal_z'](z[held]))
    initial=np.array([cap[1].sum(),-cap[0,1]])
    def response(x):return A['short_internal_z'](change(zz,ref.f,cap,x))
    def residual(x):
        d=response(x)[:,1,0]-ref.s[:,1,0];return np.r_[d.real,d.imag]
    q=least_squares(residual,initial,bounds=([1,10],[2500,2000]),x_scale=initial,max_nfev=100)
    model=response(q.x);error=abs(model-ref.s)
    report=dict(diagnostic_only=True,physical_design_accepted=False,T1_parameters_fitted=False,
        static_conductor_order=['Input','T2','T1','Regular','End'],static_C_matrix_fF=cap.tolist(),
        independent_ground_C_fF=cap.sum(axis=1).tolist(),independent_input_mutual_C_fF=float(-cap[0,1]),
        independent_side_mutual_C_fF=float(-cap[1,2]),held_out_max_S_error=float(held_error.max()),
        raw_S21_RMS=float(np.sqrt(np.mean(abs(raw[:,1,0]-ref.s[:,1,0])**2))),
        raw_S21_max_error=float(abs(raw[:,1,0]-ref.s[:,1,0]).max()),
        two_cap_calibration_fF=q.x.tolist(),two_cap_correction_percent=(100*(q.x/initial-1)).tolist(),
        calibrated_S21_RMS=float(np.sqrt(np.mean(error[:,1,0]**2))),calibrated_S21_max_error=float(error[:,1,0].max()),
        unfitted_S11_RMS=float(np.sqrt(np.mean(error[:,0,0]**2))),
        scope='Coupled EM response and static conductor coefficients, not historical ideal-ladder target parameters')
    report['coupled_response_verified']=bool(held_error.max()<.003 and report['calibrated_S21_max_error']<.01
        and report['unfitted_S11_RMS']<.003 and max(abs(np.array(report['two_cap_correction_percent'])))<1)
    (out/'analysis.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez(out/'coupled_model.npz',coefficients=coeff,powers=powers,static_C_fF=cap,f=ref.f,actual=ref.s,calibrated=model)
    print(json.dumps(report,indent=2),flush=True);return report


def search():
    """Three static-C geometry directions, never a production parameter fit."""
    out=OUT/NAME/'idc_series';report=json.loads((out/'analysis.json').read_text())
    if not report['coupled_response_verified']:raise ValueError('Coupled response has not passed verification')
    m=np.load(out/'coupled_model.npz');cap=m['static_C_fF'];co=m['coefficients'];powers=m['powers']
    from geometry.taper.staged_fitting import circuit,LT1,CT1,LR,CR,GOUT
    from geometry.taper.compact_fit import GJ
    ell=json.loads((OUT/NAME/'idc_Lcal/results.json').read_text())['L0_nH']
    f=np.linspace(4.65e9,5.4e9,1501);z=((f[:,None]/5e9)**powers@co).reshape(len(f),6,6)
    wanted=circuit(f,[ell,LT1,LR,LR],[100,CT1,CR,CR],[459,56.7,GJ,GOUT],left_ps=5)
    obs=wanted[:,1,0]
    def predict(x):
        s=A['short_internal_z'](change(z,f,cap,x[:3]));phase=np.exp(-2j*np.pi*f*x[3]*1e-12)
        s[:,1,0]*=phase;s[:,0,1]*=phase;s[:,0,0]*=phase**2;return s
    def residual(x,smooth=False):
        v=predict(x)[:,1,0]
        d=(1/v-1/obs)*np.maximum(abs(obs),.02)**2 if smooth else v-obs
        return np.r_[d.real,d.imag]
    lo=[5,100,25,-20];hi=[2000,1200,160,20];scale=[200,450,70,5]
    seeds=[[g,400,side,0] for g in (130,300,550,900) for side in (60,75,95)]
    smooth=[least_squares(lambda x:residual(x,True),x,bounds=(lo,hi),x_scale=scale,max_nfev=150) for x in seeds]
    fits=[least_squares(residual,q.x,bounds=(lo,hi),x_scale=scale,max_nfev=250)
          for q in sorted(smooth,key=lambda q:np.linalg.norm(residual(q.x)))[:4]]
    q=min(fits,key=lambda q:q.fun@q.fun);s=predict(q.x);e=abs(s[:,1,0]-obs)
    result=dict(diagnostic_only=True,physical_design_accepted=False,T1_parameters_fitted=False,
        controls=['whole_T2_ground_C_fF','input_mutual_C_fF','side_mutual_C_fF','fixture_delay_ps'],
        x=q.x.tolist(),independent_reference_fF=[float(cap[1].sum()),float(-cap[0,1]),float(-cap[1,2])],
        complex_S21_RMS=float(np.sqrt(np.mean(e**2))),max_S21_error=float(e.max()),
        magnitude_RMS=float(np.sqrt(np.mean((abs(s[:,1,0])-abs(obs))**2))),active_bounds=q.active_mask.tolist(),
        scope='Physical static-C search directions on verified coupled fields. Final modal target fit and ordinary physical EM remain required.')
    (out/'geometry_direction_search.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez(out/'geometry_direction_search.npz',f=f,target=wanted,predicted=s)
    print(json.dumps(result,indent=2),flush=True);return result


if __name__=='__main__':run()
