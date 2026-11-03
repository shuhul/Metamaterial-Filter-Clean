"""Diagnostic combination of measured plain-wire fields and measured IDC.

Only source-island C and physical IDC mutual C are fitted to the ordinary
full-assembly S21. The measured IDC series term and all downstream EM terms
remain fixed. This is a local search guide, never final design acceptance.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import numpy as np
from scipy.optimize import least_squares
from circuit_common import read_snp
from geometry.taper.plain_series import OUT,A,change_caps
from geometry.taper.optimize_response import theory
from geometry.taper.response_cases import W


def setup(f,idc_name):
    baseline=OUT/'w500_top10_wire52_gap6/series_probes'
    m=np.load(baseline/'coupled_model.npz');co=m['coefficients'];powers=m['powers'];cap=m['static_C_fF']
    z=((np.asarray(f)[:,None]/5e9)**powers@co).reshape(len(f),6,6)
    result=json.loads((OUT/idc_name/'idc_check/results.json').read_text())
    net=W['read_network'](__import__('pathlib').Path(result['source']).with_suffix('.snp'),2)
    eye=np.eye(2)[None];y=(eye-net.s)@np.linalg.inv(eye+net.s)/50
    w=2*np.pi*net.f;c=result['low_frequency_mean_C_fF']*1e-15
    inductance=((-1/y[:,0,1]).imag/w+1/(w*w*c))[net.f>4e9].mean()
    source_c=np.mean(np.array(result['source_ground_C_fF'])[:3])*1e-15
    corrected=json.loads((baseline/'two_capacitor_fit.json').read_text())['fitted_fF']
    omega=2*np.pi*np.asarray(f)
    def predict(x,side=None,extra_m=0,mutual_scale=1):
        zz=change_caps(z,f,cap,[x[0],corrected[1] if side is None else side])
        if extra_m:
            zz[:,2,4]+=1j*omega*extra_m*1e-9;zz[:,4,2]+=1j*omega*extra_m*1e-9
        for index,value in ((3,.10416484567579731),(4,.01911955440805398),(5,.006275668516268254)):
            delta=(mutual_scale-1)*value*1e-9
            zz[:,2,index]+=1j*omega*delta;zz[:,index,2]+=1j*omega*delta
        effective=zz[:,:2,:2]-zz[:,:2,2:]@np.linalg.solve(zz[:,2:,2:],zz[:,2:,:2])
        effective[:,0,0]+=1j*omega*inductance+1/(1j*omega*x[1]*1e-15)
        yy=np.linalg.inv(effective);yy[:,0,0]+=1j*omega*source_c
        s=(eye-50*yy)@np.linalg.inv(eye+50*yy)
        phase=np.exp(-1j*omega*x[2]*1e-12);s[:,1,0]*=phase;s[:,0,1]*=phase;s[:,0,0]*=phase**2
        return s
    return predict,dict(independent_IDC_C_fF=c*1e15,independent_series_term_nH=inductance*1e9,
                        independent_source_ground_C_fF=source_c*1e15,fixed_side_C_fF=corrected[1])


def fit(name):
    p=OUT/name/'idc';n=read_snp(p/'response.snp');ix=np.unique(np.linspace(0,len(n.f)-1,1001).astype(int))
    f=n.f[ix];obs=n.s[ix];predict,known=setup(f,name)
    def res(x):
        d=predict(x)[:,1,0]-obs[:,1,0];return np.r_[d.real,d.imag]
    starts=[(c,g,t) for c in (60,130,200) for g in (known['independent_IDC_C_fF'],250,700) for t in (0,5)]
    fits=[least_squares(res,x,bounds=([1,50,-20],[600,1500,20]),x_scale=[130,460,5],max_nfev=150) for x in starts]
    best=min(fits,key=lambda q:q.fun@q.fun);model=predict(best.x)
    report=dict(diagnostic_only=True,accepted=False,free_electrical_parameters=['source_ground_C_surrogate_fF','input_C_surrogate_fF'],
        x=best.x.tolist(),fixed_independently_measured=known,active_bounds=best.active_mask.tolist(),
        S21_complex_RMS=float(np.sqrt(np.mean(abs(model[:,1,0]-obs[:,1,0])**2))),
        S21_max_error=float(np.max(abs(model[:,1,0]-obs[:,1,0]))),
        unfitted_S11_RMS=float(np.sqrt(np.mean(abs(model[:,0,0]-obs[:,0,0])**2))),
        scope='Bare coupled fields plus measured IDC pi equivalent. Fitted coefficients are surrogate controls, not certified component values; their difference from independent IDC C flags imperfect geometry transfer.')
    dest=p/'cascade_diagnostic';dest.mkdir(exist_ok=True);(dest/'fit.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)
    return report


def optimize(name,*,vary_mutual=False,measured_reference=False):
    f=np.linspace(4.65e9,5.4e9,1501);target=theory(f)
    if measured_reference:
        from geometry.taper.staged_fitting import circuit,LT1,CT1,LR,CR,GOUT
        from geometry.taper.compact_fit import GJ
        ell=json.loads((OUT/name/'idc_Lcal/results.json').read_text())['L0_nH']
        target=circuit(f,[ell,LT1,LR,LR],[100,CT1,CR,CR],[459,56.63965936392767,GJ,GOUT],left_ps=5)
    predict,known=setup(f,name)
    previous=json.loads((OUT/name/'idc/cascade_diagnostic/fit.json').read_text())['x']
    seed=np.r_[previous,known['fixed_side_C_fF'],1.] if vary_mutual else np.r_[previous,known['fixed_side_C_fF']]
    lo=[1,50,-20,20]+([0] if vary_mutual else []);hi=[600,1500,20,160]+([1.5] if vary_mutual else [])
    def model(x):return predict(x[:3],side=x[3],mutual_scale=x[4] if vary_mutual else 1)
    def res(x,smooth=False):
        a=model(x)[:,1,0];b=target[:,1,0]
        d=(1/a-1/b)*np.maximum(abs(b),.02)**2 if smooth else a-b
        return np.r_[d.real,d.imag]
    rng=np.random.default_rng(42)
    seeds=[seed]
    for i in range(5):
        v=seed.copy();v[:2]*=rng.uniform(.6,1.5,2);v[3]*=rng.uniform(.75,1.4)
        if vary_mutual:v[4]=rng.uniform(.2,1.2)
        seeds.append(np.clip(v,np.array(lo)+1e-7,np.array(hi)-1e-7))
    scale=[130,460,5,70]+([1] if vary_mutual else [])
    smooth=[least_squares(lambda x:res(x,True),x,bounds=(lo,hi),x_scale=scale,max_nfev=150) for x in seeds]
    fits=[least_squares(res,q.x,bounds=(lo,hi),x_scale=scale,max_nfev=250) for q in sorted(smooth,key=lambda q:np.linalg.norm(res(q.x)))[:3]]
    best=min(fits,key=lambda q:q.fun@q.fun);pred=model(best.x)
    report=dict(diagnostic_only=True,physical_design_accepted=False,vary_mutual=vary_mutual,measured_downstream_reference=measured_reference,
        controls=['source_C_surrogate_fF','input_C_surrogate_fF','extra_delay_ps','side_C_surrogate_fF']+(['wire_mutual_scale'] if vary_mutual else []),
        x=best.x.tolist(),active_bounds=best.active_mask.tolist(),
        complex_RMS=float(np.sqrt(np.mean(abs(pred[:,1,0]-target[:,1,0])**2))),
        magnitude_RMS=float(np.sqrt(np.mean((abs(pred[:,1,0])-abs(target[:,1,0]))**2))),
        scope='Local coupled-model direction toward nominal S21, not a fit or acceptance of a physical geometry; validate with ordinary Sonnet EM')
    dest=OUT/name/'idc/cascade_diagnostic';suffix='with_mutual' if vary_mutual else 'capacitors_only'
    if measured_reference:suffix+='_measured_reference'
    (dest/f'optimize_{suffix}.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez(dest/f'optimize_{suffix}.npz',f=f,target=target,predicted=pred)
    print(json.dumps(report,indent=2),flush=True);return report


if __name__=='__main__':fit('plain_w500_wire52_idc30_span172_l504')
