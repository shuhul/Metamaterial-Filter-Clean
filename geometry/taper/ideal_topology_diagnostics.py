"""Offline topology/constant ablation; NEVER updates production T1 parameters.

The selected T1 constants were obtained with a two-resonator/three-capacitor
model. The staged model adds the feed-connected output-end LC and replaces
regular/junction constants. Test those changes separately on the same data.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import numpy as np
from scipy.optimize import least_squares
from geometry.taper.ideal_inductor_audit import OUT
from geometry.taper.staged_fitting import circuit,LT1,CT1,LR,CR,GOUT
from geometry.taper.locked_taper2 import FIXED
from geometry.taper.paired_coupling import abcd


def prediction(f,l,x,variant,free_t1=False):
    c,g,tau=x[:3]
    ct=x[3] if free_t1 else CT1
    selected=variant.startswith('selected')
    cr=FIXED['C_regular_fF'] if selected else CR
    gj=FIXED['Cg_junction_fF'] if selected else 8.40642914296533
    go=FIXED['Cg_output_fF'] if selected else GOUT
    ls=[LT1,LR,LR];cs=[ct,cr,cr];gs=[None,gj,go]
    if variant=='selected_original_topology':ls=ls[:2];cs=cs[:2];gs=gs[:2]
    n=abcd(circuit(f,ls,cs,gs))
    if variant=='selected_original_topology':
        # Original circuit ends after the output coupling capacitor; no end LC.
        z=1/(2j*np.pi*f*go*1e-15)
        n[:,:,1]+=n[:,:,0]*z[:,None]
    w=2*np.pi*f;z=1/(1j*w*g*1e-15)
    a=n[:,0,0]+z*n[:,1,0];b=n[:,0,1]+z*n[:,1,1]
    y=1j*w*c*1e-15+(0 if l is None else 1/(1j*w*l*1e-9))
    den=(a+b/50)*(1+50*y)+50*n[:,1,0]+n[:,1,1]
    return 2/den*np.exp(-1j*w*tau*1e-12)


def fit(l,variant,free_t1=False,ct_fixed=None,data=None):
    label='open' if l is None else f'L{l:g}'
    a=np.load(OUT/'same_geometry'/label/'fit/response.npz') if data is None else data
    ix=np.unique(np.linspace(0,len(a['f'])-1,1501).astype(int));f=a['f'][ix];obs=a['actual'][ix,1,0]
    lower=[0,5,-20]+([CT1*.95] if free_t1 else [])
    upper=[2000,180,30]+([CT1*1.05] if free_t1 else [])
    def predict(x,frequencies=f):
        return prediction(frequencies,l,list(x[:3])+([ct_fixed] if ct_fixed is not None else list(x[3:])),
                          variant,free_t1 or ct_fixed is not None)
    def residual(x,smooth=False):
        pred=predict(x)
        d=(1/pred-1/obs)*np.maximum(abs(obs),.03)**2 if smooth else pred-obs
        return np.r_[d.real,d.imag]
    seeds=[]
    for c in (0,100,300,800):
        for g in (54,56.7,60):
            x=[c,g,5]+([CT1] if free_t1 else [])
            q=least_squares(lambda v:residual(v,True),x,bounds=(lower,upper),x_scale='jac',max_nfev=100)
            seeds.append(q.x)
    candidates=[least_squares(residual,x,bounds=(lower,upper),x_scale='jac',max_nfev=200,
                             ftol=1e-11,xtol=1e-11,gtol=1e-11) for x in sorted(seeds,key=lambda v:np.linalg.norm(residual(v)))[:3]]
    best=min(candidates,key=lambda q:q.fun@q.fun)
    full=predict(best.x,a['f']);error=abs(full-a['actual'][:,1,0])
    r=dict(L_nH=l,variant=variant,diagnostic_T1_C_released=free_t1,production_parameters_changed=False,
           x=best.x.tolist(),fixed_diagnostic_CT1_fF=ct_fixed,active_bounds=best.active_mask.tolist(),
           S21_complex_RMS=float(np.sqrt(np.mean(error**2))),S21_complex_max=float(error.max()))
    return r,full,a


def run():
    rows=[];models={}
    for l in (None,5.,.74):
        for variant,release in [('current_staged',False),('selected_with_end_LC',False),
                                ('selected_original_topology',False),('current_staged',True)]:
            r,model,a=fit(l,variant,release)
            rows.append(r);print(json.dumps(r),flush=True)
            models[f'{l}_{variant}_{release}']=model
    # Held-out load prediction: refit T1 C on OPEN only, freeze for loaded fits.
    base=next(r for r in rows if r['L_nH'] is None and r['diagnostic_T1_C_released'])
    for l in (5.,2.,.74):
        r,model,a=fit(l,'current_staged',ct_fixed=base['x'][3])
        r['diagnostic_T1_C_source']='Fitted on no-inductor case only; held fixed on loaded cases'
        rows.append(r);print(json.dumps(r),flush=True)
        models[f'{l}_open_calibrated']=model
    report=dict(diagnostic_only=True,production_parameters_changed=False,rows=rows,
                original_topology='Two resonators, three coupling caps; no output-end LC',
                staged_topology='Adds feed-connected output-end LC; changes regular/junction constants')
    (OUT/'topology_diagnostics.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez(OUT/'topology_models.npz',**models)
    return report


if __name__=='__main__':run()
