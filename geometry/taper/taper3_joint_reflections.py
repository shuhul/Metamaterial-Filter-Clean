"""Joint reflection/transmission diagnostic on saved data. No EM/production changes."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json,time
import numpy as np
from scipy.optimize import least_squares
from geometry.taper.taper3 import OUT,digest
from geometry.taper.staged_fitting import circuit


def check(case,stage='response',*,line_z=50.):
    folder=OUT/case/stage
    fit=json.loads((folder/'fit/fit.json').read_text())
    data=np.load(folder/'fit/response.npz')
    ix=np.unique(np.linspace(0,len(data['f'])-1,1201).astype(int))
    ff=data['f'][ix];actual=data['actual'][ix]
    low=[0.,.1,.1,-100.,-100.];high=[1000.,3000.,3000.,100.,100.]
    scale=[251.,261.,153.,20.,20.]

    def model(x,f):
        cs=fit['C_fF'].copy();gs=fit['G_fF'].copy()
        cs[0]=x[0];gs[:2]=x[1:3]
        return circuit(f,fit['L_nH'],cs,gs,left_ps=x[3],right_ps=x[4],line_z=line_z)

    def residual(x,dense=False):
        f=data['f'] if dense else ff
        obs=data['actual'] if dense else actual
        delta=model(x,f)-obs
        unique=np.concatenate([delta[:,0,0],delta[:,1,0],delta[:,1,1]])
        return np.r_[unique.real,unique.imag]

    caps=[[fit['C_fF'][0],*fit['G_fF'][:2]],[251,261,153],
          [100,400,153],[450,100,153],[200,260,100],[300,200,220]]
    starts=[]
    for values in caps:
        x=np.array(values+[0.,0.]);zero=model(x,ff[:1])[0]
        for port in [0,1]:
            x[3+port]=-np.angle(actual[0,port,port]/zero[port,port])/(4*np.pi*ff[0])*1e12
        starts.append(x)
    start=time.perf_counter()
    fits=[least_squares(residual,x,bounds=(low,high),x_scale=scale,max_nfev=600,
        ftol=1e-10,xtol=1e-10,gtol=1e-10) for x in starts]
    best=min(fits,key=lambda q:q.fun@q.fun)
    q=least_squares(lambda x:residual(x,dense=True),best.x,bounds=(low,high),
        x_scale=scale,max_nfev=400,ftol=1e-10,xtol=1e-10,gtol=1e-10)
    pred=model(q.x,data['f']);delta=pred-data['actual']
    singular=np.linalg.svd(q.jac*np.array(scale),compute_uv=False)
    report=dict(case=case,stage=stage,fixture_impedance_fixed_ohm=line_z,port_reference_ohm=50.,values_C0_Ginput_Gjunction_fF=q.x[:3].tolist(),
        delays_left_right_ps=q.x[3:].tolist(),total_delay_ps=float(sum(q.x[3:])),
        S21_RMS=float(np.sqrt(np.mean(abs(delta[:,1,0])**2))),
        S11_RMS=float(np.sqrt(np.mean(abs(delta[:,0,0])**2))),
        S22_RMS=float(np.sqrt(np.mean(abs(delta[:,1,1])**2))),
        joint_RMS=float(np.sqrt(np.mean(q.fun**2)*2)),
        scaled_condition_number=float(singular[0]/singular[-1]),
        active_bounds=q.active_mask.tolist(),optimizer_success=bool(q.success),
        production_S21_only_RMS=fit['complex_RMS'],seconds=time.perf_counter()-start,
        objective='Equal-weight complex S11,S21,S22; only three new Cs and two fixture delays free; fixture impedance fixed.',
        all_L_and_older_Cs_fixed=True,source_fit_sha256=digest(folder/'fit/fit.json'),
        code_sha256=digest(__file__),production_fit_unchanged=True,
        diagnostic_only=True,limitations='Uniform-line fixtures may not describe real feeds; local condition is not confidence or uniqueness proof.')
    dest=folder/('joint_reflection_diagnostic' if line_z==50 else f'joint_reflection_diagnostic_Z{line_z:g}');dest.mkdir(exist_ok=True)
    (dest/'result.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez(dest/'response.npz',f=data['f'],actual=data['actual'],model=pred)
    print(json.dumps(report,indent=2),flush=True)
    return report
