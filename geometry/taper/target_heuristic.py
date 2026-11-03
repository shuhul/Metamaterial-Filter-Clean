"""Diagnostic all-capacitor fit; never used as target-response acceptance."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from circuit_common import read_snp
from geometry.taper.staged_fitting import LT1, LR
from geometry.taper.corrected_iteration import REGULAR_C, REGULAR_G
from geometry.taper.run_bridged_stage import digest


def fast_network(f, ls, cs, gs, left_ps=0, right_ps=0):
    """Equivalent ABCD cascade without allocating batches of 2x2 matrices."""
    w=2*np.pi*np.asarray(f)
    a=np.ones(len(w),complex);b=np.zeros(len(w),complex);c=b.copy();d=a.copy()
    for l,cap,gap in zip(ls,cs,gs):
        if gap is not None:
            z=1/(1j*w*gap*1e-15);b=b+a*z;d=d+c*z
        y=1j*w*cap*1e-15+1/(1j*w*l*1e-9);a=a+b*y;c=c+d*y
    den=a+b/50+50*c+d
    s=np.empty((len(w),2,2),complex)
    s[:,0,0]=(a+b/50-50*c-d)/den*np.exp(-2j*w*left_ps*1e-12)
    s[:,1,1]=(-a+b/50-50*c+d)/den*np.exp(-2j*w*right_ps*1e-12)
    s[:,1,0]=s[:,0,1]=2/den*np.exp(-1j*w*(left_ps+right_ps)*1e-12)
    return s


def fit(source, folder, inductance, calibration):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    net=read_snp(source);f=net.f
    indices=np.unique(np.linspace(0,len(f)-1,min(1501,len(f))).astype(int))
    scale=np.array([100,459,100,56.7,8.2865,20,20])
    low=np.array([0,.1,0,.1,.1,-100,-100])/scale
    high=np.array([1500,3000,1000,3000,100,100,100])/scale
    ls=[inductance,LT1,LR,LR]
    def network(x,ix):
        v=x*scale
        return fast_network(f[ix],ls,[v[0],v[2],REGULAR_C,REGULAR_C],
                       [v[1],v[3],v[4],REGULAR_G],left_ps=v[5],right_ps=v[6])
    def residual(x,ix=indices):
        d=network(x,ix)-net.s[ix]
        # S12 duplicates reciprocal S21 and is not counted twice.
        return np.concatenate([q for i,j in ((0,0),(1,0),(1,1))
                               for q in (d[:,i,j].real,d[:,i,j].imag)])
    seeds=[[100,459,100,56.7,8.2865,6,.7],
           [198,456,139,18.6,8.22,6,.7],
           [100,459,130,70,8.22,6,.7]]
    previous=folder.parent/'heuristic_assumed_L074/fit.json'
    if previous.exists():
        r=json.loads(previous.read_text())
        seeds.append([r['C_fF'][0],r['G_fF'][0],r['C_fF'][1],
                      r['G_fF'][1],r['G_fF'][2],r['delay_ps']-.7,.7])
    fits=[least_squares(residual,np.clip(np.array(s)/scale,low+1e-9,high-1e-9),
                       bounds=(low,high),max_nfev=250,ftol=1e-10,xtol=1e-10,gtol=1e-10)
          for s in seeds]
    best=min(fits,key=lambda q: q.fun@q.fun)
    full=np.arange(len(f));best=least_squares(lambda x:residual(x,full),best.x,
        bounds=(low,high),max_nfev=250,ftol=1e-10,xtol=1e-10,gtol=1e-10)
    v=best.x*scale;m=network(best.x,full)
    result=dict(C0_T2_fF=float(v[0]),Cg_T2_fF=float(v[1]),C0_T1_fF=float(v[2]),
                Cg_T1_fF=float(v[3]),Cg_unit_fF=float(v[4]),
                left_delay_ps=float(v[5]),right_delay_ps=float(v[6]),
                fixed_L_nH=ls,calibration=str(calibration),
                source=str(source),source_sha256=digest(source),
                S21_complex_RMS=float(np.sqrt(np.mean(abs(m[:,1,0]-net.s[:,1,0])**2))),
                active_bounds=best.active_mask.tolist(),optimizer_success=bool(best.success),
                scope='Heuristic fit of complex S11/S21/S22, with all five taper capacitors free; L and regular values fixed. Not target acceptance.')
    (folder/'heuristic.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez(folder/'response.npz',f=f,actual=net.s,model=m)
    print('HEURISTIC',json.dumps(result),flush=True)
    return result
