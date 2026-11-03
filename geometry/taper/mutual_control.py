"""Controlled coupling-first, inductor-second diagnostic of the new fixtures.

Measure Cg with the T2 wire omitted, then fit only mutual L when restoring
the independently measured wire. T1 and unit-cell values are never released.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from circuit_common import read_snp
from geometry.taper.staged_fitting import LT1,CT1,LR,CR,GOUT,save_plot
from geometry.taper.run_bridged_stage import OUT
from geometry.taper.bridged_fixtures import ROOT

GJ=8.40642914296533


def coupled(f,ls,cs,gs,mutual_nh=0,input_cg=None,left_ps=0,right_ps=0):
    """Nodal solution of parallel LCs, adjacent Cg and mutual L of first two."""
    f=np.asarray(f);w=2*np.pi*f;n=len(ls)
    if len(cs)!=n or len(gs)!=n-1:raise ValueError('Invalid ladder sizes')
    inverse=np.diag([0 if np.isinf(l) else 1/(l*1e-9) for l in ls])
    if mutual_nh:
        pair=np.array([[ls[0],mutual_nh],[mutual_nh,ls[1]]])*1e-9
        if np.linalg.eigvalsh(pair).min()<=0:raise ValueError('Inductance matrix must be positive definite')
        inverse[:2,:2]=np.linalg.inv(pair)
    offset=1 if input_cg is not None else 0;size=n+offset
    c=np.diag(([0] if offset else [])+list(cs))*1e-15
    edges=([input_cg] if offset else [])+list(gs)
    for i,g in enumerate(edges):
        v=g*1e-15;c[i,i]+=v;c[i+1,i+1]+=v;c[i,i+1]-=v;c[i+1,i]-=v
    inductive=np.zeros((size,size));inductive[offset:,offset:]=inverse
    y=1j*w[:,None,None]*c-1j*inductive[None,:,:]/w[:,None,None]
    port=[0,size-1];internal=list(range(1,size-1))
    reduced=y[:,port,:][:,:,port]
    if internal:
        reduced-=y[:,port,:][:,:,internal]@np.linalg.solve(y[:,internal,:][:,:,internal],y[:,internal,:][:,:,port])
    eye=np.eye(2)[None,:,:]
    s=(eye-50*reduced)@np.linalg.inv(eye+50*reduced)
    s[:,0,0]*=np.exp(-4j*np.pi*f*left_ps*1e-12)
    s[:,1,1]*=np.exp(-4j*np.pi*f*right_ps*1e-12)
    s[:,0,1]*=np.exp(-2j*np.pi*f*(left_ps+right_ps)*1e-12)
    s[:,1,0]*=np.exp(-2j*np.pi*f*(left_ps+right_ps)*1e-12)
    return s


def fit_one(snp,out,*,kind='coupling',fixed_g=56.7,l2=.7402929008862823,c2=100,signed_delay=False):
    net=read_snp(snp);f=net.f;observed=net.s
    coarse=np.unique(np.linspace(0,len(f)-1,min(1901,len(f))).astype(int))
    def predict(x,indices):
        g=x[0] if kind=='coupling' else fixed_g
        m=0 if kind=='coupling' else x[0]
        return coupled(f[indices],[np.inf if kind=='coupling' else l2,LT1,LR,LR],
            [c2,CT1,CR,CR],[g,GJ,GOUT],mutual_nh=m,left_ps=x[1]*.8,right_ps=x[1]*.2)
    def residual(x,indices=coarse):
        d=predict(x,indices)[:,1,0]-observed[indices,1,0]
        return np.r_[d.real,d.imag]
    bounds=([15,0],[110,25]) if kind=='coupling' else ([-.65,-30 if signed_delay else 0],[.65,30 if signed_delay else 25])
    grid=np.linspace(bounds[0][0],bounds[1][0],90)
    def amplitude_score(v):
        return np.mean((abs(predict([v,0],coarse)[:,1,0])-abs(observed[coarse,1,0]))**2)
    seeds=list(sorted(grid,key=lambda v:np.mean(residual([v,5])**2))[:3])+list(sorted(grid,key=amplitude_score)[:3])
    starts=[(v,tau) for v in seeds for tau in ((-10,5,15) if signed_delay else (5,))]
    fits=[least_squares(residual,[v,tau],bounds=bounds,x_scale=[56 if kind=='coupling' else .2,5],
        max_nfev=250,ftol=1e-10,xtol=1e-10,gtol=1e-10) for v,tau in starts]
    best=min(fits,key=lambda q:q.fun@q.fun)
    full=np.arange(len(f))
    best=least_squares(lambda x:residual(x,full),best.x,bounds=bounds,
        x_scale=[56 if kind=='coupling' else .2,5],max_nfev=150,ftol=1e-11,xtol=1e-11,gtol=1e-11)
    model=predict(best.x,full)
    key='Cg_T2_T1_fF' if kind=='coupling' else 'M_T2_T1_nH'
    report=dict(free_electrical_parameters=[key],value=float(best.x[0]),
        target=56.7 if kind=='coupling' else 0.,
        target_error_percent=float(100*(best.x[0]/56.7-1)) if kind=='coupling' else None,
        fixture_delay_ps=float(best.x[1]),active_bounds=best.active_mask.tolist(),
        signed_residual_delay_allowed=signed_delay,
        fixed=dict(L_T2_nH=None if kind=='coupling' else l2,C0_T2_provisional_fF=c2,
                   L_T1_nH=LT1,C0_T1_fF=CT1,L_regular_nH=LR,C0_regular_fF=CR,
                   Cg_T1_regular_fF=GJ,Cg_output_fF=GOUT,Cg_T2_T1_fF=None if kind=='coupling' else fixed_g),
        S21_complex_RMSE=float(np.sqrt(np.mean(abs(model[:,1,0]-observed[:,1,0])**2))),
        S21_magnitude_RMSE=float(np.sqrt(np.mean((abs(model[:,1,0])-abs(observed[:,1,0]))**2))),
        S11_complex_RMSE=float(np.sqrt(np.mean(abs(model[:,0,0]-observed[:,0,0])**2))),
        S22_complex_RMSE=float(np.sqrt(np.mean(abs(model[:,1,1]-observed[:,1,1])**2))),
        acceptance_pending=True,source=str(Path(snp).resolve()),
        note='Controlled diagnostic. Provisional T2 C0 and fixture-model sensitivity remain to be checked.')
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    (out/'fit.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez(out/'response.npz',f=f,actual=observed,model=model)
    save_plot(f,observed,model,report,out,'T2 Coupling Control' if kind=='coupling' else 'T2 Inductor Interaction Diagnostic')
    print(json.dumps(report,indent=2),flush=True)
    return report
