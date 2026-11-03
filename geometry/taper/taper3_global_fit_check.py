"""Broad numerical check of a saved T3 fit. No EM or production-fit changes."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json,time
import numpy as np
from scipy.optimize import differential_evolution,least_squares
from geometry.taper.taper3 import OUT,digest
from geometry.taper.progressive_fitting import transmission


def check(case):
    folder=OUT/case/'response'
    fit=json.loads((folder/'fit/fit.json').read_text())
    data=np.load(folder/'fit/response.npz')
    ix=np.unique(np.linspace(0,len(data['f'])-1,1001).astype(int))
    f=data['f'][ix];obs=data['actual'][ix,1,0]
    bounds=[(0.,1000.),(.1,3000.),(.1,3000.),(-100.,100.)]
    scales=[251.,261.,153.,20.]

    def residual(x,smooth=False,dense=False):
        cs=fit['C_fF'].copy();gs=fit['G_fF'].copy()
        cs[0]=x[0];gs[:2]=x[1:3]
        ff=data['f'] if dense else f
        actual=data['actual'][:,1,0] if dense else obs
        model=transmission(ff,fit['L_nH'],cs,gs,x[-1])
        delta=((1/model-1/actual)*np.maximum(abs(actual),.03)**2
               if smooth else model-actual)
        return np.r_[delta.real,delta.imag]

    start=time.perf_counter();rows=[]
    for smooth,seed in [(False,740),(True,740),(True,741),(True,742)]:
        def objective(x):
            r=residual(x,smooth=smooth)
            return r@r
        de=differential_evolution(objective,bounds,seed=seed,popsize=14,
            maxiter=150,tol=1e-8,polish=False,workers=1)
        local=least_squares(residual,de.x,bounds=np.array(bounds).T,
            x_scale=scales,max_nfev=1000,ftol=1e-11,xtol=1e-11,gtol=1e-11)
        q=least_squares(lambda x:residual(x,dense=True),local.x,
            bounds=np.array(bounds).T,x_scale=scales,max_nfev=600,
            ftol=1e-11,xtol=1e-11,gtol=1e-11)
        rows.append(dict(smooth_search=smooth,seed=seed,values=q.x.tolist(),
            complex_RMS=float(np.sqrt(np.mean(q.fun**2)*2)),
            active_bounds=q.active_mask.tolist(),global_evaluations=de.nfev))
    result=dict(case=case,saved_RMS=fit['complex_RMS'],searches=rows,
        seconds=time.perf_counter()-start,production_fit_unchanged=True,
        source_fit_sha256=digest(folder/'fit/fit.json'),code_sha256=digest(__file__),
        scope='All older values/L fixed. Diagnostic searches, not proof of global optimality.')
    (folder/'global_fit_check.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result
