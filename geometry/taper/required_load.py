"""Invert the T2 branch loading from a desired S21; no EM solves."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import numpy as np
from scipy.optimize import least_squares
from geometry.taper.optimize_response import A,OUT,theory

def required_gamma(z,target):
    # Short the T1, regular, and end gaps; leave the T2 gap accessible.
    keep=[0,1,2];other=[3,4,5]
    zz=z[:,keep,:][:,:,keep]-z[:,keep,:][:,:,other]@np.linalg.solve(z[:,other,:][:,:,other],z[:,other,:][:,:,keep])
    s=A['z_to_s'](zz);h=target-s[:,1,0]
    return h/(s[:,1,2]*s[:,2,0]+s[:,2,2]*h)

def optimize():
    m=np.load(OUT/'physical_gradient.npz');f=m['f'];z=m['z'];dz=m['derivatives']
    target=theory(f)[:,1,0];mask=abs(target)>.025
    f=f[mask];z=z[mask];dz=dz[:,mask];target=target[mask]
    n=len(dz);bounds=[(-64,80),(-150,150),(-56,24),(-200,60),(0,8),(-40,12),(-30,30)]
    def residual(x):
        zz=z+np.einsum('p,pfij->fij',x[:n],dz)
        gamma=required_gamma(zz,target*np.exp(-2j*np.pi*f*x[n]*1e-12))
        d=gamma+1
        return np.r_[d.real,d.imag]
    old=json.loads((OUT/'physical_gradient_direction.json').read_text())
    seeds=[np.zeros(n+1),np.r_[old['feature_best_delta'],0],np.array(old['complex_refined_delta_and_delay'])]
    rng=np.random.default_rng(81)
    seeds += [np.array([rng.uniform(lo,hi) for lo,hi in bounds]) for _ in range(5)]
    fits=[least_squares(residual,x,bounds=np.array(bounds).T,x_scale=[16,24,24,24,2,24,10],max_nfev=400) for x in seeds]
    q=min(fits,key=lambda a:a.fun@a.fun)
    # Refine the direct complex objective from the best inverse-load solution.
    def direct(x):
        s=A['short_internal_z'](z+np.einsum('p,pfij->fij',x[:n],dz))[:,1,0]
        d=s-target*np.exp(-2j*np.pi*f*x[n]*1e-12)
        return np.r_[d.real,d.imag]
    refined=least_squares(direct,q.x,bounds=np.array(bounds).T,x_scale=[16,24,24,24,2,24,10],max_nfev=500)
    report={'method':'Required internal reflection equals -1 for the desired response with an unchanged shorted wire',
        'inverse_solution':q.x.tolist(),'inverse_cost':float(q.fun@q.fun),
        'refined_solution':refined.x.tolist(),'peak_region_complex_RMSE':float(np.sqrt(np.mean(refined.fun**2)*2)),
        'all_starts':[{'x':r.x.tolist(),'cost':float(r.fun@r.fun)} for r in fits],
        'physical_acceptance':False}
    (OUT/'inverse_load_direction.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':optimize()
