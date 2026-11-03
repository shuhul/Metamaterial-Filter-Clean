"""Update the measured geometry Jacobian with a full-band EM secant."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from dataclasses import replace
import json,itertools
import numpy as np
from scipy.optimize import differential_evolution,least_squares,minimize_scalar
from scipy.signal import find_peaks,peak_widths
from geometry.taper.response_cases import ROOT,OUT,A,from_dict,theory
from geometry.taper.optimize_response import THEORY_VERSION
from geometry.taper.tune_wide_response import BASE,impedance

NEW='response_w400_xm160_span206_rise0_right578'

def features(f,s):
    ix,_=find_peaks(abs(s),prominence=.015)
    if len(ix)!=2:return None
    widths=peak_widths(abs(s)**2,ix,rel_height=.5)[0]*(f[1]-f[0])/1e6
    pos=[];heights=[]
    for i in ix:
        aa,bb,cc=abs(s[i-1:i+2])**2;shift=np.clip(.5*(aa-cc)/(aa-2*bb+cc),-.5,.5)
        pos.append((f[i]+shift*(f[1]-f[0]))/1e9)
        heights.append(np.sqrt(max(0,bb-.25*(aa-cc)*shift)))
    return np.r_[pos,heights,widths]

def main():
    p=from_dict(json.loads((ROOT/'geometry/taper/output'/f'{NEW}.json').read_text())['second'])
    p0=from_dict(json.loads((ROOT/'geometry/taper/output'/f'{BASE}.json').read_text())['second'])
    m=np.load(OUT/'wide_physical_gradient.npz');f=m['f'];z=impedance(NEW,f);z0=m['z'];dz=m['derivatives'].copy()
    delta=np.array([p.l_ground_rise-p0.l_ground_rise,p.taper.l_right-p0.taper.l_right,
                    p.idc_span-p0.idc_span,p.x_meander-p0.x_meander])
    scales=np.array([80,24,16,40]);u=delta/scales
    error=z-z0-np.einsum('p,pfij->fij',delta,dz)
    dz+=((u/(u@u)/scales)[:,None,None,None])*error[None,:,:,:]
    target=theory(f)[:,1,0];wanted=features(f,target)
    def prediction(x):return A['short_internal_z'](z+np.einsum('p,pfij->fij',x,dz))[:,1,0]
    def score(x):
        feat=features(f,prediction(x))
        if feat is None:return 1e4
        d=(feat-wanted)/[.015,.015,.1,.1,7,7]
        return float(d@d)
    bounds=[(0,160),(-50,20),(-40,40),(0,80)]
    q=differential_evolution(score,bounds,seed=93,popsize=10,maxiter=120,tol=1e-5,polish=True,workers=1)
    def residual(x):
        d=prediction(x[:4])-target*np.exp(-2j*np.pi*f*x[4]*1e-12)
        return np.r_[d.real,d.imag]
    fit=least_squares(residual,np.r_[q.x,0],bounds=(np.r_[[b[0] for b in bounds],-30],np.r_[[b[1] for b in bounds],30]),
        x_scale=[80,24,16,40,10],max_nfev=400)
    # Evaluate nearby actual 2 um-grid geometries, including phase alignment.
    choices=[]
    for center in (q.x,fit.x[:4]):
        rounded=2*np.floor(center/2+.5)
        for extra in itertools.product((-2,0,2),repeat=4):
            x=rounded+extra
            if any(not lo<=v<=hi for v,(lo,hi) in zip(x,bounds)):continue
            pp=replace(p,l_ground_rise=int(p.l_ground_rise+x[0]),idc_span=int(p.idc_span+x[2]),
                x_meander=int(p.x_meander+x[3]),taper=replace(p.taper,l_right=int(p.taper.l_right+x[1])))
            s=prediction(x)
            predicted_features=features(f,s)
            if predicted_features is None:continue
            def delay_cost(tau):return np.mean(abs(s-target*np.exp(-2j*np.pi*f*tau*1e-12))**2)
            tau=minimize_scalar(delay_cost,bounds=(-30,30),method='bounded')
            choices.append({'delta':x.tolist(),'complex_RMSE':float(np.sqrt(tau.fun)),
                'delay_ps':float(tau.x),'feature_score':score(x),'predicted_features':predicted_features.tolist()})
    choices.sort(key=lambda a:a['complex_RMSE'])
    report={'base':NEW,'theory_version':THEORY_VERSION,'feature_best':q.x.tolist(),'feature_score':float(q.fun),
        'complex_best':fit.x.tolist(),'complex_RMSE':float(np.sqrt(np.mean(fit.fun**2)*2)),
        'target_features':wanted.tolist(),'grid_candidates':choices[:8],'physical_acceptance':False}
    (OUT/'physical_secant_direction.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez(OUT/'physical_secant_model.npz',f=f,z=z,derivatives=dz)
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
