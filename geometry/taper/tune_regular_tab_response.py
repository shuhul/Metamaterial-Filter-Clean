"""Controlled regular-cap correction plus local T2 response refinement."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import numpy as np
from scipy.optimize import differential_evolution,least_squares
from geometry.taper.response_cases import ROOT,OUT,A,W,theory
from geometry.taper.refine_physical_response import features
from geometry.taper.tune_wide_response import impedance
from geometry.taper.optimize_response import THEORY_VERSION

BASE='response_w400_xm160_span206_rise0_right578'

def main():
    m=np.load(OUT/'physical_secant_model.npz');f=m['f'];z=m['z'];dz=list(m['derivatives'])
    s=W['read_network'](ROOT/'output_data/taper_exploration/response_best_regular_tab100/series_gap_probes/series_probe.snp',6)
    pol=np.diag([1,1,-1,-1,-1,-1]);tab_z=A['s_to_z'](pol@s.s@pol)[0]
    difference=tab_z-impedance(BASE,np.array([5e9]))[0]
    dz.append(difference[None,:,:]*(5e9/f)[:,None,None]/100)
    dz=np.array(dz);target=theory(f)[:,1,0];wanted=features(f,target)
    def predict(x):return A['short_internal_z'](z+np.einsum('p,pfij->fij',x,dz))[:,1,0]
    def score(x):
        feat=features(f,predict(x))
        if feat is None:return 1e4
        v=(feat-wanted)/[.015,.015,.10,.10,7,7]
        return float(v@v)
    bounds=[(0,160),(-50,20),(-40,40),(0,80),(0,280)]
    q=differential_evolution(score,bounds,seed=102,popsize=10,maxiter=150,tol=1e-5,polish=True,workers=1)
    def residual(x):
        d=predict(x[:5])-target*np.exp(-2j*np.pi*f*x[5]*1e-12)
        return np.r_[d.real,d.imag]
    fit=least_squares(residual,np.r_[q.x,0],bounds=(np.r_[[b[0] for b in bounds],-30],np.r_[[b[1] for b in bounds],30]),
        x_scale=[80,24,16,40,100,10],max_nfev=400)
    result={'theory_version':THEORY_VERSION,'base_case':BASE,
        'controls':['T2_ground_rise','T2_right_length','T2_IDC_span','T2_wire_offset','regular_cap_tab_width'],
        'feature_best_delta':q.x.tolist(),'feature_score':float(q.fun),
        'feature_prediction':features(f,predict(q.x)).tolist(),
        'complex_best_delta_and_delay':fit.x.tolist(),
        'complex_RMSE':float(np.sqrt(np.mean(fit.fun**2)*2)),
        'complex_prediction_features':features(f,predict(fit.x[:5])).tolist(),
        'target_features':wanted.tolist(),'physical_acceptance':False}
    np.savez(OUT/'tab_physical_model.npz',f=f,z=z,derivatives=dz)
    (OUT/'tab_physical_direction.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
