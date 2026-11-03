"""Direct response matching with measured local T2 and unit-cap sensitivities."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json,itertools
import numpy as np
from scipy.optimize import differential_evolution,least_squares,minimize_scalar
from geometry.taper.response_cases import ROOT,OUT,A,W,theory
from geometry.taper.refine_physical_response import features
from geometry.taper.tune_wide_response import impedance
from geometry.taper.optimize_response import THEORY_VERSION

BASE='response_w400_xm160_span206_rise0_right578'

def main():
    old=np.load(OUT/'physical_secant_model.npz');f=old['f'];z=old['z'];z5=impedance(BASE,np.array([5e9]))[0]
    pol=np.diag([1,1,-1,-1,-1,-1])
    def derivative(name,step,inductive=False):
        folder=ROOT/'output_data/taper_exploration'/name/'series_gap_probes'
        saved=json.loads((folder/'screen_completed.json').read_text())
        path=folder/'series_probe.snp'
        if W['digest'](path)!=saved['snp_sha256']:raise ValueError('Sensitivity data changed')
        net=W['read_network'](path,6)
        if len(net.f)!=1 or abs(net.f[0]-5e9)>1:raise ValueError('Expected 5 GHz sensitivity')
        delta=A['s_to_z'](pol@net.s@pol)[0]-z5
        scale=f/5e9 if inductive else 5e9/f
        return delta[None,:,:]*scale[:,None,None]/step
    dz=np.array([derivative('response_best_rise40',40,True),old['derivatives'][1],old['derivatives'][2],
                 derivative('response_best_idc34',4),derivative('response_best_regular_tab100_t4',100)])
    target=theory(f);wanted=features(f,target[:,1,0])
    def predict(x):return A['short_internal_z'](z+np.einsum('p,pfij->fij',x,dz))
    def feature_score(x):
        feat=features(f,predict(x)[:,1,0])
        if feat is None:return 1e4
        r=(feat-wanted)/[.01,.01,.08,.08,5,5]
        return float(r@r)
    bounds=[(0,80),(-32,16),(-32,56),(0,4),(0,240)]
    q=differential_evolution(feature_score,bounds,seed=126,popsize=10,maxiter=130,tol=1e-5,
        x0=[0,-20,40,0,100],workers=1,polish=True)
    def residual(x):
        d=predict(x[:5])[:,1,0]-target[:,1,0]*np.exp(-2j*np.pi*f*x[5]*1e-12)
        return np.r_[d.real,d.imag]
    fit=least_squares(residual,np.r_[q.x,0],bounds=(np.r_[[b[0] for b in bounds],-30],np.r_[[b[1] for b in bounds],30]),
        x_scale=[40,24,16,4,100,10],max_nfev=400)
    candidates={}
    grid=np.array([2,2,2,1,4])
    for center in (q.x,fit.x[:5]):
        opts=[np.unique([step*np.floor(v/step),step*np.ceil(v/step)]) for v,step in zip(center,grid)]
        for row in itertools.product(*opts):
            x=np.array(row)
            if any(not lo<=v<=hi for v,(lo,hi) in zip(x,bounds)):continue
            key=tuple(row)
            if key in candidates:continue
            pred=predict(x);feat=features(f,pred[:,1,0])
            if feat is None:continue
            def phase_cost(tau):return np.mean(abs(pred[:,1,0]-target[:,1,0]*np.exp(-2j*np.pi*f*tau*1e-12))**2)
            tau=minimize_scalar(phase_cost,bounds=(-30,30),method='bounded')
            rel=100*(feat[2:]/wanted[2:]-1)
            candidates[key]={'delta':x.tolist(),'predicted_features':feat.tolist(),
                'frequency_errors_MHz':((feat[:2]-wanted[:2])*1000).tolist(),
                'height_and_width_errors_percent':rel.tolist(),'feature_score':feature_score(x),
                'complex_RMSE':float(np.sqrt(tau.fun)),'delay_ps':float(tau.x)}
    rows=sorted(candidates.values(),key=lambda a:a['feature_score'])
    result={'theory_version':THEORY_VERSION,'base_case':BASE,
        'controls':['ground_rise','right_length','IDC_span','IDC_pairs_added','regular_cap_tab_width_t4'],
        'feature_best':q.x.tolist(),'feature_score':float(q.fun),
        'complex_best':fit.x.tolist(),'complex_RMSE':float(np.sqrt(np.mean(fit.fun**2)*2)),
        'target_features':wanted.tolist(),'grid_candidates':rows[:12],'physical_acceptance':False}
    np.savez(OUT/'local_tab_response_model.npz',f=f,z=z,derivatives=dz)
    (OUT/'local_tab_response_direction.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
