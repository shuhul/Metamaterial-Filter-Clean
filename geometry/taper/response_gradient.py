"""Cheap physical sensitivities: one 5 GHz solve per geometry perturbation.

This is search guidance only. A proposed geometry must be solved across the
band and then checked with an ordinary unprobed Sonnet sweep.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from dataclasses import replace,asdict
import json
import numpy as np
from scipy.optimize import least_squares,differential_evolution
from scipy.signal import find_peaks
from geometry.taper.response_cases import case,W,A,ROOT,OUT,theory,from_dict

BASE='response_uniform30_gap2_rise358_right536'
CONTROLS=[('idc_span',16),('l_ground_rise',24),('w_ground_riser',-24),('l_right',24),('d_t2_t1',2),('x_meander',-24)]

def perturb(p,key,step):
    return replace(p,taper=replace(p.taper,l_right=p.taper.l_right+step)) if key=='l_right' else replace(p,**{key:getattr(p,key)+step})

def matrices():
    p=from_dict(json.loads((ROOT/'geometry/taper/output'/f'{BASE}.json').read_text())['second'])
    data=np.load(ROOT/'output_data/taper_exploration'/BASE/'response_screen.npz')
    powers=data['powers'];coeff=data['coefficients'];f=np.linspace(4.5e9,5.5e9,2001)
    z=(((f/5e9)[:,None]**powers)@coeff).reshape(len(f),6,6)
    z5=((1.**powers)@coeff).reshape(6,6)
    polarity=np.diag([1,1,-1,-1,-1,-1]);derivatives=[]
    for key,step in CONTROLS:
        name=f'{BASE}_grad_{key}'
        print('Sensitivity:',name,flush=True)
        candidate=replace(p,w_ground_riser=56,x_meander=p.x_meander+step) if key=='x_meander' else perturb(p,key,step)
        reference=z5
        if key=='x_meander':
            # Use the already measured narrow pedestal so an offset wire fits.
            candidate=replace(candidate,w_ground_riser=56)
            prev=ROOT/'output_data/taper_exploration'/f'{BASE}_grad_w_ground_riser/series_gap_probes/series_probe.snp'
            width_net=W['read_network'](prev,6)
            reference=A['s_to_z'](polarity@width_net.s@polarity)[0]
        net=case(name,candidate,frequencies=(5.,),analyze_screen=False)
        if len(net.f)!=1:raise ValueError('Expected one frequency for derivative')
        delta=A['s_to_z'](polarity@net.s@polarity)[0]-reference
        # Capacitance-dominated changes scale as 1/f, length/L changes as f.
        # This approximation is always tested in a later full-band solve.
        scale=(f/5e9) if key in ('l_ground_rise','x_meander') else (5e9/f)
        derivatives.append(delta[None,:,:]*scale[:,None,None]/step)
    np.savez(OUT/'physical_gradient.npz',f=f,z=z,derivatives=np.array(derivatives),
        keys=np.array([k for k,_ in CONTROLS]),steps=np.array([s for _,s in CONTROLS]))
    return p,f,z,np.array(derivatives)

def optimize():
    p,f,z,dz=matrices();target=theory(f)[:,1,0]
    def prediction(delta):return A['short_internal_z'](z+np.einsum('p,pfij->fij',delta,dz))[:,1,0]
    target_peaks,_=find_peaks(abs(target),prominence=.02)
    def peak_features(s):
        ix,_=find_peaks(abs(s),prominence=.01)
        if len(ix)<2:return None
        ix=np.sort(ix[np.argsort(abs(s[ix]))[-2:]])
        positions=[];heights=[]
        for i in ix:
            # Parabolic sub-grid estimate for a smooth search objective.
            a,b,c=abs(s[i-1:i+2])**2
            shift=np.clip(.5*(a-c)/(a-2*b+c),-.5,.5)
            positions.append((f[i]+shift*(f[1]-f[0]))/1e9)
            heights.append(np.sqrt(max(0,b-.25*(a-c)*shift)))
        return np.r_[positions,heights]
    wanted=peak_features(target)
    def objective(delta):
        half_width=(80+delta[2])/2
        if delta[5]-half_width<=-54 or delta[5]+half_width>=58:return 1e4
        features=peak_features(prediction(delta))
        if features is None:return 1e4
        r=(features-wanted)/[.02,.02,.2,.2]
        return float(r@r)
    bounds=[(-64,80),(-150,150),(-56,24),(-200,60),(0,8),(-40,12)]
    q=differential_evolution(objective,bounds,seed=37,popsize=10,maxiter=100,tol=.0001,polish=True,workers=1)
    # Refine the full complex response, allowing only a reference-plane delay.
    def residual(x):
        pred=prediction(x[:6]);expected=target*np.exp(-2j*np.pi*f*x[6]*1e-12)
        r=pred-expected
        return np.r_[r.real,r.imag]
    fit=least_squares(residual,np.r_[q.x,0],bounds=(np.r_[[b[0] for b in bounds],-30],np.r_[[b[1] for b in bounds],30]),
        x_scale=[16,24,24,24,2,24,10],max_nfev=350)
    report={'base_case':BASE,'controls':CONTROLS,'feature_best_delta':q.x.tolist(),
            'feature_objective':float(q.fun),'feature_prediction':peak_features(prediction(q.x)).tolist(),
            'complex_refined_delta_and_delay':fit.x.tolist(),
            'complex_refined_RMSE':float(np.sqrt(np.mean(fit.fun**2)*2)),
            'target_features':wanted.tolist(),'physical_acceptance':False}
    (OUT/'physical_gradient_direction.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':optimize()
