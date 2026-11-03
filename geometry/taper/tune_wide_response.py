"""Local physical response search after separating the T2 and T1 inductors."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from dataclasses import replace
import json
import numpy as np
from scipy.optimize import differential_evolution,least_squares
from scipy.signal import find_peaks,peak_widths
from geometry.taper.response_cases import ROOT,OUT,W,A,case,from_dict,theory
from geometry.taper.optimize_response import THEORY_VERSION

BASE='response_w400_wireleft120_uniform30_gap6_pad24_rise250'
OTHER='response_w400_wireleft120_uniform30_gap6_pad24'

def impedance(name,f):
    m=np.load(ROOT/'output_data/taper_exploration'/name/'response_screen.npz')
    return (((f/5e9)[:,None]**m['powers'])@m['coefficients']).reshape(len(f),6,6)

def main():
    p=from_dict(json.loads((ROOT/'geometry/taper/output'/f'{BASE}.json').read_text())['second'])
    f=np.linspace(4.5e9,5.5e9,2001);z=impedance(BASE,f);z5=impedance(BASE,np.array([5e9]))[0]
    dz=[(impedance(OTHER,f)-z)/(358-250)]
    pol=np.diag([1,1,-1,-1,-1,-1])
    for key,step in [('l_right',24),('idc_span',16),('x_meander',40)]:
        pp=replace(p,taper=replace(p.taper,l_right=p.taper.l_right+step)) if key=='l_right' else replace(p,**{key:getattr(p,key)+step})
        name='response_wide250_grad_'+key
        print('Wide sensitivity',name,flush=True)
        net=case(name,pp,frequencies=(5.,),analyze_screen=False)
        delta=A['s_to_z'](pol@net.s@pol)[0]-z5
        scale=f/5e9 if key=='x_meander' else 5e9/f
        dz.append(delta[None,:,:]*scale[:,None,None]/step)
    dz=np.array(dz)
    target=theory(f)[:,1,0]
    def features(s):
        ix,_=find_peaks(abs(s),prominence=.015)
        if len(ix)!=2:return None
        widths=peak_widths(abs(s)**2,ix,rel_height=.5)[0]*(f[1]-f[0])/1e6
        pos=[];heights=[]
        for i in ix:
            aa,bb,cc=abs(s[i-1:i+2])**2
            shift=np.clip(.5*(aa-cc)/(aa-2*bb+cc),-.5,.5)
            pos.append((f[i]+shift*(f[1]-f[0]))/1e9)
            heights.append(np.sqrt(max(0,bb-.25*(aa-cc)*shift)))
        return np.r_[pos,heights,widths]
    wanted=features(target)
    def predict(x):return A['short_internal_z'](z+np.einsum('p,pfij->fij',x,dz))[:,1,0]
    def score(x):
        v=features(predict(x))
        if v is None:return 1e4
        r=(v-wanted)/[.025,.025,.18,.18,10,10]
        return float(r@r)
    bounds=[(-250,150),(-136,64),(-48,80),(-40,80)]
    q=differential_evolution(score,bounds,seed=54,popsize=10,maxiter=130,tol=1e-5,polish=True,workers=1)
    def residual(x):
        s=predict(x[:4]);t=target*np.exp(-2j*np.pi*f*x[4]*1e-12);d=s-t
        return np.r_[d.real,d.imag]
    fit=least_squares(residual,np.r_[q.x,0],bounds=(np.r_[[b[0] for b in bounds],-30],np.r_[[b[1] for b in bounds],30]),
        x_scale=[80,24,16,40,10],max_nfev=400)
    report={'base_case':BASE,'theory_version':THEORY_VERSION,'controls':['ground_rise','right_length','IDC_span','wire_offset'],
            'feature_best_delta':q.x.tolist(),'feature_score':float(q.fun),
            'predicted_features':features(predict(q.x)).tolist(),'target_features':wanted.tolist(),
            'complex_best_delta_and_delay':fit.x.tolist(),'complex_RMSE':float(np.sqrt(np.mean(fit.fun**2)*2)),
            'physical_acceptance':False}
    np.savez(OUT/'wide_physical_gradient.npz',f=f,z=z,derivatives=dz)
    (OUT/'wide_physical_direction.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)

if __name__=='__main__':main()
