"""Use driven-node ABCD invariants to separate downstream Cg from input LC.

After removing the feed and the known downstream chain, an ideal directly
driven shunt has matrix [[1,0],[Y,1]]. The three invariant entries determine
one downstream coupling without assuming the directly driven node's C0.
The remaining entry is checked independently against known L and one C0.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from geometry.taper.staged_fitting import circuit
from geometry.taper.paired_coupling import abcd


def normalized(m):
    m=m.copy();m[:,0,1]/=50;m[:,1,0]*=50
    return m


def feed(f,ps,z):
    angle=2*np.pi*f*ps*1e-12
    m=np.zeros((len(f),2,2),complex)
    m[:,0,0]=m[:,1,1]=np.cos(angle)
    m[:,0,1]=1j*(z/50)*np.sin(angle)
    m[:,1,0]=1j*(50/z)*np.sin(angle)
    return m


def extract(f,s,ls,cs,couplings,known_input_l_nh,target_cg,out,*,left_z=56.8872):
    """couplings[0] is the ONLY electrical fit variable; all others fixed."""
    f=np.asarray(f);measured=normalized(abcd(s))
    # Keep all samples; a narrow resonance may be the most informative data.
    def isolate(x):
        downstream=normalized(abcd(circuit(f,ls,cs,[x[0]]+list(couplings[1:]),
            right_ps=x[2],line_z=56.8872)))
        inverse=downstream[:,::-1,::-1].copy()
        inverse[:,0,1]=-downstream[:,0,1];inverse[:,1,0]=-downstream[:,1,0]
        q=feed(f,-x[1],left_z)@measured@inverse
        return q,downstream
    def residual(x):
        q,n=isolate(x)
        scale=np.maximum(1,np.sqrt(np.sum(abs(n)**2,axis=(1,2))))
        delta=np.c_[q[:,0,0]-1,q[:,0,1],q[:,1,1]-1]/scale[:,None]
        return np.r_[delta.real.ravel(),delta.imag.ravel()]
    lower=[target_cg*.5,0,0];upper=[target_cg*1.5,12,8]
    fits=[least_squares(residual,[g,3.5,1.0],bounds=(lower,upper),
          x_scale=[target_cg,4,1],max_nfev=180,ftol=1e-10,xtol=1e-10,gtol=1e-10)
          for g in (target_cg*.8,target_cg,target_cg*1.2)]
    best=min(fits,key=lambda q:q.fun@q.fun)
    q,_=isolate(best.x);y=q[:,1,0]/50;w=2*np.pi*f
    adjusted=y.imag+1/(w*known_input_l_nh*1e-9)
    c=float(np.sum(w*adjusted)/np.sum(w*w)*1e15)
    ccurve=adjusted/w*1e15
    wanted=1j*w*c*1e-15+1/(1j*w*known_input_l_nh*1e-9)
    singular=np.linalg.svd(best.jac*np.array([target_cg,4,1])[None,:],compute_uv=False)
    report=dict(free_electrical_parameters=['downstream_Cg_fF'],Cg_fF=float(best.x[0]),
        Cg_error_percent=float(100*(best.x[0]/target_cg-1)),
        feed_nuisance_parameters=dict(left_delay_ps=float(best.x[1]),right_delay_ps=float(best.x[2])),
        left_feed_Z_fixed_ohm=left_z,right_feed_Z_fixed_ohm=56.8872,
        weighted_invariant_RMS=float(np.sqrt(np.mean(best.fun**2))),
        independently_fitted_C0_fF=c,C0_single_parameter_given_L_nH=known_input_l_nh,
        C0_frequency_range_10_90_percentile_fF=np.quantile(ccurve,[.1,.9]).tolist(),
        node_Y_complex_RMS_S=float(np.sqrt(np.mean(abs(y-wanted)**2))),
        node_Y_real_RMS_S=float(np.sqrt(np.mean(y.real**2))),
        active_bounds=best.active_mask.tolist(),jacobian_condition=float(singular[0]/singular[-1]),
        samples=len(f),acceptance_pending=True,
        limitation='Requires a separable feed, local shunt and downstream block; residuals and feed-model sensitivity must pass before use')
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    (out/'driven_node_extraction.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez(out/'driven_node_data.npz',f=f,q=q,node_y=y,C0_curve_fF=ccurve)
    print(json.dumps(report,indent=2),flush=True)
    return report
