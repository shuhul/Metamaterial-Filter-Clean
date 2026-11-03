"""Offline mutual-inductance hypothesis, never a final two-variable extraction.

The production fit remains unchanged. M is deliberately released ONLY here to
test whether adding an inductor next to T1 introduced an omitted coupling path.
An M inferred from this same S21 cannot subsequently be called an independently
known constant to disguise a three-parameter fit as a two-parameter fit.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import numpy as np
from scipy.optimize import least_squares
from geometry.taper.locked_taper2 import ROOT,FIXED,response
from circuit_common import read_snp


def coupled_response(f,x,l2,*,cross_capacitance=None):
    """x = C0_t2 [fF], Cg_input [fF], signed M_t2_t1 [nH]."""
    w=2*np.pi*f
    inductance=np.diag([l2,FIXED['L_taper_nH'],FIXED['L_regular_nH']])*1e-9
    inductance[0,1]=inductance[1,0]=x[2]*1e-9
    if np.linalg.eigvalsh(inductance).min()<=0:raise ValueError('Nonphysical inductance matrix')
    # Nodes 0/4 are the ports; 1/2/3 are T2/T1/regular.
    y=np.zeros((len(f),5,5),complex)
    for i,c in enumerate([x[1],FIXED['Cg_input_fF'],FIXED['Cg_junction_fF'],FIXED['Cg_output_fF']]):
        v=1j*w*c*1e-15
        y[:,i,i]+=v;y[:,i+1,i+1]+=v;y[:,i,i+1]-=v;y[:,i+1,i]-=v
    if cross_capacitance is not None:
        i,j,c=cross_capacitance;v=1j*w*c*1e-15
        y[:,i,i]+=v;y[:,j,j]+=v;y[:,i,j]-=v;y[:,j,i]-=v
    y[:,1:4,1:4]+=np.linalg.inv(inductance)[None,:,:]/(1j*w[:,None,None])
    for i,c in enumerate([x[0],FIXED['C_taper_fF'],FIXED['C_regular_fF']],1):y[:,i,i]+=1j*w*c*1e-15
    pp=y[:,[0,4],:][:,:,[0,4]]
    pi=y[:,[0,4],1:4];ip=y[:,1:4,:][:,:,[0,4]]
    effective=pp-pi@np.linalg.solve(y[:,1:4,1:4],ip)
    eye=np.tile(np.eye(2),(len(f),1,1))
    s=(eye-50*effective)@np.linalg.inv(eye+50*effective)
    left=FIXED['fixture_delay_left_ps']*1e-12;right=FIXED['fixture_delay_right_ps']*1e-12
    s[:,0,0]*=np.exp(-2j*w*left);s[:,1,1]*=np.exp(-2j*w*right)
    s[:,1,0]*=np.exp(-1j*w*(left+right));s[:,0,1]*=np.exp(-1j*w*(left+right))
    return s


def diagnose(name):
    out=ROOT/'output_data/taper_exploration'/name
    trial=json.loads((out/'trial_locked.json').read_text());n=read_snp(out/'response.snp')
    f=n.f;l2=trial['fit']['L_t2_nH'];obs=n.s[:,1,0]
    if not np.allclose(coupled_response(f,[100,459,0],l2),response(f,[100,459],l2),atol=1e-8):
        raise ValueError('Zero-M model failed production-chain equivalence')
    def residual(x):
        d=coupled_response(f,x,l2)[:,1,0]-obs
        return np.r_[d.real,d.imag]
    grid=[np.array([c,g,m]) for c in (10,100,250) for g in (100,300,459,800)
          for m in np.linspace(-.2,.2,81)]
    grid.sort(key=lambda x:float(residual(x)@residual(x)))
    starts=grid[:10]+[np.r_[trial['fit']['x'],0],np.array([100,459,0])]
    fits=[least_squares(residual,x,bounds=([0,30,-.5],[400,1500,.5]),
                        x_scale=[100,459,.05],max_nfev=350) for x in starts]
    best=min(fits,key=lambda r:r.fun@r.fun);d=coupled_response(f,best.x,l2)[:,1,0]-obs
    result={'purpose':'Three-parameter mutual-coupling diagnostic, NOT an accepted two-variable fit',
            'case':name,'free_parameters':['C0_t2_fF','Cg_t2_fF','M_t2_t1_nH'],
            'x':best.x.tolist(),'active_bounds':best.active_mask.tolist(),
            'zero_M_identity_checked':True,'complex_rmse':float(np.sqrt(np.mean(abs(d)**2))),
            'max_complex_error':float(max(abs(d)))}
    (out/'mutual_inductance_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result


if __name__=='__main__':
    from geometry.taper.run_locked_taper2 import CASE_NAME
    diagnose(CASE_NAME)
