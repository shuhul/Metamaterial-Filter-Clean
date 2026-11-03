"""Offline diagnostic: does the saved T1 two-port transfer under changed loading?

Reference EM versus its ideal LC fit is compared without refitting T1. A
hybrid cascade is diagnostic only: its input fixture capacitance is not
deembedded, so its fitted T2 values must not be treated as final C0/Cg.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json
import numpy as np
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from circuit_common import read_snp
from geometry.taper.locked_taper2 import ROOT,VERIFIED,FIXED
from geometry.taper.fitting import chain


def s_to_abcd(s,z=50):
    a,b,c,d=s[:,0,0],s[:,0,1],s[:,1,0],s[:,1,1]
    m=np.empty_like(s)
    m[:,0,0]=((1+a)*(1-d)+b*c)/(2*c)
    m[:,0,1]=z*((1+a)*(1+d)-b*c)/(2*c)
    m[:,1,0]=((1-a)*(1-d)-b*c)/(2*z*c)
    m[:,1,1]=((1-a)*(1+d)+b*c)/(2*c)
    return m


def abcd_to_s(m,z=50):
    a,b,c,d=m[:,0,0],m[:,0,1],m[:,1,0],m[:,1,1]
    den=a+b/z+c*z+d;s=np.empty_like(m)
    s[:,1,0]=2/den;s[:,0,1]=2*(a*d-b*c)/den
    s[:,0,0]=(a+b/z-c*z-d)/den;s[:,1,1]=(-a+b/z-c*z+d)/den
    return s


def strip_delay(f,s):
    s=s.copy();w=2*np.pi*f
    left=FIXED['fixture_delay_left_ps']*1e-12;right=FIXED['fixture_delay_right_ps']*1e-12
    s[:,0,0]*=np.exp(2j*w*left);s[:,1,1]*=np.exp(2j*w*right)
    s[:,1,0]*=np.exp(1j*w*(left+right));s[:,0,1]*=np.exp(1j*w*(left+right))
    return s


def check(name):
    out=ROOT/'output_data/taper_exploration'/name
    trial=json.loads((out/'trial_locked.json').read_text());net=read_snp(trial['snp'])
    ref=read_snp(VERIFIED['artifacts']['touchstone'])
    bare=strip_delay(ref.f,ref.s);matrix=s_to_abcd(bare)
    if not np.allclose(abcd_to_s(matrix),bare,rtol=1e-7,atol=1e-9):raise ValueError('Two-port conversion failed')
    fixed=chain(2*np.pi*ref.f,[FIXED['L_taper_nH'],FIXED['L_regular_nH']],
                [FIXED['C_taper_fF'],FIXED['C_regular_fF']],
                [FIXED['Cg_input_fF'],FIXED['Cg_junction_fF'],FIXED['Cg_output_fF']])
    mask=(net.f>=max(ref.f[0],4.55e9))&(net.f<=min(ref.f[-1],5.5e9))
    f=net.f[mask];observed=net.s[mask,1,0]
    downstream=np.empty((len(f),2,2),complex)
    for i in range(2):
        for j in range(2):downstream[:,i,j]=np.interp(f,ref.f,matrix[:,i,j].real)+1j*np.interp(f,ref.f,matrix[:,i,j].imag)
    def response(x):
        w=2*np.pi*f;m=np.tile(np.eye(2,dtype=complex),(len(f),1,1));m[:,0,1]=1/(1j*w*x[1]*1e-15)
        shunt=np.tile(np.eye(2,dtype=complex),(len(f),1,1));shunt[:,1,0]=1j*w*x[0]*1e-15+1/(1j*w*trial['fit']['L_t2_nH']*1e-9)
        s=abcd_to_s(m@shunt@downstream)
        return s[:,1,0]*np.exp(-1j*w*FIXED['effective_fixture_delay_ps']*1e-12)
    def residual(x):
        d=response(x)-observed;return np.r_[d.real,d.imag]
    starts=[[c,g] for c in np.linspace(0,400,11) for g in np.linspace(100,1000,11)]
    starts=sorted(starts,key=lambda x:np.sum(residual(x)**2))[:8]+[[100,459]]
    fits=[least_squares(residual,x,bounds=([0,100],[500,1000]),x_scale='jac',max_nfev=250) for x in starts]
    best=min(fits,key=lambda r:r.fun@r.fun);pred=response(best.x)
    result={'purpose':'Diagnostic hybrid two-port cascade; input fixture not deembedded; NOT a final capacitance extraction',
            'T1_parameters_refitted':False,'hybrid_x':best.x.tolist(),'active_bounds':best.active_mask.tolist(),
            'hybrid_s21_complex_rmse':float(np.sqrt(np.mean(abs(pred-observed)**2))),
            'hybrid_s21_max_complex_error':float(np.max(abs(pred-observed))),
            'reference_ideal_errors':{label:float(np.sqrt(np.mean(abs(bare[:,i,j]-fixed[:,i,j])**2)))
                                      for label,i,j in [('S21',1,0),('S11',0,0),('S22',1,1)]}}
    (out/'reference_transfer_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')
    fig,axes=plt.subplots(2,1,figsize=(10,8))
    axes[0].plot(f/1e9,abs(observed),'r',label='Full EM');axes[0].plot(f/1e9,abs(pred),'k--',label='Hybrid diagnostic');axes[0].legend();axes[0].set_ylabel('|S21|')
    axes[1].plot(ref.f/1e9,abs(bare[:,0,0]-fixed[:,0,0]),label='Reference S11 mismatch')
    axes[1].plot(ref.f/1e9,abs(bare[:,1,0]-fixed[:,1,0]),label='Reference S21 mismatch')
    axes[1].legend();axes[1].set(xlabel='Frequency (GHz)',ylabel='Complex error')
    for ax in axes:ax.grid(alpha=.3)
    fig.suptitle('Reference network transfer diagnostic — not an accepted fit');fig.tight_layout();fig.savefig(out/'reference_transfer_diagnostic.png',dpi=160);plt.close(fig)
    print(json.dumps(result,indent=2));return result


if __name__=='__main__':check('straight4_idc30_span180_ground160_rise300_right584')
