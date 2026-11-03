"""Full-band side-capacitor control with the T2 wire removed.

Fit only the remaining source island's parasitic C and side coupling Cg.
T1 and all downstream electrical constants remain fixed. Delay is explicit.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from circuit_common import read_snp
from geometry.taper.compact_fit import GJ
from geometry.taper.paired_coupling import abcd
from geometry.taper.staged_fitting import circuit,LT1,CT1,LR,CR,GOUT
from geometry.taper.fit_restored_t2 import features
from geometry.taper.run_bridged_stage import digest
import matplotlib.pyplot as plt


def predictor(f):
    w=2*np.pi*np.asarray(f)
    n=abcd(circuit(f,[LT1,LR,LR],[CT1,CR,CR],[None,GJ,GOUT]))
    def predict(x,full=False):
        cap,g,tau=x;zz=1/(1j*w*g*1e-15);yy=1j*w*cap*1e-15
        a=n[:,0,0]+zz*n[:,1,0];b=n[:,0,1]+zz*n[:,1,1]
        c=n[:,1,0]+yy*a;d=n[:,1,1]+yy*b
        den=a+b/50+c*50+d
        transmission=2/den*np.exp(-1j*w*tau*1e-12)
        if not full:return transmission
        s=np.empty_like(n);s[:,1,0]=s[:,0,1]=transmission
        s[:,0,0]=(a+b/50-c*50-d)/den*np.exp(-2j*w*tau*1e-12)
        s[:,1,1]=(-a+b/50-c*50+d)/den
        return s
    return predict


def fit(source,out):
    source=Path(source);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    n=read_snp(source);ix=np.unique(np.linspace(0,len(n.f)-1,min(1501,len(n.f))).astype(int))
    pred=predictor(n.f[ix]);obs=n.s[ix,1,0]
    def res(x):
        error=pred(x)-obs;return np.r_[error.real,error.imag]
    grid=[]
    for c in (0,30,80,160,300):
        for g in np.arange(35,80,.1):
            x=[c,g,5];r=res(x);grid.append((r@r,x))
    grid.sort(key=lambda q:q[0]);bounds=([0,15,-30],[2000,150,30])
    fits=[least_squares(res,x,bounds=bounds,x_scale=[100,57,5],max_nfev=250) for _,x in grid[:8]]
    best=min(fits,key=lambda q:q.fun@q.fun)
    full=predictor(n.f)
    def residual(x):
        d=full(x)-n.s[:,1,0];return np.r_[d.real,d.imag]
    best=least_squares(residual,best.x,bounds=bounds,x_scale=[100,57,5],max_nfev=150)
    model=full(best.x,True);error=abs(model[:,1,0]-n.s[:,1,0])
    j=best.jac*np.array([100,57,5])[None,:];sv=np.linalg.svd(j,compute_uv=False)
    result=dict(source=str(source.resolve()),source_sha256=digest(source),
        free_electrical_parameters=['source_island_ground_C_fF','side_Cg_fF'],wire_removed=True,
        C_source_fF=float(best.x[0]),Cg_fF=float(best.x[1]),delay_ps=float(best.x[2]),
        Cg_error_percent=float(100*(best.x[1]/56.7-1)),active_bounds=best.active_mask.tolist(),
        S21_complex_RMS=float(np.sqrt(np.mean(error**2))),max_S21_error=float(error.max()),
        unfitted_S11_RMS=float(np.sqrt(np.mean(abs(model[:,0,0]-n.s[:,0,0])**2))),
        unfitted_S22_RMS=float(np.sqrt(np.mean(abs(model[:,1,1]-n.s[:,1,1])**2))),
        scaled_Jacobian_condition=float(sv[0]/sv[-1]),
        peaks=features(n.f,n.s),model_peaks=features(n.f,model),
        target_T2_design_accepted=False,
        scope='Capacitor-only control. Source C is not the final IDC-equipped T2 C0. Transfer into the grounded-wire network must be checked independently.')
    (out/'fit.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez(out/'response.npz',f=n.f,actual=n.s,model=model)
    fig,ax=plt.subplots(figsize=(12,7));ax.plot(n.f/1e9,abs(n.s[:,1,0]),lw=2.6,color='#b22222',label='Sonnet: wire removed')
    ax.plot(n.f/1e9,abs(model[:,1,0]),'k--',lw=2.2,label='Fixed T1 + capacitor control fit')
    ax.text(.03,.96,f"Capacitor-only control\nCg = {best.x[1]:.3f} fF ({result['Cg_error_percent']:+.2f}%)\nSource parasitic C = {best.x[0]:.2f} fF\nComplex S21 RMS = {result['S21_complex_RMS']:.4f}",transform=ax.transAxes,va='top',fontsize=13,bbox=dict(boxstyle='round',fc='white',ec='gray'))
    ax.set_ylim(0,1.6);ax.set_title('Plain Taper 2: Side-Capacitor Control',fontsize=20,pad=14)
    ax.set_xlabel('Frequency (GHz)',fontsize=17);ax.set_ylabel(r'$|S_{21}|$',fontsize=20)
    ax.tick_params(labelsize=15);ax.legend(fontsize=12,loc='upper right');ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(out/'s21_control_linear.png',dpi=210);plt.close(fig)
    print(json.dumps(result,indent=2),flush=True);return result
