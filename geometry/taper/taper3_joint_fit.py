"""Shared-component fit to independently simulated direct and IDC loading.

All L fixed. Older capacitances move at most +/-1% of selected T2 model.
No EM launched. Saves both residuals, bounds and source hashes; never selects.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import numpy as np
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from circuit_common import read_snp
from geometry.taper import taper3 as work
from geometry.taper.progressive_fitting import transmission
from geometry.taper.staged_fitting import circuit


def run(case, older_fraction=.01):
    out=work.OUT/case/f'joint_older{100*older_fraction:g}percent'
    out.mkdir(exist_ok=True)
    paths=[work.NATIVE/case/s/'generated_ncell_4.snp' for s in ('direct','response')]
    nets=[read_snp(p) for p in paths]
    refs=work.REFERENCE_FIT
    ls=[json.loads((work.OUT/case/'inductance.json').read_text())['L0_nH']]+refs['L_nH']
    old=np.array(refs['C_fF']+refs['G_fF'][1:])
    # Three new Cs, seven bounded older Cs, two independent fixture delays.
    lo=np.r_[0,.1,.1,old*(1-older_fraction),-100,-100]
    hi=np.r_[1000,3000,3000,old*(1+older_fraction),100,100]
    scale=np.r_[251,261,153,old,20,20]
    def expand(x,stage):
        return [x[0]]+x[3:7].tolist(),[None if stage==0 else x[1],x[2]]+x[7:10].tolist()
    ix=[np.unique(np.linspace(0,len(n.f)-1,1501).astype(int)) for n in nets]
    def residual(x,full=False,smooth=False):
        blocks=[]
        for i,n in enumerate(nets):
            indices=np.arange(len(n.f)) if full else ix[i]
            c,g=expand(x,i)
            m=transmission(n.f[indices],ls,c,g,x[10+i]);actual=n.s[indices,1,0]
            d=((1/m-1/actual)*np.maximum(abs(actual),.03)**2 if smooth else m-actual)
            blocks.extend([d.real,d.imag])
        return np.concatenate(blocks)
    starts=[]
    for c in (80,176,251,500,900):
        for gin in (100,261,500):
            for gj in (100,153,300):
                x=np.r_[c,gin,gj,old,5,5].astype(float)
                for i,n in enumerate(nets):
                    peak=np.argmax(abs(n.s[:,1,0]));cp,gp=expand(x,i)
                    zero=transmission(n.f[peak:peak+1],ls,cp,gp,0)[0]
                    x[10+i]=-np.angle(n.s[peak,1,0]/zero)/(2*np.pi*n.f[peak])*1e12
                starts.append(np.clip(x,lo+1e-8,hi-1e-8))
    # Keep low-residual widely separated parameter starts rather than only
    # nearby seeds, because the direct fixture has a competing high-C branch.
    starts.sort(key=lambda x:np.dot(residual(x),residual(x)))
    chosen=[]
    for x in starts:
        if not chosen or all(np.linalg.norm((x[:3]-q[:3])/scale[:3])>.25 for q in chosen):chosen.append(x)
        if len(chosen)==12:break
    smooth=[least_squares(lambda x:residual(x,smooth=True),x,bounds=(lo,hi),x_scale=scale,max_nfev=500) for x in chosen]
    seeds=chosen+[q.x for q in smooth];seeds.sort(key=lambda x:residual(x)@residual(x))
    candidates=[least_squares(residual,x,bounds=(lo,hi),x_scale=scale,max_nfev=600,
        ftol=1e-11,xtol=1e-11,gtol=1e-11) for x in seeds[:8]]
    best=min(candidates,key=lambda q:q.fun@q.fun)
    best=least_squares(lambda x:residual(x,True),best.x,bounds=(lo,hi),
        x_scale=scale,max_nfev=600,ftol=1e-11,xtol=1e-11,gtol=1e-11)
    report=dict(case=case,method='Shared Cs across direct and IDC-connected EM; separate matched fixture delays',
        source_sha256={str(p):work.digest(p) for p in paths},L_nH=ls,
        C_fF=expand(best.x,1)[0],G_fF=expand(best.x,1)[1],delay_ps=best.x[10:].tolist(),
        older_relative_bounds=older_fraction,older_relative_changes=(best.x[3:10]/old-1).tolist(),
        active_bounds=best.active_mask.tolist(),optimizer_success=bool(best.success),
        new_C_fF=best.x[:3].tolist(),targets_fF=[251,261,153],
        new_C_errors_percent=(100*(best.x[:3]/[251,261,153]-1)).tolist(),stages={})
    for i,n in enumerate(nets):
        name=('direct','response')[i];c,g=expand(best.x,i)
        m=circuit(n.f,ls,c,g,left_ps=best.x[10+i]/2,right_ps=best.x[10+i]/2)
        delta=abs(m[:,1,0]-n.s[:,1,0]);rms=float(np.sqrt(np.mean(delta**2)))
        report['stages'][name]=dict(complex_RMS=rms,complex_max=float(delta.max()))
        np.savez(out/(name+'.npz'),f=n.f,actual=n.s,model=m)
        fig,ax=plt.subplots(figsize=(13,7));ax.plot(n.f/1e9,abs(n.s[:,1,0]),lw=2.5,label='Sonnet')
        ax.plot(n.f/1e9,abs(m[:,1,0]),'k--',lw=2,label='Shared-component fit')
        ax.text(.025,.96,f'C0 T3 = {best.x[0]:.3f} fF\nInput Cg = {best.x[1]:.3f} fF\nT3-T2 Cg = {best.x[2]:.3f} fF\nFit RMS = {rms:.5f}',
            transform=ax.transAxes,va='top',fontsize=14,bbox=dict(fc='white',ec='gray',boxstyle='round'))
        ax.set(xlabel='Frequency (GHz)',ylabel=r'$|S_{21}|$',title='T3 Shared-Component Check: '+name,
            ylim=(0,max(abs(n.s[:,1,0]).max(),abs(m[:,1,0]).max())*1.6))
        ax.tick_params(labelsize=14);ax.xaxis.label.set_size(18);ax.yaxis.label.set_size(20)
        ax.grid(alpha=.2);ax.legend();fig.tight_layout();fig.savefig(out/(name+'.png'),dpi=180);plt.close(fig)
    report['accepted']=False
    report['scope']='Diagnostic only. Bounds and stage consistency must be reviewed before accepting extracted values.'
    (out/'fit.json').write_text(json.dumps(report,indent=2)+'\n')
    print('JOINT RESULT',json.dumps(report),flush=True)
    return report


if __name__=='__main__':run('w84_dr256_rightpad148_i36_j22')
