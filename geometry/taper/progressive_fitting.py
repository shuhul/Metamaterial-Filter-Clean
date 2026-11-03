"""Complex-S21 staged fits with explicitly named free capacitances and fixed L.

Geometry remains independent of fitting. Saved fits always record which values
were free, bound activity, independent reflection checks and parameter changes.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import numpy as np
from scipy.optimize import least_squares
from circuit_common import read_snp
from geometry.taper.staged_fitting import circuit
from geometry.taper.run_bridged_stage import digest
from geometry.specifications import S21_COMPLEX_RMS_CUTOFF


def transmission(f,ls,cs,gs,delay):
    w=2*np.pi*np.asarray(f);a=np.ones(len(w),complex);b=np.zeros(len(w),complex);c=b.copy();d=a.copy()
    for l,cap,gap in zip(ls,cs,gs):
        if gap is not None:
            z=1/(1j*w*gap*1e-15);b=b+a*z;d=d+c*z
        y=1j*w*cap*1e-15+1/(1j*w*l*1e-9);a=a+b*y;c=c+d*y
    return 2/(a+b/50+50*c+d)*np.exp(-1j*w*delay*1e-12)


def fit(snp,out,ls,cs,gs,free,*,title='Progressive Taper Extraction',label='Sonnet',max_nfev=400,
        capacitance_seeds=None, capacitance_bounds=None):
    """free=[('C',index),('G',index)]; total matched phase delay is also free."""
    out.mkdir(parents=True,exist_ok=True);net=read_snp(snp);f=net.f;obs=net.s[:,1,0]
    if len(f)<101:
        raise ValueError('At least 101 frequency samples required; sparse screening points cannot constrain this fit')
    provenance=dict(source=str(snp),sha256=digest(snp),code_sha256=digest(__file__),
                    fixed_L_nH=ls,initial_C_fF=cs,initial_G_fF=gs,free_capacitances=free,
                    delay_bounds_ps=[-100,100])
    if capacitance_seeds is not None:
        provenance['additional_capacitance_seeds']=capacitance_seeds
    if capacitance_bounds is not None:
        provenance['capacitance_bounds_fF']=capacitance_bounds
    marker=out/'fit.json'
    if marker.exists():
        old=json.loads(marker.read_text())
        if old['provenance']==json.loads(json.dumps(provenance)):return old
    initial=[(cs if k=='C' else gs)[i] for k,i in free]+[5.]
    scale=[max(v,8) for v in initial[:-1]]+[20.]
    lower=[0. if k=='C' else .1 for k,i in free]+[-100.]
    upper=[1000. if k=='C' else 3000. for k,i in free]+[100.]
    if capacitance_bounds is not None:
        if len(capacitance_bounds)!=len(free):
            raise ValueError('One bound pair required per free capacitance')
        for j,(lo,hi) in enumerate(capacitance_bounds):
            if not lower[j]<=lo<hi<=upper[j]:
                raise ValueError('Invalid capacitance bounds')
            lower[j],upper[j]=lo,hi
    def expand(x):
        cp=list(cs);gp=list(gs)
        for (k,i),v in zip(free,x): (cp if k=='C' else gp)[i]=v
        return cp,gp
    ix=np.unique(np.linspace(0,len(f)-1,min(len(f),2001)).astype(int))
    def residual(x,indices=ix,smooth=False):
        cp,gp=expand(x);m=transmission(f[indices],ls,cp,gp,x[-1]);actual=obs[indices]
        delta=(1/m-1/actual)*np.maximum(abs(actual),.03)**2 if smooth else m-actual
        return np.r_[delta.real,delta.imag]
    rng=np.random.default_rng(740)
    peak=int(np.argmax(abs(obs)))
    def seed_phase(seed):
        cp,gp=expand(seed)
        zero=transmission(f[peak:peak+1],ls,cp,gp,0)[0]
        seed[-1]=-np.angle(obs[peak]/zero)/(2*np.pi*f[peak])*1e12
        return np.clip(seed,np.array(lower)+1e-8,np.array(upper)-1e-8)
    starts=[seed_phase(np.array(initial))]
    for values in capacitance_seeds or []:
        if len(values)!=len(free):raise ValueError('A seed must specify every free capacitance')
        starts.append(seed_phase(np.array(list(values)+[5.])))
    for _ in range(7):
        seed=np.array(initial);seed[:-1]*=np.exp(rng.normal(0,.16,len(free)))
        starts.append(seed_phase(seed))
    smooth=[least_squares(lambda x:residual(x,smooth=True),x,bounds=(lower,upper),
                          x_scale=scale,max_nfev=max_nfev) for x in starts]
    seeds=sorted(starts+[q.x for q in smooth],key=lambda x:np.linalg.norm(residual(x)))[:5]
    fits=[least_squares(residual,x,bounds=(lower,upper),x_scale=scale,max_nfev=max_nfev,
                        ftol=1e-11,xtol=1e-11,gtol=1e-11) for x in seeds]
    best=min(fits,key=lambda q:q.fun@q.fun)
    full=np.arange(len(f));best=least_squares(lambda x:residual(x,full),best.x,bounds=(lower,upper),
                 x_scale=scale,max_nfev=max_nfev,ftol=1e-11,xtol=1e-11,gtol=1e-11)
    cp,gp=expand(best.x);model=circuit(f,ls,cp,gp,left_ps=best.x[-1]/2,right_ps=best.x[-1]/2)
    error=abs(model[:,1,0]-obs)
    r=dict(provenance=provenance,L_nH=ls,C_fF=list(map(float,cp)),G_fF=[None if q is None else float(q) for q in gp],
           delay_ps=float(best.x[-1]),free_capacitances=free,optimizer_success=bool(best.success),
           frequency_samples=len(f),frequency_range_GHz=[float(f[0]/1e9),float(f[-1]/1e9)],
           active_bounds=best.active_mask.tolist(),complex_RMS=float(np.sqrt(np.mean(error**2))),
           complex_max=float(error.max()),
           S11_magnitude_RMS=float(np.sqrt(np.mean((abs(model[:,0,0])-abs(net.s[:,0,0]))**2))),
           S22_magnitude_RMS=float(np.sqrt(np.mean((abs(model[:,1,1])-abs(net.s[:,1,1]))**2))),
           capacitances_accepted=False,scope='Fixed-inductance complex S21 fit; reflection magnitudes are not fitted')
    r['curve_acceptance']=dict(metric='complex_S21_RMS',maximum=S21_COMPLEX_RMS_CUTOFF,
                               policy='User-approved cutoff; maximum error and peak errors are diagnostics.')
    r['curve_pass']=bool(r['optimizer_success'] and r['complex_RMS']<=S21_COMPLEX_RMS_CUTOFF)
    r['interior_free_parameters']=not any(r['active_bounds'])
    marker.write_text(json.dumps(r,indent=2)+'\n')
    np.savez(out/'response.npz',f=f,actual=net.s,model=model)
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(12,7))
    ax.plot(f/1e9,abs(obs),color='#b22222',lw=2.8,label=label)
    ax.plot(f/1e9,abs(model[:,1,0]),'k--',lw=2.2,label='Complex circuit fit')
    text='\n'.join([f"C[{i}] = {cp[i]:.3f} fF" if k=='C' else f"Cg[{i}] = {gp[i]:.3f} fF" for k,i in free])
    ax.text(.03,.97,text+f"\nComplex S21 RMS = {r['complex_RMS']:.5f}\nFixed L = {', '.join(f'{v:.4f}' for v in ls)} nH",
            transform=ax.transAxes,va='top',fontsize=14,bbox=dict(boxstyle='round,pad=.5',fc='white',ec='gray'))
    ax.set_ylim(0,max(abs(obs).max(),abs(model[:,1,0]).max())*1.65)
    ax.set(title=title);ax.title.set_size(21);ax.set_xlabel('Frequency (GHz)',fontsize=18);ax.set_ylabel(r'$|S_{21}|$',fontsize=22)
    ax.tick_params(labelsize=15);ax.legend(fontsize=12,loc='upper right');ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(out/'s21_fit.png',dpi=190);plt.close(fig)
    print('FIT',str(out),r['complex_RMS'],r['complex_max'],r['C_fF'],r['G_fF'],flush=True)
    return r
