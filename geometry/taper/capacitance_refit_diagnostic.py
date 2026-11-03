"""Cached-data diagnostic: zero T2 shunt C, fixed L, constrained vs free C.

The all-C fit is exploratory. Its parameters do not replace calibrated values.
Run: python -m geometry.taper.capacitance_refit_diagnostic [meander width]
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import numpy as np
from scipy.optimize import least_squares
from circuit_common import read_snp
from geometry.taper.minimum_inductance_review import OUT,ROOT,PROFILES,assessment
from geometry.taper.staged_fitting import LT1,LR
from geometry.taper.run_bridged_stage import digest

PROFILE='full_T1_staged_calibration'
NAMES=['C0_T1_fF','C0_regular_fF','C0_output_end_fF',
       'Cg_T2_T1_fF','Cg_T1_regular_fF','Cg_output_fF','delay_ps']


def predict(f,l,x):
    """Four shunt branches, T2 C exactly zero; matched total fixture delay."""
    ct,cr,ce,gin,gj,go,delay=x
    w=2*np.pi*np.asarray(f)
    a=np.ones(len(w),complex);b=np.zeros(len(w),complex)
    c=b.copy();d=a.copy()
    for inductance,cap,gap in zip((l,LT1,LR,LR),(0,ct,cr,ce),(None,gin,gj,go)):
        if gap is not None:
            z=1/(1j*w*gap*1e-15)
            b=b+a*z;d=d+c*z
        y=1j*w*cap*1e-15+1/(1j*w*inductance*1e-9)
        a=a+b*y;c=c+d*y
    return 2/(a+b/50+50*c+d)*np.exp(-1j*w*delay*1e-12)


def fit_data(f,actual,l,free_all):
    p=PROFILES[PROFILE]
    base=np.array([p['ct'],p['cr'],p['ce'],56.7,p['gj'],p['go'],5.])
    free=np.arange(7) if free_all else np.array([3,6])
    lower=np.array([1,1,1,5,.2,.2,-20.])[free]
    upper=np.array([500,500,500,180,40,40,30.])[free]
    scales=np.array([100,190,190,57,8,8,5.])[free]
    ix=np.unique(np.linspace(0,len(f)-1,min(len(f),1501)).astype(int))
    def expand(x):
        full=base.copy();full[free]=x;return full
    def residual(x,frequency,obs):
        delta=predict(frequency,l,expand(x))-obs
        return np.r_[delta.real,delta.imag]
    rng=np.random.default_rng(2520)
    seeds=[base[free]]
    if free_all:
        for _ in range(5):
            seed=base.copy();seed[:6]*=np.exp(rng.normal(0,.12,6))
            seeds.append(np.clip(seed[free],lower+1e-8,upper-1e-8))
    solutions=[least_squares(residual,seed,args=(f[ix],actual[ix]),bounds=(lower,upper),
                             x_scale=scales,max_nfev=600,ftol=1e-11,xtol=1e-11,gtol=1e-11)
               for seed in seeds]
    best=min(solutions,key=lambda r:r.fun@r.fun)
    best=least_squares(residual,best.x,args=(f,actual),bounds=(lower,upper),x_scale=scales,
                       max_nfev=400,ftol=1e-11,xtol=1e-11,gtol=1e-11)
    x=expand(best.x);model=predict(f,l,x);delta=abs(model-actual)
    # Dimensionless parameter-scaled sensitivity: large ratio flags weak identifiability.
    sv=np.linalg.svd(best.jac*scales[None,:],compute_uv=False)
    r=dict(mode='all_capacitances_free' if free_all else 'coupling_only',
           fixed_L_nH=dict(T2=l,T1=LT1,regular=LR,output_end=LR),fixed_C0_T2_fF=0,
           parameters=dict(zip(NAMES,map(float,x))),free_parameters=[NAMES[i] for i in free],
           seed_parameters=dict(zip(NAMES,map(float,base))),reference=PROFILE,
           percent_change_from_seed={NAMES[i]:float(100*(x[i]/base[i]-1)) for i in range(6)},
           complex_RMS=float(np.sqrt(np.mean(delta**2))),complex_max=float(delta.max()),
           active_bounds=best.active_mask.tolist(),optimizer_success=bool(best.success),
           scaled_jacobian_singular_values=sv.tolist(),
           scaled_jacobian_condition=float(sv[0]/sv[-1]) if sv[-1]>0 else None,
           starts=[dict(cost=float(s.fun@s.fun),success=bool(s.success),
                        parameters=dict(zip(NAMES,map(float,expand(s.x))))) for s in solutions],
           parameter_extraction_accepted=False,
           caveat='S21-only exploratory refit; good transmission agreement does not establish unique physical capacitances.')
    r['assessment']=assessment(f,actual,model,r)
    return r,model


def run(width=80):
    case=OUT/f'meander_w{width}'
    l=json.loads((case/'Lcal/results.json').read_text())['L0_nH']
    sources=dict(physical=ROOT/f'sonnet_files/taper_exploration/minimum_inductance_2026-10-23/meander_w{width}/direct/generated_ncell_3.snp',
                 ideal=case/'paired_ideal/ideal_loaded.s2p')
    folder=case/'capacitance_refit_zero_T2_C';folder.mkdir(parents=True,exist_ok=True)
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(16,7),sharex=True,sharey=True)
    report={}
    for ax,(kind,source) in zip(axes,sources.items()):
        net=read_snp(source);actual=net.s[:,1,0];rows={}
        ax.plot(net.f/1e9,abs(actual),color='#b22222',lw=2.8,label=f'{kind.capitalize()} response')
        for free_all,color,style in ((False,'#444444','--'),(True,'#1565c0','-.')):
            r,model=fit_data(net.f,actual,l,free_all)
            r['source_sha256']=digest(source);r['code_sha256']=digest(__file__)
            rows[r['mode']]=r
            (folder/f"{kind}_{r['mode']}.json").write_text(json.dumps(r,indent=2)+'\n')
            np.savez(folder/f"{kind}_{r['mode']}.npz",f=net.f,actual=actual,model=model)
            ax.plot(net.f/1e9,abs(model),ls=style,color=color,lw=2,
                    label=f"{r['mode'].replace('_',' ')}: RMS {r['complex_RMS']:.5f}")
            print(kind,r['mode'],r['complex_RMS'],r['parameters'],flush=True)
        report[kind]=rows
        ax.set_title(f'{kind.capitalize()} inductor: {l:.5f} nH',fontsize=19)
        ax.set_xlabel('Frequency (GHz)',fontsize=17);ax.tick_params(labelsize=14)
        ax.legend(fontsize=11);ax.grid(alpha=.25)
    axes[0].set_ylabel(r'$|S_{21}|$',fontsize=22)
    fig.suptitle('Fixed Inductances; T2 Shunt Capacitance = 0',fontsize=23)
    fig.tight_layout();fig.savefig(folder/'s21_comparison.png',dpi=190);plt.close(fig)
    (folder/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    return report


if __name__=='__main__':
    import sys
    run(int(sys.argv[1]) if len(sys.argv)>1 else 80)
