"""Search for the fewest capacitor changes needed to meet target S21 RMS.

Circuit prediction only. All measured L and unselected C remain fixed at the
adopted original-2um extraction. Only matched target fixture delay is nuisance.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
from pathlib import Path
import json,itertools
import numpy as np
from scipy.optimize import least_squares,minimize_scalar,brentq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from geometry.taper.taper2_grid4 import OLD_FIT
from geometry.taper.progressive_fitting import transmission
from geometry.specifications import TAPER_2,TAPER_1,UNIT_CELL
from circuit_common import read_snp

OUT=OLD_FIT.parent.parent/'minimal_taper_capacitance_changes'
ALLOWED=[0,1,4,5,6]  # T2 C0, T1 C0, input/T2 Cg, T2/T1 Cg, T1/regular Cg.
NAMES=['T2 C0','T1 C0','Regular C0','Output-end C0','T2 input Cg',
       'T2-to-T1 Cg','T1-to-regular Cg','Regular-to-output Cg']


def main():
    fit=json.loads((OLD_FIT.parent.parent/'extraction_three_caps_adopted_reference/fit.json').read_text())
    net=read_snp(Path(fit['provenance']['source']));f=net.f
    ls=fit['L_nH'];base=np.array(fit['C_fF']+fit['G_fF']);delay=fit['delay_ps']
    ts=[TAPER_2,TAPER_1,UNIT_CELL,UNIT_CELL]
    target=transmission(f,[t.L0_nH for t in ts],[t.C0_fF for t in ts],[t.Cg_input_fF for t in ts],0)
    ix=np.unique(np.linspace(0,len(f)-1,1501).astype(int));full=np.arange(len(f))
    low=np.r_[np.full(4,1e-5),np.full(4,.1)];high=np.r_[np.full(4,1000.),np.full(4,3000.)]
    def model(c,indices):return transmission(f[indices],ls,c[:4],c[4:],delay)
    def phase_score(c,indices=full):
        s=model(c,indices);t=target[indices];freq=f[indices]
        # Search one period of the central frequency then refine the winner.
        obj=lambda tau:float(np.mean(abs(s-t*np.exp(-2j*np.pi*freq*tau*1e-12))**2))
        grid=np.linspace(-100,100,41);k=int(np.argmin([obj(v) for v in grid]))
        q=minimize_scalar(obj,bounds=(max(-100,grid[k]-6),min(100,grid[k]+6)),method='bounded')
        return float(np.sqrt(q.fun)),float(q.x)
    baseline,tau0=phase_score(base)
    cache={}
    # Reuse completed single/pair calculations only when every starting value
    # and fixed inductance matches. The earlier broader search is preserved.
    prior=OLD_FIT.parent.parent/'minimal_capacitance_changes/results.json'
    if prior.exists():
        saved=json.loads(prior.read_text())
        if (saved['baseline_C_fF']==base[:4].tolist() and saved['baseline_G_fF']==base[4:].tolist()
                and saved['fixed_L_nH']==ls):
            for row in saved['single_changes']+saved['pair_changes']:cache[tuple(row['indices'])]=row
    def solve(free):
        free=list(free);lower=np.r_[np.log(low[free]/base[free]),-5.]
        upper=np.r_[np.log(high[free]/base[free]),5.]
        def expand(z):
            c=base.copy();c[free]=base[free]*np.exp(z[:-1]);return c
        def residual(z,indices=ix):
            d=model(expand(z),indices)-target[indices]*np.exp(-2j*np.pi*f[indices]*z[-1]*20e-12)
            return np.r_[d.real,d.imag]
        seeds=[np.r_[np.zeros(len(free)),tau0/20]]
        if len(free)==1:
            # Wide capacitance scan plus dense local values catches narrow
            # basins and avoids relying on one nonlinear least-squares start.
            grid=np.unique(np.r_[np.linspace(-.35,.35,81),np.linspace(lower[0],upper[0],101)])
            screened=[]
            for v in grid:
                c=base.copy();c[free[0]]=base[free[0]]*np.exp(v)
                rms,tau=phase_score(c,ix);screened.append((rms,v,tau))
            screened.sort();seeds += [np.array([v,tau/20]) for _,v,tau in screened[:6]]
        else:
            rng=np.random.default_rng(9731+sum(free))
            for _ in range(9):seeds.append(np.r_[rng.normal(0,.18,len(free)),tau0/20])
            for subset,row in cache.items():
                if set(subset).issubset(free):
                    cp=np.array(row['C_fF']+row['G_fF'])
                    seeds.append(np.r_[np.log(cp[free]/base[free]),row['target_delay_ps']/20])
        fits=[least_squares(residual,np.clip(z,lower+1e-9,upper-1e-9),bounds=(lower,upper),max_nfev=350,
                           ftol=1e-11,xtol=1e-11,gtol=1e-11) for z in seeds]
        fits.sort(key=lambda q:np.linalg.norm(q.fun))
        fine=[least_squares(lambda z:residual(z,full),q.x,bounds=(lower,upper),max_nfev=450,
              ftol=1e-12,xtol=1e-12,gtol=1e-12) for q in fits[:3]]
        best=min(fine,key=lambda q:np.linalg.norm(q.fun));c=expand(best.x)
        rms,tau=phase_score(c)
        result=dict(indices=free,names=[NAMES[i] for i in free],old_fF=base[free].tolist(),
            new_fF=c[free].tolist(),changes_percent=(100*(c[free]/base[free]-1)).tolist(),
            C_fF=c[:4].tolist(),G_fF=c[4:].tolist(),complex_RMS=rms,target_delay_ps=tau,
            passes=bool(rms<=.025),optimizer_success=bool(best.success))
        if len(free)==1 and rms<.025:
            index=free[0];end=float(np.log(c[index]/base[index]))
            def at_fraction(v):
                cp=base.copy();cp[index]=base[index]*np.exp(end*v);return cp
            # First acceptable point along this direction, with all other Cs
            # fixed. This is the minimum excursion along the located basin.
            grid=np.linspace(0,1,101);values=[phase_score(at_fraction(v))[0]-.025 for v in grid]
            k=next(k for k,v in enumerate(values) if v<=0)
            root=brentq(lambda v:phase_score(at_fraction(v))[0]-.025,grid[k-1],grid[k],xtol=1e-10)
            threshold=at_fraction(root);rr,tt=phase_score(threshold)
            result['threshold_change']=dict(value_fF=float(threshold[index]),
                change_percent=float(100*(threshold[index]/base[index]-1)),complex_RMS=rr,target_delay_ps=tt,
                scope='First threshold crossing from baseline toward this located minimum; negligible safety margin.')
        return result
    OUT.mkdir(exist_ok=True)
    result=dict(baseline_model_target_RMS=baseline,cutoff=.025,fixed_L_nH=ls,
        baseline_C_fF=base[:4].tolist(),baseline_G_fF=base[4:].tolist(),
        source_fit=str(OLD_FIT.parent.parent/'extraction_three_caps_adopted_reference/fit.json'),
        scope='Circuit prediction; no geometry or EM optimization; no claim of unique extraction.',
        cap_lower_bounds_fF=low.tolist(),cap_upper_bounds_fF=high.tolist(),
        allowed_indices=ALLOWED,allowed_names=[NAMES[i] for i in ALLOWED],
        fixed_capacitance_indices=[i for i in range(8) if i not in ALLOWED],
        minimum_capacitance_count=None,search_by_count={})
    available=[]
    for count in range(1,len(ALLOWED)+1):
        group=[]
        for subset in itertools.combinations(ALLOWED,count):
            row=cache.get(tuple(subset))
            if row is None:row=solve(subset);cache[tuple(subset)]=row
            cp=np.array(row['C_fF']+row['G_fF'])
            locked=[i for i in range(8) if i not in subset]
            if not np.array_equal(cp[locked],base[locked]):raise ValueError('An unselected capacitance changed')
            group.append(row);result['search_by_count'][str(count)]=group
            (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
            print(row['names'],row['complex_RMS'],row['new_fF'],flush=True)
        available=group
        if any(row['passes'] for row in group):
            result['minimum_capacitance_count']=count;break
    best=min(available,key=lambda r:r['complex_RMS']);result['best']=best
    (OUT/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    cp=np.array(best['C_fF']+best['G_fF']);pred=model(cp,full)
    fig,ax=plt.subplots(figsize=(14,8))
    ax.plot(f/1e9,abs(target),'k',lw=2.6,label='Specified target circuit')
    ax.plot(f/1e9,abs(model(base,full)),color='#888',lw=1.8,label='Current extracted circuit')
    ax.plot(f/1e9,abs(pred),'--',color='#087f65',lw=2.6,label='Minimal-capacitance circuit change')
    text='\n'.join(f'{name}: {old:.4f} -> {new:.4f} fF ({pct:+.2f}%)' for name,old,new,pct in zip(best['names'],best['old_fF'],best['new_fF'],best['changes_percent']))
    text+=f'\nPredicted target RMS: {best["complex_RMS"]:.6f}\nCurrent model target RMS: {baseline:.6f}\nAll inductances and other capacitances fixed'
    ax.text(.025,.97,text,transform=ax.transAxes,va='top',fontsize=14,bbox=dict(boxstyle='round,pad=.5',fc='white',ec='gray'))
    ax.set(xlabel='Frequency (GHz)',ylabel=r'$|S_{21}|$',title='Circuit Prediction: Fewest Capacitance Changes',
           xlim=(4.65,5.4),ylim=(0,max(abs(pred).max(),abs(target).max())*1.5))
    ax.xaxis.label.set_size(20);ax.yaxis.label.set_size(23);ax.title.set_size(22);ax.tick_params(labelsize=16)
    ax.grid(alpha=.2);ax.legend(fontsize=12,loc='upper right');fig.tight_layout()
    for ext in ('png','pdf'):fig.savefig(OUT/f's21_comparison.{ext}',dpi=190)
    plt.close(fig)
    print('MINIMUM COUNT',result['minimum_capacitance_count'],'BEST',json.dumps(best),flush=True)


if __name__=='__main__':main()
