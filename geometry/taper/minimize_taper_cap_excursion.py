"""Reduce fractional changes of the passing four-cap circuit prediction."""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import json
import numpy as np
from scipy.optimize import minimize,brentq,minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from geometry.taper.minimal_cap_changes import OUT
from geometry.taper.progressive_fitting import transmission
from geometry.specifications import TAPER_2,TAPER_1,UNIT_CELL
from circuit_common import read_snp
from pathlib import Path


def main():
    source=json.loads((OUT/'results.json').read_text());best=source['best']
    fit=json.loads(Path(source['source_fit']).read_text());net=read_snp(Path(fit['provenance']['source']));f=net.f
    base=np.array(source['baseline_C_fF']+source['baseline_G_fF']);free=best['indices']
    ls=source['fixed_L_nH'];delay=fit['delay_ps'];limit=.0245
    ts=[TAPER_2,TAPER_1,UNIT_CELL,UNIT_CELL]
    target=transmission(f,[t.L0_nH for t in ts],[t.C0_fF for t in ts],[t.Cg_input_fF for t in ts],0)
    full=np.arange(len(f));ix=np.unique(np.linspace(0,len(f)-1,2001).astype(int))
    lower=np.r_[np.log(np.array(source['cap_lower_bounds_fF'])[free]/base[free]),-5.]
    upper=np.r_[np.log(np.array(source['cap_upper_bounds_fF'])[free]/base[free]),5.]
    def expand(z):
        cp=base.copy();cp[free]=base[free]*np.exp(z[:-1]);return cp
    def model(z,indices):
        cp=expand(z);return transmission(f[indices],ls,cp[:4],cp[4:],delay)
    def mse(z,indices):
        d=model(z,indices)-target[indices]*np.exp(-2j*np.pi*f[indices]*z[-1]*20e-12)
        return float(np.mean(abs(d)**2))
    endpoint=np.r_[np.log(np.array(best['new_fF'])/base[free]),best['target_delay_ps']/20]
    objective=lambda z:float(z[:-1]@z[:-1])
    # Keep a verified feasible point even if the constrained optimizer leaves
    # the narrow feasible branch. First find a smaller feasible interpolation.
    def path_point(v):
        z=endpoint.copy();z[:-1]*=v
        grid=np.linspace(-5,5,41)
        values=[mse(np.r_[z[:-1],t],full) for t in grid];k=int(np.argmin(values))
        q=minimize_scalar(lambda t:mse(np.r_[z[:-1],t],full),
            bounds=(max(-5,grid[k]-.3),min(5,grid[k]+.3)),method='bounded')
        z[-1]=q.x;return z
    grid=np.linspace(0,1,61);path=[path_point(v) for v in grid]
    k=next(i for i,z in enumerate(path) if mse(z,full)<=limit**2)
    root=brentq(lambda v:mse(path_point(v),full)-limit**2,grid[k-1],grid[k],xtol=1e-10)
    candidates=[endpoint,path_point(root)]
    q=minimize(objective,candidates[-1],bounds=list(zip(lower,upper)),
        constraints=[dict(type='ineq',fun=lambda z:(limit**2-mse(z,full))*1e4)],
        method='SLSQP',options=dict(maxiter=1000,ftol=1e-11,eps=1e-7))
    if np.sqrt(mse(q.x,full))<.025:candidates.append(q.x)
    chosen=min(candidates,key=objective)
    cp=expand(chosen);rms=float(np.sqrt(mse(chosen,full)))
    if rms>.025:raise ValueError('Candidate did not meet actual full-grid criterion')
    fixed=[i for i in range(8) if i not in free]
    if not np.array_equal(cp[fixed],base[fixed]):raise ValueError('An unselected capacitor changed')
    result=dict(indices=free,names=best['names'],old_fF=base[free].tolist(),new_fF=cp[free].tolist(),
        changes_percent=(100*(cp[free]/base[free]-1)).tolist(),C_fF=cp[:4].tolist(),G_fF=cp[4:].tolist(),
        fixed_L_nH=ls,complex_RMS=rms,target_delay_ps=float(chosen[-1]*20),
        objective='Reduce squared log fractional changes within the located feasible branch; verified path fallback. Not a certified global minimum.',
        optimizer_success=bool(q.success),message=q.message,scope='Circuit-only prediction; no EM or geometry changes',
        original_target_RMS=source['baseline_model_target_RMS'],criterion=.025)
    (OUT/'smaller_four_cap_changes.json').write_text(json.dumps(result,indent=2)+'\n')
    pred=model(chosen,full);current=transmission(f,ls,base[:4],base[4:],delay)
    fig,ax=plt.subplots(figsize=(14,8))
    ax.plot(f/1e9,abs(target),'k',lw=2.5,label='Specified target circuit')
    ax.plot(f/1e9,abs(current),color='#888',lw=1.8,label='Current extracted circuit')
    ax.plot(f/1e9,abs(pred),'--',color='#087f65',lw=2.5,label='Four-capacitance prediction')
    text='\n'.join(f'{n}: {a:.3f} -> {b:.3f} fF ({p:+.2f}%)' for n,a,b,p in zip(result['names'],result['old_fF'],result['new_fF'],result['changes_percent']))
    text+=f'\nPredicted target RMS: {rms:.5f}\nRegular/output values, T1 C0 and all L fixed'
    ax.text(.025,.97,text,transform=ax.transAxes,va='top',fontsize=13,bbox=dict(boxstyle='round,pad=.5',fc='white',ec='gray'))
    ax.set(xlim=(4.65,5.4),ylim=(0,max(abs(target).max(),abs(pred).max(),abs(current).max())*1.5),
           xlabel='Frequency (GHz)',ylabel=r'$|S_{21}|$',title='Circuit Prediction: Smaller Changes to Four Taper Capacitances')
    ax.xaxis.label.set_size(20);ax.yaxis.label.set_size(23);ax.title.set_size(21);ax.tick_params(labelsize=16)
    ax.grid(alpha=.2);ax.legend(fontsize=12,loc='upper right');fig.tight_layout()
    for ext in ('png','pdf'):fig.savefig(OUT/f'smaller_four_cap_changes.{ext}',dpi=190)
    plt.close(fig);print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
