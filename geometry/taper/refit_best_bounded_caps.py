"""Fit saved best-run S21: three T2 Cs free, older Cs within a chosen bound.

No geometry generation or Sonnet calculation. Run as a module from the root.
All inductances stay fixed; one matched fixture delay remains a nuisance term.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[key]='1'
from pathlib import Path
from hashlib import sha256
import json
import numpy as np
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from circuit_common import read_snp
from geometry.taper.progressive_fitting import transmission

ROOT=Path(__file__).resolve().parents[2]
CASE='idc_rail8_right30_dl256_pad148_w148_p200_i42n70_j50_r74tg2'
RESPONSE=ROOT/'output_data/taper_exploration/three_tapers_v3_T2_plain_2026-10-28'/CASE/'response'
BOUND_PERCENT=2.0  # Edit here; every bound gets its own output folder.
NAMES=['T2 C0','T1 C0','Regular C0','Output-end C0',
       'T2 input Cg','T2-to-T1 Cg','T1-to-regular Cg','Regular-to-output Cg']
UNRESTRICTED=[0,4,5]
BOUNDED=[1,2,3,6,7]


def main(bound_percent=BOUND_PERCENT):
    if not 0<bound_percent<100:raise ValueError('Bound must be between 0 and 100 percent')
    fraction=bound_percent/100
    OUT=RESPONSE/f'extraction_three_caps_older_pm{bound_percent:g}pct'
    previous=json.loads((RESPONSE/'extraction_three_caps/fit.json').read_text())
    old=json.loads((RESPONSE/'extraction/fit.json').read_text())
    provenance=previous['provenance'];source=Path(provenance['source'])
    source_hash=sha256(source.read_bytes()).hexdigest()
    if source_hash!=provenance['sha256']:raise ValueError('Original response changed')
    net=read_snp(source);f=net.f;observed=net.s[:,1,0];ls=previous['L_nH']
    reference=np.array(provenance['initial_C_fF']+provenance['initial_G_fF'])
    lower=np.r_[np.zeros(4),np.full(4,.1),-100.]
    upper=np.r_[np.full(4,1000.),np.full(4,3000.),100.]
    lower[BOUNDED]=reference[BOUNDED]*(1-fraction);upper[BOUNDED]=reference[BOUNDED]*(1+fraction)
    scale=np.r_[np.maximum(reference,8),20.]
    lo=lower/scale;hi=upper/scale
    ix=np.unique(np.linspace(0,len(f)-1,2001).astype(int));full=np.arange(len(f))
    def model(z,indices):
        x=z*scale
        return transmission(f[indices],ls,x[:4],x[4:8],x[8])
    def residual(z,indices=ix,smooth=False):
        m=model(z,indices);actual=observed[indices]
        d=(1/m-1/actual)*np.maximum(abs(actual),.03)**2 if smooth else m-actual
        return np.r_[d.real,d.imag]
    starts=[]
    for q in (previous,old,json.loads((ROOT/'investigations/t2_j50_sensitivity/all_capacitances/fit.json').read_text())):
        x=np.array(q['C_fF']+q['G_fF']+[q['delay_ps']])
        starts.append(np.clip(x/scale,lo+1e-8,hi-1e-8))
    five_percent=RESPONSE/'extraction_three_caps_older_pm5pct/fit.json'
    if five_percent.exists():
        q=json.loads(five_percent.read_text())
        x=np.array(q['C_fF']+q['G_fF']+[q['delay_ps']])
        starts.append(np.clip(x/scale,lo+1e-8,hi-1e-8))
    rng=np.random.default_rng(529001)
    for _ in range(12):
        x=reference.copy();x[UNRESTRICTED]*=np.exp(rng.normal(0,.6,3))
        x[BOUNDED]*=rng.uniform(1-.98*fraction,1+.98*fraction,len(BOUNDED))
        starts.append(np.clip(np.r_[x,previous['delay_ps']]/scale,lo+1e-8,hi-1e-8))
    smooth=[least_squares(lambda z:residual(z,smooth=True),z,bounds=(lo,hi),max_nfev=300) for z in starts]
    seeds=sorted(starts+[q.x for q in smooth],key=lambda z:np.linalg.norm(residual(z)))[:8]
    fits=[least_squares(residual,z,bounds=(lo,hi),max_nfev=700,ftol=1e-12,xtol=1e-12,gtol=1e-12) for z in seeds]
    fits.sort(key=lambda q:np.linalg.norm(q.fun))
    final=[least_squares(lambda z:residual(z,full),q.x,bounds=(lo,hi),max_nfev=700,
                         ftol=1e-12,xtol=1e-12,gtol=1e-12) for q in fits[:3]]
    best=min(final,key=lambda q:np.linalg.norm(q.fun));x=best.x*scale;fitted=model(best.x,full)
    rms=float(np.sqrt(np.mean(abs(fitted-observed)**2)))
    if rms>previous['complex_RMS']+1e-8:raise ValueError('Expanded fit should not be worse')
    if not np.all(x[BOUNDED]>=lower[BOUNDED]-1e-7) or not np.all(x[BOUNDED]<=upper[BOUNDED]+1e-7):
        raise ValueError('A constrained capacitance exceeded its allowed range')
    rows=[dict(name=NAMES[i],value_fF=float(x[i]),reference_fF=float(reference[i]),
               change_percent=float(100*(x[i]/reference[i]-1)),
               constraint='free' if i in UNRESTRICTED else f'+/-{bound_percent:g}% of previously fixed value',
               lower_fF=float(lower[i]),upper_fF=float(upper[i]),active_bound=int(best.active_mask[i]),
               near_lower_limit=bool(abs(x[i]-lower[i])<max(reference[i],1)*1e-7),
               near_upper_limit=bool(abs(x[i]-upper[i])<max(reference[i],1)*1e-7)) for i in range(8)]
    result=dict(provenance=dict(source=str(source),sha256=source_hash,reference_fit=str(RESPONSE/'extraction_three_caps/fit.json'),
        code_sha256=sha256(Path(__file__).read_bytes()).hexdigest()),
        L_nH=ls,C_fF=x[:4].tolist(),G_fF=x[4:8].tolist(),delay_ps=float(x[8]),
        parameters=rows,older_capacitance_bound_percent=bound_percent,
        complex_RMS=rms,previous_three_cap_RMS=previous['complex_RMS'],
        previous_five_cap_RMS=old['complex_RMS'],complex_max=float(np.max(abs(fitted-observed))),
        optimizer_success=bool(best.success),active_bounds=best.active_mask.tolist(),
        multistart_full_solutions=[dict(values=(q.x*scale).tolist(),RMS=float(np.sqrt(np.mean(abs(model(q.x,full)-observed)**2)))) for q in final],
        samples=len(f),frequency_range_GHz=[float(f[0]/1e9),float(f[-1]/1e9)],
        metric='RMS of complex S21 residual; real and imaginary jointly fitted',
        no_new_EM=True,capacitance_uniqueness_tested=False,
        scope='Fit to unchanged saved Sonnet data, not an improved geometry or a new target-curve result.')
    OUT.mkdir(exist_ok=True);(OUT/'fit.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez(OUT/'response.npz',f=f,actual=observed,fit=fitted)
    fig,ax=plt.subplots(figsize=(14,8))
    ax.plot(f/1e9,abs(observed),color='#2768aa',lw=2.8,label='Original Sonnet data (2 um mesh)')
    ax.plot(f/1e9,abs(fitted),'k--',lw=2.3,label='Bounded complex circuit fit')
    text='\n'.join(f'{r["name"]} = {r["value_fF"]:.3f} fF'+
        (f' ({r["change_percent"]:+.2f}%)' if i in BOUNDED else ' (free)') for i,r in enumerate(rows))
    text+=f'\nComplex fit RMS = {rms:.6f}\nPrevious three-C fit = {previous["complex_RMS"]:.6f}'
    ax.text(.025,.97,text,transform=ax.transAxes,va='top',fontsize=13,
            bbox=dict(boxstyle='round,pad=.55',fc='white',ec='gray'))
    ax.set(xlabel='Frequency (GHz)',ylabel=r'$|S_{21}|$',title=f'T2 Transmission: Older Capacitances Allowed +/-{bound_percent:g}%',
           xlim=(f[0]/1e9,f[-1]/1e9),ylim=(0,max(abs(observed).max(),abs(fitted).max())*1.65))
    ax.xaxis.label.set_size(19);ax.yaxis.label.set_size(22);ax.title.set_size(22);ax.tick_params(labelsize=15)
    ax.grid(alpha=.2);ax.legend(fontsize=12,loc='upper right');fig.tight_layout()
    for extension in ('png','pdf'):fig.savefig(OUT/f's21_fit.{extension}',dpi=190)
    plt.close(fig);print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
