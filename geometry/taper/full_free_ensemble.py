"""Broad full-filter-only multistart fits; rank small-circuit predictions afterwards."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
from datetime import datetime,timezone
import json
import numpy as np
from scipy.optimize import least_squares,minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from circuit_common import read_snp
from geometry.taper.full_filter import ladder
from geometry.taper.progressive_fitting import transmission
ROOT=Path(__file__).resolve().parents[2];OUT=ROOT/'output_data/taper_exploration/full_free_ensemble_2026-11-02'
FULL=ROOT/'output_data/taper_exploration/full_filter/full_filter_10_regular_3_tapers_each_side_T2_cleanup_2026-11-01'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def save(n,v):(OUT/n).write_text(json.dumps(v,indent=2)+'\n',encoding='utf-8')
def rms(z):return float(np.sqrt(np.mean(abs(z)**2)))
def ri(z):return np.r_[z.real,z.imag]
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 sel=read(ROOT/'output_data/taper_exploration/T2_cleanup_2026-11-01/selected_cleanup.json');a=read(sel['fit_path']);net=read_snp(Path(sel['project']).with_suffix('.snp'));sf=net.f;so=net.s[:,1,0]
 raw=(FULL/'em.log').read_text(encoding='utf-8-sig');(OUT/'em_snapshot.log').write_text(raw,encoding='utf-8');rows=[]
 for line in raw.splitlines():
  try:v=list(map(float,line.split()))
  except ValueError:continue
  if len(v)==9 and 4.65<=v[0]<=5.4:rows.append(v)
 ar=np.array(sorted(rows,key=lambda v:v[0]));f=ar[:,0]*1e9;obs=ar[:,3]*np.exp(1j*np.deg2rad(ar[:,4]));base=np.r_[a['C_fF'][:4],a['G_fF']];bl=np.array(a['L_nH'][:4]);rng=np.random.default_rng(116)
 np.savez(OUT/'snapshot.npz',full_f=f,full_s21=obs,small_f=sf,small_s21=so)
 def expand(x):return base*np.exp(x[:9]),bl*np.exp(x[9:13])
 def model(x,ff=f):
  c,l=expand(x);cs=c[:4];gs=c[4:]
  return ladder(ff,list(l[:3])+[l[3]]*10+list(l[2::-1]),list(cs[:3])+[cs[3]]*10+list(cs[2::-1]),list(gs[:4])+[gs[4]]*9+list(gs[3::-1]),x[13])[:,1,0]
 def delay(z,ff,actual):
  grid=np.linspace(-250,250,501);phase=np.exp(-2j*np.pi*ff[:,None]*grid[None,:]*1e-12);errs=np.mean(abs(z[:,None]*phase-actual[:,None])**2,axis=0);d=grid[np.argmin(errs)]
  objective=lambda t:np.mean(abs(z*np.exp(-2j*np.pi*ff*t*1e-12)-actual)**2)
  q=minimize_scalar(objective,bounds=(max(-250,d-1),min(250,d+1)),method='bounded');return float(q.x)
 def record(x,label,limits):
  c,l=expand(x);fp=model(x);sp=transmission(sf,list(l)+[a['L_nH'][4]],list(c[:4])+[a['C_fF'][4]],c[4:],0);sd=delay(sp,sf,so);sp*=np.exp(-2j*np.pi*sf*sd*1e-12)
  return dict(label=label,limits_ratio=limits,full_complex_RMS=rms(fp-obs),full_magnitude_RMS=rms(abs(fp)-abs(obs)),small_complex_RMS=rms(sp-so),small_magnitude_RMS=rms(abs(sp)-abs(so)),C0_fF=c[:4].tolist(),Cg_fF=c[4:].tolist(),L_nH=l.tolist(),L_change_percent=((l/bl-1)*100).tolist(),full_delay_ps=float(x[13]),small_delay_ps=sd,x=x.tolist())
 old=read(ROOT/'output_data/taper_exploration/full_first_extraction_2026-11-02/report.json');warm=[]
 for run in old['runs']:
  warm.extend(sorted(run['starts'],key=lambda q:q['full_RMS'])[:4])
 warm.sort(key=lambda q:q['full_RMS']);parents=[np.array(q['x']) for q in warm[:16]]
 records=[];report=dict(started_utc=datetime.now(timezone.utc).isoformat(),full_points=len(f),anchor=sel['fit_path'],policy='Only full-filter complex S21 fits component values. All nine capacitances and four shared inductances vary, positive with broad numerical ratio limits. Small-circuit prediction is evaluated after fitting and used to rank an ensemble; it is therefore a selection dataset, not independent validation. Common component values are never refitted to the small circuit. Only its phase delay is adjusted. The small-only output fixture C0 and L remain at saved values because absent from full circuit. No EM, calibration, selection, or geometry changes.',records=records)
 for lower,upper in [(.25,4),(.1,10)]:
  lo=np.r_[np.full(13,np.log(lower)),-250];hi=np.r_[np.full(13,np.log(upper)),250];scale=np.r_[np.ones(13),30]
  starts=[np.clip(x,lo+1e-9,hi-1e-9) for x in parents[:12]]
  # Preserve approximate LC resonance locations in broad random starts.
  for j in range(32):
   x=parents[j%len(parents)].copy();delta=rng.normal(0,.10 if j<16 else .5,4);x[9:13]+=delta;x[:4]-=delta;x[4:9]+=rng.normal(0,.12 if j<16 else .5,5)
   x=np.clip(x,lo+1e-9,hi-1e-9);z=model(x)*np.exp(2j*np.pi*f*x[13]*1e-12);x[13]=delay(z,f,obs);starts.append(x)
  fitted=[]
  for j,seed in enumerate(starts):
   save('status.json',dict(state='full_only_broad_multistart',ratio_limits=[lower,upper],start=j+1,total=len(starts),pid=os.getpid()))
   # Do not overwrite warm starts with magnitude optimization.
   if j>=12 and j%3==0:
    mag=least_squares(lambda x:abs(model(x))-abs(obs),seed,bounds=(lo,hi),x_scale=scale,max_nfev=250);seed=mag.x;seed[13]=delay(model(seed)*np.exp(2j*np.pi*f*seed[13]*1e-12),f,obs)
   q=least_squares(lambda x:ri(model(x)-obs),seed,bounds=(lo,hi),x_scale=scale,max_nfev=1000,ftol=1e-11,xtol=1e-11,gtol=1e-11);fitted.append(q)
   rec=record(q.x,f'Full-only {lower}x to {upper}x, start {j+1}',[lower,upper]);rec.update(success=bool(q.success),nfev=int(q.nfev),active_bounds=q.active_mask.tolist());records.append(rec)
   save('report.json',report)
   if (j+1)%8==0:print('PROGRESS',lower,upper,j+1,'best full',min(z['full_complex_RMS'] for z in records),'best small',min(z['small_complex_RMS'] for z in records),flush=True)
  fitted.sort(key=lambda q:q.fun@q.fun);parents=[q.x for q in fitted[:12]]
 # Keep full-fit best and a separately selected transferable candidate.
 best=min(records,key=lambda q:q['full_complex_RMS']);strict=[r for r in records if r['full_complex_RMS']<.025];eligible=strict or [r for r in records if r['full_complex_RMS']<=best['full_complex_RMS']+.005];chosen=min(eligible,key=lambda q:q['small_complex_RMS'])
 report.update(best_full=best,best_small_among_eligible=chosen,selection_rule='Full RMS below 0.025 if any; otherwise within 0.005 absolute RMS of best full fit, explicitly not called a passing full fit.',number_full_below_0025=len(strict),number_both_below_0025=sum(r['full_complex_RMS']<.025 and r['small_complex_RMS']<.025 for r in records),number_both_below_0005=sum(r['full_complex_RMS']<.005 and r['small_complex_RMS']<.005 for r in records))
 dense=np.linspace(f.min(),f.max(),8001);fig,axes=plt.subplots(3,1,figsize=(14,14))
 axes[0].plot(f/1e9,abs(obs),'o',label=f'Full direct Sonnet ({len(f)} points)',color='black')
 axes[1].plot(sf/1e9,abs(so),label='Small Sonnet, not fitted',color='black')
 for r,style in [(best,'-'),(chosen,'--')]:
  x=np.array(r['x']);c,l=expand(x);label='Best full' if r is best else 'Best small among eligible full fits';p=transmission(sf,list(l)+[a['L_nH'][4]],list(c[:4])+[a['C_fF'][4]],c[4:],r['small_delay_ps'])
  axes[0].plot(dense/1e9,abs(model(x,dense)),style,label=f"{label}: full RMS {r['full_complex_RMS']:.4f}")
  axes[1].plot(sf/1e9,abs(p),style,label=f"{label}: small RMS {r['small_complex_RMS']:.4f}")
 for ax in axes[:2]:ax.set(xlabel='Frequency (GHz)',ylabel='|S21|',ylim=(-.02,1.08));ax.grid(alpha=.2);ax.legend(fontsize=10)
 axes[0].set_title('Broad full-filter-only parameter ensemble; curves between direct points are predictions')
 axes[1].set_title('Transfer to small circuit without refitting common component values')
 axes[2].scatter([q['full_complex_RMS'] for q in records],[q['small_complex_RMS'] for q in records],s=20,alpha=.7);axes[2].axvline(.025,color='gray',ls=':');axes[2].axhline(.025,color='gray',ls=':');axes[2].set(xlabel='Full-filter complex RMS',ylabel='Small-circuit prediction complex RMS',xscale='log',yscale='log');axes[2].grid(alpha=.2)
 fig.tight_layout();fig.savefig(OUT/'ensemble_and_transfer.png',dpi=165);plt.close(fig);save('report.json',report);save('status.json',dict(state='complete',report=str(OUT/'report.json'),no_EM_actions=True));print('COMPLETE',best['full_complex_RMS'],best['small_complex_RMS'],chosen['full_complex_RMS'],chosen['small_complex_RMS'],flush=True)
if __name__=='__main__':main()
