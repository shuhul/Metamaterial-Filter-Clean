"""Fit completed full-filter direct fields; validate on saved staged response. No EM controls."""
import os
for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[k]='1'
from pathlib import Path
import json,hashlib
from datetime import datetime,timezone
import numpy as np
from scipy.optimize import least_squares,minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from circuit_common import read_snp
from geometry.taper.full_filter import ladder
from geometry.taper.progressive_fitting import transmission
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output_data/taper_exploration/full_first_extraction_2026-11-02'
FULL=ROOT/'output_data/taper_exploration/full_filter/full_filter_10_regular_3_tapers_each_side_T2_cleanup_2026-11-01'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def save(n,x):(OUT/n).write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8')
def rms(z):return float(np.sqrt(np.mean(abs(z)**2)))
def ri(z):return np.r_[z.real,z.imag]
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 chosen=read(ROOT/'output_data/taper_exploration/T2_cleanup_2026-11-01/selected_cleanup.json');anchor=read(chosen['fit_path']);net=read_snp(Path(chosen['project']).with_suffix('.snp'))
 text=(FULL/'em.log').read_text(encoding='utf-8-sig');(OUT/'em_snapshot.log').write_text(text,encoding='utf-8');rows=[]
 for line in text.splitlines():
  try:v=list(map(float,line.split()))
  except ValueError:continue
  if len(v)==9 and 4.65<=v[0]<=5.4:rows.append(v)
 a=np.array(sorted(rows,key=lambda v:v[0]));f=a[:,0]*1e9;obs=a[:,3]*np.exp(1j*np.deg2rad(a[:,4]));sf=net.f;so=net.s[:,1,0]
 np.savez(OUT/'snapshot.npz',full_f=f,full_s21=obs,staged_f=sf,staged_s21=so)
 base=np.r_[anchor['C_fF'][:4],anchor['G_fF']];bl=np.array(anchor['L_nH'][:4]);rng=np.random.default_rng(113)
 # x: 4 C0, 5 Cg, 4 L, full delay ps. All caps/L expressed as fractional log ratios.
 def unpack(x):return base*np.exp(x[:9]),bl*np.exp(x[9:13])
 def model(x,freq=f):
  c,l=unpack(x);cs=c[:4];gs=c[4:]
  return ladder(freq,list(l[:3])+[l[3]]*10+list(l[2::-1]),list(cs[:3])+[cs[3]]*10+list(cs[2::-1]),list(gs[:4])+[gs[4]]*9+list(gs[3::-1]),x[13])[:,1,0]
 def phase(x,frequency=f,actual=obs):
  z=model(x,frequency)*np.exp(2j*np.pi*frequency*x[13]*1e-12)
  objective=lambda d:np.mean(abs(z*np.exp(-2j*np.pi*frequency*d*1e-12)-actual)**2)
  grid=np.linspace(-250,250,501);d=grid[np.argmin([objective(q) for q in grid])];q=minimize_scalar(objective,bounds=(max(-250,d-1),min(250,d+1)),method='bounded');return float(q.x)
 def transfer(x):
  c,l=unpack(x);ls=list(l)+[anchor['L_nH'][4]];cs=list(c[:4])+[anchor['C_fF'][4]];gs=c[4:]
  z=transmission(sf,ls,cs,gs,0)
  objective=lambda d:np.mean(abs(z*np.exp(-2j*np.pi*sf*d*1e-12)-so)**2)
  grid=np.linspace(-250,250,501);d=grid[np.argmin([objective(q) for q in grid])];q=minimize_scalar(objective,bounds=(max(-250,d-1),min(250,d+1)),method='bounded');p=z*np.exp(-2j*np.pi*sf*q.x*1e-12)
  return p,float(q.x)
 report=dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),full_points=len(f),full_source_log=str(FULL/'em.log'),staged_project=chosen['project'],anchor_fit=chosen['fit_path'],policy='Full S21 only is used to fit component values. Small-circuit data is used only after each fit to assess transfer; only its phase delay is adjusted. Symmetric tapers and identical regular cells share parameters. Output test fixture C0 and L are absent from full model and remain at saved values. Calibrations, geometry, selections, and live solver unchanged.',parameter_order=['T3 C0','T2 C0','T1 C0','regular C0','T3 input Cg','T3-T2 Cg','T2-T1 Cg','T1-regular Cg','regular-regular Cg','T3 L','T2 L','T1 L','regular L','full delay'],runs=[])
 def summarize(x,cap_pct,lpct,label):
  c,l=unpack(x);p,sd=transfer(x);pred=model(x)
  return dict(label=label,cap_bound_percent=cap_pct,L_bound_percent=lpct,C_fF=c[:4].tolist(),G_fF=c[4:].tolist(),L_nH=l.tolist(),L_change_percent=((l/bl-1)*100).tolist(),full_delay_ps=float(x[13]),full_complex_RMS=rms(pred-obs),full_magnitude_RMS=rms(abs(pred)-abs(obs)),staged_prediction_complex_RMS=rms(p-so),staged_prediction_magnitude_RMS=rms(abs(p)-abs(so)),staged_delay_ps=sd,output_fixture_C0_fF=anchor['C_fF'][4],output_fixture_L_nH=anchor['L_nH'][4],x=x.tolist())
 bestx=np.zeros(14);bestx[13]=phase(bestx);report['baseline']=summarize(bestx,5,0,'Current extracted parameters')
 previous=[bestx];stored=[]
 # All fixed-L fits precede any fitted-L tests.
 for cap_pct,lpct in [(10,0),(30,0),(50,0),(50,1),(50,3),(50,10),(50,20)]:
  free=np.r_[np.arange(9),np.arange(9,13) if lpct else np.array([],int),13].astype(int)
  lower=np.r_[np.full(9,np.log(1-cap_pct/100)),np.full(4,np.log(1-lpct/100) if lpct else 0),-250]
  upper=np.r_[np.full(9,np.log(1+cap_pct/100)),np.full(4,np.log(1+lpct/100) if lpct else 0),250]
  lo=lower[free];hi=upper[free];scale=np.ones(len(free));scale[-1]=30
  def expand(y):
   x=np.zeros(14);x[free]=y;return x
  def residual(y):return ri(model(expand(y))-obs)
  starts=[]
  for x in previous[:6]:starts.append(np.clip(x[free],lo+1e-9,hi-1e-9))
  for j in range(28):
   parent=previous[j%len(previous)].copy();parent[:9]+=rng.normal(0,.012 if j<14 else .07,9)
   if lpct:parent[9:13]+=rng.normal(0,lpct/350,4)
   parent[free]=np.clip(parent[free],lo+1e-9,hi-1e-9);parent[13]=phase(parent);starts.append(parent[free])
  results=[]
  for j,seed in enumerate(starts):
   save('status.json',dict(state='full_only_multistart',cap_bounds_percent=cap_pct,L_bounds_percent=lpct,start=j+1,total=len(starts),pid=os.getpid()))
   # Magnitude warm-up helps place narrow modes before complex optimization.
   if j%3==0:
    mag=least_squares(lambda y:abs(model(expand(y)))-abs(obs),seed,bounds=(lo,hi),x_scale=scale,max_nfev=200)
    xx=expand(mag.x);xx[13]=phase(xx);seed=xx[free]
   q=least_squares(residual,seed,bounds=(lo,hi),x_scale=scale,max_nfev=650,ftol=1e-10,xtol=1e-10,gtol=1e-10);results.append(q)
  results.sort(key=lambda q:q.fun@q.fun);previous=[expand(q.x) for q in results[:6]];bestx=previous[0]
  rec=summarize(bestx,cap_pct,lpct,f'Full only: caps {cap_pct}%, L {lpct}%');rec['active_bounds_free']=results[0].active_mask.tolist();rec['free_indices']=free.tolist();rec['starts']=[dict(full_RMS=float(np.sqrt(np.mean(abs(model(expand(q.x))-obs)**2))),x=expand(q.x).tolist(),success=bool(q.success)) for q in results]
  sv=np.linalg.svd(results[0].jac*np.r_[np.full(len(free)-1,.01),1],compute_uv=False);rec['local_singular_values']=sv.tolist();rec['local_condition_number']=float(sv[0]/sv[-1]);report['runs'].append(rec);stored.append(bestx.copy());save('report.json',report);print('RESULT',cap_pct,lpct,rec['full_complex_RMS'],rec['staged_prediction_complex_RMS'],flush=True)
 # Choose by full data alone. Staged prediction must not influence selection.
 idx=int(np.argmin([r['full_complex_RMS'] for r in report['runs']]));x=stored[idx];report['best_by_full_only']=report['runs'][idx];p,sd=transfer(x)
 # Held-out full points: keep every fourth frequency out, optimize starting at full solution (optimistic diagnostic).
 mask=np.arange(len(f))%4!=1;idxs=np.where(mask)[0];rec=report['best_by_full_only'];cp=rec['cap_bound_percent'];lp=rec['L_bound_percent'];free=np.r_[np.arange(9),np.arange(9,13) if lp else np.array([],int),13].astype(int)
 lower=np.r_[np.full(9,np.log(1-cp/100)),np.full(4,np.log(1-lp/100) if lp else 0),-250][free];upper=np.r_[np.full(9,np.log(1+cp/100)),np.full(4,np.log(1+lp/100) if lp else 0),250][free]
 def expand2(y):
  z=np.zeros(14);z[free]=y;return z
 q=least_squares(lambda y:ri(model(expand2(y))[mask]-obs[mask]),np.clip(x[free],lower+1e-9,upper-1e-9),bounds=(lower,upper),x_scale=np.r_[np.ones(len(free)-1),30],max_nfev=800)
 h=model(expand2(q.x));report['exploratory_holdout']=dict(training_points=int(mask.sum()),heldout_points=int((~mask).sum()),training_RMS=rms((h-obs)[mask]),heldout_RMS=rms((h-obs)[~mask]),caveat='Initialization and model selection used all full points; this is optimistic, not independent validation.')
 dense=np.linspace(f.min(),f.max(),6001);fig,axes=plt.subplots(2,1,figsize=(14,10))
 axes[0].plot(f/1e9,abs(obs),'o',color='#235d95',label=f'Completed full Sonnet points ({len(f)})')
 fixedidx=int(np.argmin([r['full_complex_RMS'] if r['L_bound_percent']==0 else np.inf for r in report['runs']]))
 axes[0].plot(dense/1e9,abs(model(stored[fixedidx],dense)),'--',label=f"Fixed L fit: RMS {report['runs'][fixedidx]['full_complex_RMS']:.4f}")
 axes[0].plot(dense/1e9,abs(model(x,dense)),label=f"Best full-only fit: RMS {rec['full_complex_RMS']:.4f}")
 axes[1].plot(sf/1e9,abs(so),color='black',label='Small-circuit Sonnet (not fitted)')
 axes[1].plot(sf/1e9,abs(p),'--',label=f"Prediction from full fit: RMS {rec['staged_prediction_complex_RMS']:.4f}")
 for ax in axes:
  ax.set(xlabel='Frequency (GHz)',ylabel='|S21|',ylim=(-.02,1.08));ax.grid(alpha=.2);ax.legend(fontsize=11)
 axes[0].set_title('Full-filter-only fit: sparse completed points; dense lines are circuit predictions')
 axes[1].set_title('Reverse check: component values transferred without refitting')
 fig.tight_layout();fig.savefig(OUT/'full_first_fit_and_prediction.png',dpi=180);plt.close(fig)
 np.savez(OUT/'best_curves.npz',full_f=f,full_actual=obs,dense_f=dense,full_model_dense=model(x,dense),small_f=sf,small_actual=so,small_prediction=p)
 save('report.json',report);save('status.json',dict(state='complete',no_EM_actions=True,report=str(OUT/'report.json')));print('COMPLETE',OUT,flush=True)
if __name__=='__main__':main()
