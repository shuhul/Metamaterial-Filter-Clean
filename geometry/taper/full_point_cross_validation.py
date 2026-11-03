"""Read-only saved-data cross-validation. Never launches or controls Sonnet."""
import os
for k in ['OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS']:os.environ[k]='1'
from pathlib import Path
import json,hashlib
from datetime import datetime,timezone
import numpy as np
from scipy.optimize import least_squares,minimize_scalar
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from circuit_common import read_snp
from geometry.taper.progressive_fitting import transmission
from geometry.taper.full_filter import ladder
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'output_data/taper_exploration/full_point_cross_validation_2026-11-02'
FULL=ROOT/'output_data/taper_exploration/full_filter/full_filter_10_regular_3_tapers_each_side_T2_cleanup_2026-11-01'
def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def save(name,obj):
 (OUT/name).write_text(json.dumps(obj,indent=2)+'\n',encoding='utf-8')
def packed(z):return np.r_[z.real,z.imag]
def rms(z):return float(np.sqrt(np.mean(abs(z)**2)))
def main():
 OUT.mkdir(parents=True,exist_ok=True)
 chosen=read(ROOT/'output_data/taper_exploration/T2_cleanup_2026-11-01/selected_cleanup.json');anchor=read(chosen['fit_path'])
 source=Path(chosen['project']).with_suffix('.snp');net=read_snp(source);f=net.f;obs=net.s[:,1,0]
 log=(FULL/'em.log').read_text(encoding='utf-8-sig');(OUT/'em_snapshot.log').write_text(log,encoding='utf-8')
 rows=[]
 for line in log.splitlines():
  try:v=list(map(float,line.split()))
  except ValueError:continue
  if len(v)==9 and 4.65<=v[0]<=5.4:rows.append(v)
 a=np.array(sorted(rows,key=lambda v:v[0]));ff=a[:,0]*1e9;fo=a[:,3]*np.exp(1j*np.deg2rad(a[:,4]))
 np.savez(OUT/'data_snapshot.npz',staged_f=f,staged_s21=obs,full_f=ff,full_s21=fo)
 ls=np.array(anchor['L_nH']);base=np.array([v for pair in zip(anchor['C_fF'],anchor['G_fF']) for v in pair])
 ix=np.unique(np.linspace(0,len(f)-1,1501).astype(int));train=np.arange(len(ff))%4!=1;test=~train
 def stage(x,indices=ix):return transmission(f[indices],ls,x[:10:2],x[1:10:2],x[10])
 def full(x,freq=ff,delay=0):
  cs=x[:10:2];gs=x[1:10:2]
  return ladder(freq,list(ls[:3])+[ls[3]]*10+list(ls[2::-1]),list(cs[:3])+[cs[3]]*10+list(cs[2::-1]),list(gs[:4])+[gs[4]]*9+list(gs[3::-1]),delay)[:,1,0]
 def delay_check(x,mask=None):
  mask=np.ones(len(ff),bool) if mask is None else mask
  z=full(x);grid=np.linspace(-200,200,401)
  fun=lambda d:np.mean(abs(z[mask]*np.exp(-2j*np.pi*ff[mask]*d*1e-12)-fo[mask])**2)
  d0=grid[np.argmin([fun(d) for d in grid])];opt=minimize_scalar(fun,bounds=(max(-200,d0-1),min(200,d0+1)),method='bounded')
  pred=z*np.exp(-2j*np.pi*ff*opt.x*1e-12)
  return float(opt.x),pred
 report=dict(timestamp_utc=datetime.now(timezone.utc).isoformat(),full_direct_points=len(ff),staged_source=str(source),staged_sha256=hashlib.sha256(source.read_bytes()).hexdigest(),anchor_path=chosen['fit_path'],fixed_L_nH=ls.tolist(),policy='Snapshot only. No EM controls or selected design changes. Staged-only fits predict full points without component refitting; full phase delay is a nuisance fitted separately. Joint fit uses 75% of sorted direct frequencies, every fourth held out. Adaptive samples are nonuniform and held-out points are exploratory, not an independent EM simulation.',candidates=[])
 rng=np.random.default_rng(1102);candidates=[]
 def record(x,label,pct):
  d,pred=delay_check(x);dt,pt=delay_check(x,train)
  return dict(label=label,bounds_percent=pct,C_fF=x[:10:2].tolist(),G_fF=x[1:10:2].tolist(),staged_delay_ps=float(x[10]),staged_RMS=rms(stage(x,np.arange(len(f)))-obs),full_complex_RMS=rms(pred-fo),full_magnitude_RMS=rms(abs(pred)-abs(fo)),full_delay_ps=d,heldout_complex_RMS=rms((pt-fo)[test]),train_delay_ps=dt)
 baseline=np.r_[base,anchor['delay_ps']];report['baseline']=record(baseline,'Current bounded fit',5)
 for pct in [10,20,50]:
  save('status.json',dict(state='staged_multistart',bounds_percent=pct,pid=os.getpid()))
  lo=np.r_[base*(1-pct/100),-100];hi=np.r_[base*(1+pct/100),100];scale=np.r_[base,20]
  starts=[baseline]+[np.clip(baseline+rng.normal(size=11)*(hi-lo)*.20,lo+1e-7,hi-1e-7) for _ in range(7)]
  solutions=[]
  for seed in starts:
   q=least_squares(lambda x:packed(stage(x)-obs[ix]),seed,bounds=(lo,hi),x_scale=scale,max_nfev=400,ftol=1e-9,xtol=1e-9,gtol=1e-9);solutions.append(q)
  ordered=sorted(solutions,key=lambda q:q.fun@q.fun)
  for j,q in enumerate(ordered[:3]):
   if j==0:q=least_squares(lambda x:packed(stage(x,np.arange(len(f)))-obs),q.x,bounds=(lo,hi),x_scale=scale,max_nfev=300,ftol=1e-10,xtol=1e-10,gtol=1e-10)
   x=q.x;rec=record(x,f'Staged only {pct}% #{j+1}',pct);rec['active_bounds']=q.active_mask.tolist();rec['success']=bool(q.success);report['candidates'].append(rec);candidates.append((rec,x))
  save('report.json',report);print('STAGED',pct,report['candidates'][-3],flush=True)
 # Joint fit: components now deliberately trained on both circuits, held-out full points excluded.
 save('status.json',dict(state='joint_fit_and_identifiability',pid=os.getpid()))
 lo=np.r_[base*.5,-100,-200];hi=np.r_[base*1.5,100,200];scale=np.r_[base,20,40]
 def jr(x,include_full=True):
  z=packed(stage(x)-obs[ix])/np.sqrt(len(ix))
  if not include_full:return z
  zz=packed(full(x,ff[train],x[11])-fo[train])/np.sqrt(train.sum())
  return np.r_[z,zz]
 seeds=sorted(candidates,key=lambda v:v[0]['full_complex_RMS'])[:4]+[({},baseline)]
 for j in range(24):
  parent=baseline if j%2==0 else candidates[j%len(candidates)][1]
  trial=np.clip(parent+rng.normal(size=11)*np.r_[base*(.01 if j<12 else .05),4],lo[:11]+1e-7,hi[:11]-1e-7)
  seeds.append(({},trial))
 joint=[]
 for rec,x in seeds:
  d,_=delay_check(x,train)
  q=least_squares(jr,np.r_[x,d],bounds=(lo,hi),x_scale=scale,max_nfev=800,ftol=1e-11,xtol=1e-11,gtol=1e-11);joint.append(q)
 best=min(joint,key=lambda q:q.fun@q.fun);x=best.x;pred=full(x,ff,x[11])
 # Relative sensitivity columns: a unit step means 1% of anchor capacitance (1 ps delay).
 steps=np.r_[base*.01,1,1]
 J=best.jac*steps;Js=J[:2*len(ix),:11]
 def spectrum(j):
  u,s,v=np.linalg.svd(j,full_matrices=False)
  return dict(singular_values=s.tolist(),condition_number=float(s[0]/s[-1]),weakest_direction=v[-1].tolist(),note='Local sensitivity only; no measurement-noise model or uniqueness guarantee.')
 report['joint']=dict(C_fF=x[:10:2].tolist(),G_fF=x[1:10:2].tolist(),staged_delay_ps=float(x[10]),full_delay_ps=float(x[11]),staged_RMS=rms(stage(x,np.arange(len(f)))-obs),full_train_complex_RMS=rms((pred-fo)[train]),full_heldout_complex_RMS=rms((pred-fo)[test]),full_all_complex_RMS=rms(pred-fo),full_magnitude_RMS=rms(abs(pred)-abs(fo)),active_bounds=best.active_mask.tolist(),optimizer_success=bool(best.success),staged_local_sensitivity=spectrum(Js),joint_local_sensitivity=spectrum(J),capacitance_parameter_order=['T3 C0','T3 input Cg','T2 C0','T2 input Cg','T1 C0','T1 input Cg','regular C0','T1-regular Cg','test output C0','regular-regular Cg'])
 # Profile each cap +/-5% around joint solution, reoptimize remaining parameters.
 profiles=[]
 for k in range(10):
  for direction in [-1,1]:
   fixed=float(np.clip(x[k]+direction*.05*base[k],lo[k]+1e-7,hi[k]-1e-7));mask=np.arange(12)!=k
   def expand(y):
    z=x.copy();z[k]=fixed;z[mask]=y;return z
   q=least_squares(lambda y:jr(expand(y)),x[mask],bounds=(lo[mask],hi[mask]),x_scale=scale[mask],max_nfev=300,ftol=1e-9,xtol=1e-9,gtol=1e-9);z=expand(q.x);p=full(z,ff,z[11])
   profiles.append(dict(cap_index=k,fixed_fF=fixed,delta_anchor_percent=float((fixed-x[k])/base[k]*100),joint_objective=float(q.fun@q.fun),objective_ratio=float((q.fun@q.fun)/(best.fun@best.fun)),staged_RMS=rms(stage(z,np.arange(len(f)))-obs),heldout_RMS=rms((p-fo)[test])))
 report['joint_cap_profiles']=profiles;save('report.json',report)
 dense=np.linspace(ff.min(),ff.max(),3001);fig,axes=plt.subplots(2,1,figsize=(13,10))
 axes[0].plot(f/1e9,abs(obs),color='black',label='Staged Sonnet')
 axes[0].plot(f/1e9,abs(stage(x,np.arange(len(f)))),'--',label=f"Joint fit, RMS {report['joint']['staged_RMS']:.4f}")
 axes[1].plot(ff[train]/1e9,abs(fo[train]),'o',label='Full direct points: training')
 axes[1].plot(ff[test]/1e9,abs(fo[test]),'s',label='Full direct points: held out')
 axes[1].plot(dense/1e9,abs(full(baseline,dense)),':',label='Original bounded prediction')
 ranked=sorted(candidates,key=lambda v:v[0]['staged_RMS']);rec,bx=ranked[0]
 axes[1].plot(dense/1e9,abs(full(bx,dense)),'--',label='Best staged-only prediction')
 axes[1].plot(dense/1e9,abs(full(x,dense)),label='Joint fit prediction')
 for ax in axes:
  ax.set(xlabel='Frequency (GHz)',ylabel='|S21|',ylim=(-.02,1.08));ax.grid(alpha=.25);ax.legend(fontsize=10)
 axes[0].set_title('Saved-data fit: calibrated inductances fixed')
 axes[1].set_title(f'Full-filter snapshot: {len(ff)} completed direct points; sweep still running')
 fig.tight_layout();fig.savefig(OUT/'cross_validation.png',dpi=170);plt.close(fig)
 save('status.json',dict(state='complete',report=str(OUT/'report.json'),no_EM_actions=True));print('COMPLETE',OUT,report['joint'],flush=True)
if __name__=='__main__':main()
