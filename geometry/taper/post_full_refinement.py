"""Wait for complete full-filter analysis, then refine saved staged fits only.

No EM launch, termination or geometry edits occur in this script. Measured
inductances remain immutable; optional fitted L experiments are labeled.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
from datetime import datetime,timezone
import json,time
import numpy as np
from scipy.optimize import least_squares
from geometry.taper.progressive_fitting import fit,transmission
from circuit_common import read_snp
from geometry.specifications import TAPER_3,TAPER_2,TAPER_1,UNIT_CELL

ROOT=Path(__file__).resolve().parents[2]
CAMPAIGN=ROOT/'output_data/taper_exploration/T2_cleanup_2026-11-01'
DEST=ROOT/'output_data/taper_exploration/post_full_refinement_2026-11-02'

def read(p):return json.loads(Path(p).read_text(encoding='utf-8-sig'))
def save(p,x):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(x,indent=2)+'\n',encoding='utf-8')

def main():
 DEST.mkdir(parents=True,exist_ok=True)
 save(DEST/'status.json',dict(state='waiting_for_full_simulation_and_analysis',pid=os.getpid(),started_utc=datetime.now(timezone.utc).isoformat(),no_geometry_or_EM_actions=True))
 while True:
  state=read(CAMPAIGN/'workflow_status.json')
  if state['phase']=='complete' and state.get('full_completed'):break
  if state['phase']=='failed':
   save(DEST/'status.json',dict(state='waiting_on_failed_full_workflow',reason=state.get('reason'),new_fits_started=False));return
  time.sleep(50)
 full=Path(state['full_output'])
 assert read(full/'run_status.json')['status']=='complete'
 assert read(full/'response_integrity.json')['basic_integrity_pass']
 assert (full/'full_filter_analysis.json').exists()
 chosen=read(CAMPAIGN/'selected_cleanup.json');anchor=read(chosen['fit_path'])
 source=Path(chosen['project']).with_suffix('.snp')
 save(DEST/'anchor_fit.json',anchor)
 report=dict(full_simulation_completed_before_fits=True,full_analysis=str(full/'full_filter_analysis.json'),source=str(source),baseline_RMS=anchor['complex_RMS'],capacitance_tests=[],inductance_tests=[],policy='No original geometry/calibration/selection overwritten. All cap bounds share immutable current anchor. Fitted L is a last-resort diagnostic, not a replacement for measured L. Geometry tuning remains a separate subsequent step.')
 save(DEST/'status.json',dict(state='capacitance_only_fits',full_completed=True,pid=os.getpid()))
 free=[(k,i) for i in range(5) for k in ('C','G')]
 values=[anchor['C_fF' if k=='C' else 'G_fF'][i] for k,i in free]
 previous=None
 for pct in [1,2,3,5,10,20]:
  seeds=[values]
  if previous:seeds.append([previous['C_fF' if k=='C' else 'G_fF'][i] for k,i in free])
  fitted=fit(source,DEST/f'caps_{pct}pct',anchor['L_nH'],anchor['C_fF'],anchor['G_fF'],free,capacitance_bounds=[(v*(1-pct/100),v*(1+pct/100)) for v in values],capacitance_seeds=seeds,max_nfev=900,title=f'Post-Full Check: Capacitances Within {pct}% of Current Fit')
  report['capacitance_tests'].append(dict(bounds_percent=pct,RMS=fitted['complex_RMS'],fit_path=str(DEST/f'caps_{pct}pct/fit.json')))
  save(DEST/'report.json',report);previous=fitted
  if fitted['complex_RMS']<.005:break
 if previous['complex_RMS']>=.005:
  # All capacitance-only tests have completed before any inductance is freed.
  net=read_snp(source);obs=net.s[:,1,0];ix=np.unique(np.linspace(0,len(net.f)-1,min(2001,len(net.f))).astype(int))
  cap_seed=np.array([previous['C_fF' if k=='C' else 'G_fF'][i] for k,i in free]);ls=np.array(anchor['L_nH'])
  for pct in [.25,.5,1,2]:
   save(DEST/'status.json',dict(state='last_resort_fitted_L_diagnostic',L_bound_percent=pct,measured_L_unchanged=True))
   lo=np.r_[np.array(values)*.8,ls*(1-pct/100),-100.];hi=np.r_[np.array(values)*1.2,ls*(1+pct/100),100.]
   initial=np.r_[cap_seed,ls,previous['delay_ps']];scale=np.r_[np.maximum(cap_seed,8),ls,20.]
   def residual(x,indices):
    cs=x[:10:2];gs=x[1:10:2];model=transmission(net.f[indices],x[10:15],cs,gs,x[-1]);d=model-obs[indices];return np.r_[d.real,d.imag]
   rng=np.random.default_rng(1202);starts=[initial]
   for _ in range(2):starts.append(np.clip(initial+rng.normal(0,1,16)*(hi-lo)*.02,lo+1e-9,hi-1e-9))
   fits=[least_squares(lambda x:residual(x,ix),seed,bounds=(lo,hi),x_scale=scale,max_nfev=900,ftol=1e-10,xtol=1e-10,gtol=1e-10) for seed in starts]
   best=min(fits,key=lambda q:q.fun@q.fun);best=least_squares(lambda x:residual(x,np.arange(len(net.f))),best.x,bounds=(lo,hi),x_scale=scale,max_nfev=900,ftol=1e-11,xtol=1e-11,gtol=1e-11)
   model=transmission(net.f,best.x[10:15],best.x[:10:2],best.x[1:10:2],best.x[-1]);rms=float(np.sqrt(np.mean(abs(model-obs)**2)))
   record=dict(L_bound_percent=pct,capacitance_bound_percent=20,RMS=rms,L_nH=best.x[10:15].tolist(),C_fF=best.x[:10:2].tolist(),G_fF=best.x[1:10:2].tolist(),delay_ps=float(best.x[-1]),optimizer_success=bool(best.success),active_bounds=best.active_mask.tolist(),scope='Fitted-L diagnostic only; calibration unchanged')
   save(DEST/f'fitted_L_{pct}pct.json',record);report['inductance_tests'].append(record);save(DEST/'report.json',report)
   if rms<.005:break
 targets=[TAPER_3,TAPER_2,TAPER_1,UNIT_CELL,UNIT_CELL];cs=[t.C0_fF for t in targets];gs=[t.Cg_input_fF for t in targets];tv=[(cs if k=='C' else gs)[i] for k,i in free]
 target_fit=fit(source,DEST/'all_caps_target_1pct',anchor['L_nH'],cs,gs,free,capacitance_bounds=[(.99*v,1.01*v) for v in tv],capacitance_seeds=[values],max_nfev=900,title='Current Geometry: All Capacitances Within 1% of Targets')
 report['target_1pct_fit_RMS']=target_fit['complex_RMS'];report['target_1pct_fit_path']=str(DEST/'all_caps_target_1pct/fit.json')
 save(DEST/'report.json',report)
 save(DEST/'status.json',dict(state='saved_data_refinement_complete_geometry_iteration_next',completed_utc=datetime.now(timezone.utc).isoformat(),report=str(DEST/'report.json'),no_geometry_or_calibration_overwritten=True))
 print('POST_FULL_REFINEMENT_COMPLETE',DEST/'report.json',flush=True)

if __name__=='__main__':
 try:main()
 except BaseException as e:save(DEST/'status.json',dict(state='failed',reason=str(e)));raise
