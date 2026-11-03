"""Read-only progress snapshot for the active cleanup/full-filter controller."""
from pathlib import Path
import json,re
root=Path('output_data/taper_exploration/T2_cleanup_2026-11-01')
s=json.loads((root/'workflow_status.json').read_text(encoding='utf-8-sig'))
r={k:s.get(k) for k in ['phase','case','calibration_stage','full_project','full_completed','reason'] if s.get(k) is not None}
if s.get('full_output') and s['phase']=='full_filter_EM_running':logs=[Path(s['full_output'])/'em.log']
else:logs=list(root.glob('*/Lcal_explicit_*/*.log'))+list(root.glob('*/response/em.log'))
logs=[p for p in logs if p.exists()]
if logs:
 p=max(logs,key=lambda x:x.stat().st_mtime);t=p.read_text(encoding='utf-8-sig',errors='replace');f=re.findall(r'Frequency:\s+([\d.]+) GHZ',t)
 r.update(log=str(p),direct_frequency_headers=len(f),latest_frequency_GHz=float(f[-1]) if f else None,last_log_lines=t.strip().splitlines()[-3:])
 q=Path(str(p)+'.resources.json')
 if q.exists():
  z=json.loads(q.read_text(encoding='utf-8-sig'));r.update(solver_state=z.get('state'),solver_pid=z.get('child_pid'),thread_limit=z.get('thread_limit'))
  a=z.get('samples',[])
  if a:r.update(elapsed_minutes=round(a[-1]['elapsed_seconds']/60,1),active_cores=round((a[-1].get('child_core_equivalents') or 0),2))
print(json.dumps(r))
