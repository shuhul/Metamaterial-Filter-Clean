"""Read-only solver monitoring plus persistent full-filter task state; never starts EM."""
from datetime import datetime,timezone
from pathlib import Path
import json,re
import numpy as np
from geometry.taper.full_filter import OUT

def snapshot():
    log=(OUT/'em.log').read_text(encoding='utf-8-sig')
    resources=json.loads((OUT/'em.log.resources.json').read_text(encoding='utf-8-sig'))
    status=json.loads((OUT/'run_status.json').read_text(encoding='utf-8-sig'))
    sample=resources['samples'][-1]
    rows=[]
    pattern=r'^\s+(\d+\.\d+)\s+([\d.]+)\s+([-\d.]+)\s+([\d.]+)\s+([-\d.]+)\s+([\d.]+)\s+([-\d.]+)\s+([\d.]+)\s+([-\d.]+)\s*$'
    for match in re.finditer(pattern,log,re.M):
        x=list(map(float,match.groups()))
        rows.append(dict(f_GHz=x[0],S11_mag=x[1],S21_magnitude=x[3],S21_phase_deg=x[4],S21_dB=float(20*np.log10(max(x[3],1e-15)))))
    frequencies=re.findall(r'Frequency:\s+([\d.]+) GHZ',log)
    current=float(frequencies[-1]) if frequencies else None
    plan_path=Path('geometry/taper/full_filter_work_plan.json')
    plan=json.loads(plan_path.read_text(encoding='utf-8-sig'))
    plan.update(completed_direct_samples=rows,full_state=status['status'],full_current_frequency_GHz=current,
                full_monitor_last_seen_utc=datetime.now(timezone.utc).isoformat(),
                full_curve_convergence_pending=status['status']!='complete')
    plan_path.write_text(json.dumps(plan,indent=2)+'\n',encoding='utf-8')
    return dict(elapsed_min=round(sample['elapsed_seconds']/60,1),active_cores=round(sample['child_core_equivalents'],2),
                completed_frequencies=len(rows),current_frequency_GHz=current,state=status['status'],
                native_log_length=len(log),tail=log[-1400:])

if __name__=='__main__':print(json.dumps(snapshot()))
