"""Reuse archived EM data to rank unchanged-T1 cases against the specification."""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks
from circuit_common import read_snp
from geometry.taper.optimize_response import ROOT,OUT,theory
from geometry.taper.second_layout import from_dict
from geometry.taper.locked_taper2 import require_locked

rows=[]
for p in (ROOT/'output_data/taper_exploration').glob('*/response.snp'):
    meta=ROOT/'geometry/taper/output'/f'{p.parent.name}.json'
    try:
        m=json.loads(meta.read_text());require_locked(from_dict(m['second']))
        n=read_snp(p);mask=(n.f>=4.5e9)&(n.f<=5.5e9)
        if n.f[0]>4.5e9 or n.f[-1]<5.5e9:continue
        f=n.f[mask];actual=n.s[mask,1,0];target=theory(f)[:,1,0]
        phase=lambda tau:float(np.mean(abs(actual-target*np.exp(-2j*np.pi*f*tau*1e-12))**2))
        delay=minimize_scalar(phase,bounds=(-40,40),method='bounded')
        ix,_=find_peaks(abs(actual),prominence=.02)
        rows.append({'case':p.parent.name,'complex_RMSE':float(np.sqrt(phase(0))),
            'delay_aligned_complex_RMSE':float(np.sqrt(delay.fun)),'extra_reference_delay_ps':float(delay.x),
            'magnitude_RMSE':float(np.sqrt(np.mean((abs(actual)-abs(target))**2))),
            'peaks':list(zip((f[ix]/1e9).tolist(),abs(actual[ix]).tolist()))})
    except (KeyError,ValueError,FileNotFoundError):continue
rows.sort(key=lambda q:q['magnitude_RMSE'])
(OUT/'archive_comparison.json').write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(rows,indent=2))
