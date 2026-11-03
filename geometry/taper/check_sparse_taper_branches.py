"""Check smaller subsets using starts projected from the passing large move.

An extra multistart check avoids calling a local near-baseline minimum global.
No regular/output C or L can move; no EM calculation.
"""
import os
for key in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS'):os.environ[key]='1'
import json,itertools
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from geometry.taper.minimal_cap_changes import OUT,ALLOWED
from geometry.taper.progressive_fitting import transmission
from geometry.specifications import TAPER_1,TAPER_2,UNIT_CELL
from circuit_common import read_snp


def main():
    saved=json.loads((OUT/'results.json').read_text());fit=json.loads(Path(saved['source_fit']).read_text())
    f=read_snp(Path(fit['provenance']['source'])).f;base=np.array(saved['baseline_C_fF']+saved['baseline_G_fF'])
    ls=saved['fixed_L_nH'];delay=fit['delay_ps'];ts=[TAPER_2,TAPER_1,UNIT_CELL,UNIT_CELL]
    target=transmission(f,[t.L0_nH for t in ts],[t.C0_fF for t in ts],[t.Cg_input_fF for t in ts],0)
    ix=np.unique(np.linspace(0,len(f)-1,1501).astype(int));full=np.arange(len(f))
    cache={tuple(r['indices']):r for group in saved['search_by_count'].values() for r in group}
    anchors=sorted(saved['search_by_count']['4'],key=lambda r:r['complex_RMS'])[:3]
    rows=[]
    for count in (2,3):
        for group in itertools.combinations(ALLOWED,count):
            free=list(group);old=cache[group]
            lower=np.r_[np.log(np.array(saved['cap_lower_bounds_fF'])[free]/base[free]),-5.]
            upper=np.r_[np.log(np.array(saved['cap_upper_bounds_fF'])[free]/base[free]),5.]
            def expand(z):
                cp=base.copy();cp[free]=base[free]*np.exp(z[:-1]);return cp
            def residual(z,indices=ix):
                cp=expand(z);pred=transmission(f[indices],ls,cp[:4],cp[4:],delay)
                dd=pred-target[indices]*np.exp(-2j*np.pi*f[indices]*z[-1]*20e-12)
                return np.r_[dd.real,dd.imag]
            starts=[]
            for row in [old]+anchors:
                cp=np.array(row['C_fF']+row['G_fF'])
                starts.append(np.clip(np.r_[np.log(cp[free]/base[free]),row['target_delay_ps']/20],lower+1e-8,upper-1e-8))
            fits=[least_squares(residual,z,bounds=(lower,upper),max_nfev=500,
                   ftol=1e-11,xtol=1e-11,gtol=1e-11) for z in starts]
            q=min(fits,key=lambda q:np.linalg.norm(q.fun))
            q=least_squares(lambda z:residual(z,full),q.x,bounds=(lower,upper),max_nfev=500,
                            ftol=1e-12,xtol=1e-12,gtol=1e-12)
            cp=expand(q.x);rms=float(np.sqrt(np.mean(abs(residual(q.x,full))**2)*2))
            row=dict(indices=free,names=old['names'],old_fF=base[free].tolist(),new_fF=cp[free].tolist(),
                changes_percent=(100*(cp[free]/base[free]-1)).tolist(),C_fF=cp[:4].tolist(),G_fF=cp[4:].tolist(),
                complex_RMS=rms,target_delay_ps=float(q.x[-1]*20),passes=bool(rms<.025),previous_RMS=old['complex_RMS'])
            rows.append(row);print(row['names'],rms,flush=True)
            (OUT/'additional_branch_checks.json').write_text(json.dumps(rows,indent=2)+'\n')
        if any(r['passes'] for r in rows):break
    print('BEST',json.dumps(min(rows,key=lambda r:r['complex_RMS'])),flush=True)


if __name__=='__main__':main()
