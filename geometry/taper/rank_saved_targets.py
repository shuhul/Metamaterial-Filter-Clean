"""Cheap historical guidance only; archived layouts may violate current rules."""
import json
import numpy as np
from geometry.taper.target_iteration import ROOT, OUT, START, STOP, target_response
from circuit_common import read_snp


def run():
    rows=[]
    for path in (ROOT/'sonnet_files/taper_exploration').rglob('generated_ncell_3.snp'):
        try:
            net=read_snp(path)
            if net.s.shape[1:]!=(2,2) or len(net.f)<101:continue
            mask=(net.f>=START*1e9)&(net.f<=STOP*1e9)
            f=net.f[mask];s=net.s[mask,1,0]
            if len(f)<101 or f[0]>START*1e9+1e6 or f[-1]<STOP*1e9-1e6:continue
            rms=float(np.sqrt(np.mean((abs(s)-abs(target_response(f)))**2)))
            rows.append(dict(path=str(path.relative_to(ROOT)),magnitude_RMS=rms,
                             scope='Historical guidance; topology and current geometry constraints not validated'))
        except (ValueError,IndexError):continue
    rows.sort(key=lambda r:r['magnitude_RMS'])
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/'historical_magnitude_ranking.json').write_text(json.dumps(rows,indent=2)+'\n')
    print(json.dumps(rows[:12],indent=2))


if __name__=='__main__':run()
