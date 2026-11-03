"""Preliminary geometry guide from completed direct points of a running sweep.

Never accepts a design. The complete dense sweep is required for acceptance.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import re,json
import numpy as np
from scipy.optimize import least_squares
from geometry.taper.three_tapers import ROOT,CAMPAIGN
from geometry.taper.staged_fitting import circuit


def analyze(name,*,campaign=CAMPAIGN,input_end_cell=True):
    out=ROOT/'output_data/taper_exploration'/campaign/name
    lines=(out/'response/em.log').read_text().splitlines()
    rows=[]
    for line in lines:
        if not re.match(r'^\s+\d+\.\d{8}\s',line):continue
        v=list(map(float,line.split()))
        if len(v)==9:rows.append(v)
    if len(rows)<3:return
    data=np.array(rows);f=data[:,0]*1e9
    complex_values=data[:,1::2]*np.exp(1j*np.deg2rad(data[:,2::2]))
    observed=complex_values[:,[0,1,3]]
    l=json.loads((out/'Lcal/results.json').read_text())['L0_nH']
    def pred(x):
        ls=[l,l,4.922723945,4.922723945] if input_end_cell else [l,4.922723945,4.922723945]
        cs=[x[0],x[0],188.878116386,188.878116386] if input_end_cell else [x[0],188.878116386,188.878116386]
        gs=[None,x[1],x[2],8.241620220] if input_end_cell else [x[1],x[2],8.241620220]
        s=circuit(f,ls,cs,gs,left_ps=x[3],right_ps=x[4])
        return np.stack([s[:,0,0],s[:,1,0],s[:,1,1]],axis=1)
    def residual(x):
        delta=pred(x)-observed
        return np.r_[delta.real.ravel(),delta.imag.ravel()]
    rng=np.random.default_rng(91)
    starts=[[160,14,8,5,5]]+[[rng.uniform(140,220),rng.uniform(8,22),rng.uniform(5,12),5,5] for _ in range(19)]
    fits=[least_squares(residual,x,bounds=([80,1,1,-100,-100],[300,50,20,100,100]),
                        x_scale=[150,12,8,20,20],max_nfev=300) for x in starts]
    best=min(fits,key=lambda q:q.fun@q.fun)
    result={'scope':'PRELIMINARY sparse direct-point guide; not design acceptance or final extraction',
            'frequencies_GHz':data[:,0].tolist(),'parameters_C0_Gin_Gunit_leftps_rightps':best.x.tolist(),
            'joint_complex_RMS':float(np.sqrt(np.mean(abs(pred(best.x)-observed)**2))),
            'fixed_L_nH':l,'input_end_LC_present':input_end_cell}
    (out/'preliminary_direct_point_guide.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result,flush=True)
    return result


if __name__=='__main__':analyze('w224_pad8_d20')
