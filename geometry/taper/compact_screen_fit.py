"""Analyze two-point geometry screens; never claim design acceptance.

The wire-removed control estimates the junction Cg first. With that estimate
fixed, the grounded-wire case tests whether a nonnegative bare shunt C can
explain both observations. A dense sweep remains mandatory before acceptance.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import numpy as np
from scipy.optimize import least_squares
from circuit_common import read_snp
from geometry.taper.compact_stages import OUT
from geometry.taper.compact_fit import direct_predictor,GJ
from geometry.taper.staged_fitting import circuit,LT1,CT1,LR,CR,GOUT


def coupling(name):
    p=OUT/name/'direct';n=read_snp(p/'screen.snp');f=n.f;a=n.s[:,1,0]
    def pred(x):
        return circuit(f,[LT1,LR,LR],[CT1,CR,CR],[x[0],GJ,GOUT],left_ps=x[1])[:,1,0]
    def res(x):
        d=pred(x)-a;return np.r_[d.real,d.imag]
    fits=[least_squares(res,[g,0],bounds=([5,-30],[180,30]),x_scale=[57,5],max_nfev=100)
          for g in np.arange(10,170,2)]
    best=min(fits,key=lambda q:q.fun@q.fun)
    result=dict(diagnostic_only=True,accepted=False,two_points_insufficient_for_acceptance=True,
                Cg_fF=float(best.x[0]),delay_ps=float(best.x[1]),
                complex_RMS=float(np.sqrt(np.mean(abs(pred(best.x)-a)**2))),
                frequencies_GHz=(f/1e9).tolist(),magnitude=abs(a).tolist(),phase_deg=np.angle(a,deg=True).tolist())
    (p/'two_point_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')
    print(name,json.dumps(result),flush=True)
    return result


def assembly(name,wire_absent_name=None):
    g=coupling(wire_absent_name or name+'_wire_absent')['Cg_fF']
    p=OUT/name/'direct';n=read_snp(p/'screen.snp');a=n.s[:,1,0]
    l=json.loads((OUT/name/'Lcal/results.json').read_text())['L0_nH']
    predict=direct_predictor(n.f,l)
    def res(x):
        d=predict([x[0],g,x[1]])-a;return np.r_[d.real,d.imag]
    fits=[least_squares(res,[c,0],bounds=([0,-30],[2000,30]),x_scale=[100,5],max_nfev=100)
          for c in (0,30,100,500,1000)]
    best=min(fits,key=lambda q:q.fun@q.fun)
    result=dict(diagnostic_only=True,accepted=False,two_points_insufficient_for_acceptance=True,
        fixed_Cg_fF=g,L_T2_nH=l,C_bare_fF=float(best.x[0]),delay_ps=float(best.x[1]),
        complex_RMS=float(np.sqrt(np.mean(abs(predict([best.x[0],g,best.x[1]])-a)**2))),
        active_bounds=best.active_mask.tolist(),observed_magnitude=abs(a).tolist(),
        predicted_magnitude=abs(predict([best.x[0],g,best.x[1]])).tolist())
    (p/'paired_two_point_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')
    print(name,json.dumps(result),flush=True)
    return result
