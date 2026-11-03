"""Offline differential circuit hints, never geometry acceptance or extraction.

Multiply the measured curve by a relative circuit response to retain its
existing residual. The result only suggests which physical dimensions to test.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import numpy as np
from scipy.optimize import least_squares
from geometry.taper import taper2_plain as w


def hint(case, L2=None, *, mode='ratio'):
    if mode not in ('ratio','inverse'):raise ValueError('Unknown differential model')
    folder=w.ROOT/'output_data/taper_exploration'/w.CAMPAIGN/case/'response'
    a=np.load(folder/'comparison.npz');fit=json.loads((folder/'extraction/fit.json').read_text())
    idx=np.arange(0,len(a['f']),5);f=a['f'][idx];obs=a['observed'][idx]
    old=w.common.transmission(f,fit['L_nH'],fit['C_fF'],fit['G_fF'],fit['delay_ps'])
    specs=[w.common.TAPER_2,w.common.TAPER_1,w.common.UNIT_CELL,w.common.UNIT_CELL]
    target=w.common.transmission(f,[s.L0_nH for s in specs],[s.C0_fF for s in specs],
                                  [s.Cg_input_fF for s in specs],0)
    ls=list(fit['L_nH'])
    if L2 is not None:ls[0]=L2
    def curve(q):
        c=list(fit['C_fF']);g=list(fit['G_fF'])
        c[0],c[1],g[0],g[1],g[2],tau=q
        new=w.common.transmission(f,ls,c,g,fit['delay_ps'])
        predicted=obs*new/old if mode=='ratio' else 1/(1/new+1/obs-1/old)
        return predicted*np.exp(-2j*np.pi*f*tau*1e-12)
    def residual(q):
        d=curve(q)-target
        return np.r_[d.real,d.imag]
    initial=[*fit['C_fF'][:2],*fit['G_fF'][:3],0]
    starts=[initial]
    for c2 in (80,115,150,200):
        for g0 in (100,180,280):
            starts.append([c2,148,g0,15.3,8.2865,0])
    results=[least_squares(residual,seed,bounds=([10,135,40,5,6,-100],[300,170,500,40,12,100]),
                         x_scale='jac',max_nfev=400) for seed in starts]
    result=min(results,key=lambda r:r.fun@r.fun)
    report=dict(scope='Offline differential hint only; not extracted or accepted values',source=case,
        differential_model=mode,assumed_L2_nH=ls[0],parameters=['C2_fF','C1_fF','input_Cg_fF','T2_T1_Cg_fF','T1_unit_Cg_fF','delay_ps'],
        source_values=initial,predicted_values=result.x.tolist(),predicted_RMS=float(np.sqrt(np.mean(abs(curve(result.x)-target)**2))),
        starting_points=len(starts),local_solutions=[dict(values=r.x.tolist(),RMS=float(np.sqrt(np.mean(abs(curve(r.x)-target)**2)))) for r in results])
    filename='electrical_dimension_hint.json' if mode=='ratio' else 'electrical_inverse_hint.json'
    (folder/filename).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k!='local_solutions'},indent=2));return report


if __name__=='__main__':
    hint('idc_rail8_dl256_pad148_w152_i32n70_j42_r78tg2')
