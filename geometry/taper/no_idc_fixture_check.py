"""Offline check: distinguish 50-ohm normalization from physical feed impedance.

The Sonnet port log gives ~56.887 ohms and effective epsilon 6.5538. Those
are fixed here; they are not extra fitted circuit parameters. This is a
diagnostic approximation to the feed, not an automatic model replacement.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import numpy as np
from scipy.optimize import least_squares
from geometry.taper.no_idc_verification import OUT
from geometry.taper.staged_fitting import circuit


def run(*,free_impedance=False):
    source=OUT/'physical'
    r=json.loads((source/'fit.json').read_text());data=np.load(source/'response.npz')
    f=data['f'];obs=data['actual'][:,1,0];ix=np.arange(0,len(f),3)
    # Remaining port-to-node lengths after the existing 400 um deembedding.
    delays=np.array([836-400,2562-1985-400])*1e-6*np.sqrt(6.5538)/299792458*1e12
    def predict(x,frequency):
        return circuit(frequency,r['L_nH'],[0,x[0],r['C_fF'][2],r['C_fF'][3]],
            [None,x[1],x[2],r['G_fF'][3]],left_ps=delays[0],right_ps=delays[1],
            line_z=x[4] if free_impedance else 56.8874)[:,1,0]*np.exp(-2j*np.pi*frequency*x[3]*1e-12)
    def residual(x,smooth=False,indices=ix):
        actual=obs[indices];m=predict(x,f[indices])
        d=(1/m-1/actual)*np.maximum(abs(actual),.03)**2 if smooth else m-actual
        return np.r_[d.real,d.imag]
    seed=np.array([r['C_fF'][1],r['G_fF'][1],r['G_fF'][2],0.])
    bounds=([0,.1,.1,-100],[1000,3000,3000,100]);scale=[100,30,8,20]
    if free_impedance:
        seed=np.r_[seed,100.];bounds=(bounds[0]+[30.],bounds[1]+[250.]);scale+=[100.]
    q=least_squares(lambda x:residual(x,True),seed,bounds=bounds,x_scale=scale,max_nfev=400)
    q=least_squares(residual,q.x,bounds=bounds,x_scale=scale,max_nfev=400)
    q=least_squares(lambda x:residual(x,indices=np.arange(len(f))),q.x,bounds=bounds,
                    x_scale=scale,max_nfev=400)
    s=predict(q.x,f);err=abs(s-obs)
    out=OUT/('free_feed_impedance_diagnostic' if free_impedance else 'fixed_feed_impedance_check');out.mkdir(exist_ok=True)
    result=dict(C0_T1_fF=q.x[0],Cg_T2_T1_fF=q.x[1],Cg_T1_regular_fF=q.x[2],extra_delay_ps=q.x[3],
        line_Z_ohms=float(q.x[4]) if free_impedance else 56.8874,line_impedance_fitted=free_impedance,
        fixed_line_delays_ps=delays.tolist(),port_normalization_ohms=50,
        complex_RMS=float(np.sqrt(np.mean(err**2))),complex_max=float(err.max()),
        scope='Diagnostic physical-feed approximation; not independently calibrated deembedding.')
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez(out/'response.npz',f=f,actual=obs,model=s)
    print(json.dumps(result,indent=2))
    return result


if __name__=='__main__':run()
