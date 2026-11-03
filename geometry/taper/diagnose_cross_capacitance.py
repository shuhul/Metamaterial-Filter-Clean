"""Offline test of unmodelled direct input coupling; not a production fit."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import numpy as np
from scipy.optimize import least_squares
from circuit_common import read_snp
from geometry.taper.locked_taper2 import ROOT
from geometry.taper.diagnose_mutual_inductance import coupled_response


def diagnose(name):
    out=ROOT/'output_data/taper_exploration'/name
    trial=json.loads((out/'trial_locked.json').read_text());net=read_snp(out/'response.snp')
    f=net.f;l2=trial['fit']['L_t2_nH'];observed=net.s[:,1,0];results={}
    for node,upper in ((2,10),(3,1)):
        def predict(x):return coupled_response(f,[x[0],x[1],0],l2,cross_capacitance=(0,node,x[2]))[:,1,0]
        def residual(x):
            d=predict(x)-observed;return np.r_[d.real,d.imag]
        grid=[np.array([c,g,cross]) for c in (10,100,250) for g in (100,300,459,800)
              for cross in np.linspace(0,upper,41)]
        grid.sort(key=lambda x:float(residual(x)@residual(x)))
        fits=[least_squares(residual,x,bounds=([0,30,0],[400,1500,upper]),
                            x_scale=[100,459,upper/5],max_nfev=350)
              for x in grid[:10]+[np.r_[trial['fit']['x'],0]]]
        best=min(fits,key=lambda r:r.fun@r.fun);d=predict(best.x)-observed
        results[f'input_to_node_{node}']={'x':best.x.tolist(),'active_bounds':best.active_mask.tolist(),
                'complex_rmse':float(np.sqrt(np.mean(abs(d)**2))),'max_complex_error':float(max(abs(d)))}
    result={'purpose':'Three-parameter direct-cross-coupling hypothesis ONLY; not accepted extraction',
            'free_parameters':['C0_t2_fF','Cg_t2_fF','cross_capacitance_fF'],'case':name,'hypotheses':results}
    (out/'cross_capacitance_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True);return result


if __name__=='__main__':diagnose('compact_idc30_span172_top148_w156_ground80_rise370_right536')
