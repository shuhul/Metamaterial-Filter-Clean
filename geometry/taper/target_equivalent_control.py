"""Model-only guide: compensate reference values while matching target S21.

No simulation, no extraction from EM, and no geometry acceptance occurs here.
This checks whether the small fixed regular/reference differences require
different taper capacitances to reproduce the specified target response.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import numpy as np
from scipy.optimize import least_squares
from geometry.taper.target_iteration import OUT, target_response
from geometry.taper.progressive_fitting import transmission
from geometry.taper.staged_fitting import LT1, LR
from geometry.taper.corrected_iteration import REGULAR_C, REGULAR_G


def run(l2=.7347878325866108, *, fixed_unit_coupling=False, c1_floor=0):
    f=np.linspace(4.65e9,5.4e9,7501);goal=target_response(f)
    scale=np.array([100,459,100,56.7,8.2865,20])
    lo=np.array([0,.1,c1_floor,.1,.1,-100])/scale
    hi=np.array([1000,3000,1000,1000,100,100])/scale
    def model(x):
        c2,gi,c1,g1,gu,t=x*scale
        return transmission(f,[l2,LT1,LR,LR],[c2,c1,REGULAR_C,REGULAR_C],
                            [gi,g1,gu,REGULAR_G],t)
    def residual(x):
        d=model(x)-goal
        return np.r_[d.real,d.imag]
    start=np.array([100,459,100,56.7,8.2865,0])/scale
    free=[0,1,2,3,5] if fixed_unit_coupling else list(range(6))
    def expand(x):
        q=start.copy();q[free]=x
        if fixed_unit_coupling:q[4]=8.22465346074856/scale[4]
        return q
    best=least_squares(lambda x:residual(expand(x)),start[free],
                       bounds=(lo[free],hi[free]),max_nfev=300,ftol=1e-12,xtol=1e-12,gtol=1e-12)
    optimum=expand(best.x)
    result=dict(scope='Model-only compensation, not EM extraction or a tested design',
                values=dict(zip(['C0_T2_fF','Cg_T2_fF','C0_T1_fF','Cg_T1_fF','Cg_unit_fF','delay_ps'],
                                (optimum*scale).tolist())),
                fixed_L_nH=[l2,LT1,LR,LR],fixed_regular_C_fF=REGULAR_C,
                fixed_output_Cg_fF=REGULAR_G,
                target_complex_RMS=float(np.sqrt(np.mean(abs(model(optimum)-goal)**2))),
                unit_coupling_fixed=fixed_unit_coupling,T1_C0_lower_bound_fF=c1_floor,
                optimizer_success=bool(best.success),active_bounds=best.active_mask.tolist())
    OUT.mkdir(parents=True,exist_ok=True)
    name='model_only_target_equivalent'+('_fixed_unit_Cg' if fixed_unit_coupling else '')+(f'_c1floor{c1_floor}' if c1_floor else '')
    (OUT/(name+'.json')).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
    return result


if __name__=='__main__':run()
