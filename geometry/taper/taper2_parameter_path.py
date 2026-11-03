"""Predict a next trial from two completed T2 parameter variants, without EM.

Interpolate inverse complex transmission (the circuit denominator), rather
than mixing peak magnitudes. This is a local heuristic: the proposed geometry
must still pass a fresh full Sonnet simulation before it is accepted.
"""
import json
import numpy as np
from scipy.optimize import least_squares
from geometry.taper.taper2_v3 import ROOT,CAMPAIGN,transmission
from geometry.specifications import TAPER_1,TAPER_2,UNIT_CELL


def suggest(low='t2_w152_pad532_idc28_n72',high='t2_w152_pad24x586_d12_idc34_n72'):
    base=ROOT/'output_data/taper_exploration'/CAMPAIGN
    a=np.load(base/low/'response/comparison.npz');b=np.load(base/high/'response/comparison.npz')
    f=a['f']
    if not np.array_equal(f,b['f']):raise ValueError('Frequency grids differ')
    first=a['observed'];last=b['observed']
    target=transmission(f,[TAPER_2.L0_nH,TAPER_1.L0_nH,UNIT_CELL.L0_nH,UNIT_CELL.L0_nH],
                        [TAPER_2.C0_fF,TAPER_1.C0_fF,UNIT_CELL.C0_fF,UNIT_CELL.C0_fF],
                        [TAPER_2.Cg_input_fF,TAPER_1.Cg_input_fF,UNIT_CELL.Cg_input_fF,UNIT_CELL.Cg_input_fF],0)
    def predict(x):
        return 1/((1-x[0])/first+x[0]/last)*np.exp(-2j*np.pi*f*x[1]*1e-12)
    def residual(x):
        delta=predict(x)-target
        return np.r_[delta.real,delta.imag]
    fits=[least_squares(residual,[u,0],bounds=([-1,-100],[2,100]),x_scale=[1,20],max_nfev=200)
          for u in (-1,0,.5,1,2)]
    best=min(fits,key=lambda q:q.fun@q.fun)
    result={'method':'Local interpolation of inverse complex S21; proposal only',
            'low':low,'high':high,'high_fraction':float(best.x[0]),'delay_ps':float(best.x[1]),
            'predicted_target_RMS':float(np.sqrt(np.mean(abs(predict(best.x)-target)**2))),
            'requires_fresh_EM_validation':True}
    (base/'parameter_path_proposal.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2));return result


if __name__=='__main__':suggest()
