"""Complex S21 extraction with previously measured taper-1 self C held fixed.

New taper-2 C0 and all couplings remain free. This is an additional diagnostic
of whether the independently characterized downstream stage transfers.
"""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks
from circuit_common import read_snp
from geometry.taper.fitting import chain
from geometry.taper.second_experiment import L1,LR


def fit(snp,l2,output):
    net=read_snp(snp)
    peaks,_=find_peaks(abs(net.s[:,1,0]),prominence=.02)
    peaks=[i for i in peaks if net.f[i]<6e9]
    if len(peaks)!=2:raise ValueError('Expected two lower passband peaks')
    mask=(net.f>=net.f[peaks[0]]-.35e9)&(net.f<=net.f[peaks[-1]]+.35e9)
    f=net.f[mask];observed=net.s[mask,1,0]
    # C02, Cgin2, Cg21, Cg1r, Cgout, Cr, delay ps.
    lower=[10,20,5,2,2,150,0];upper=[500,2000,200,25,25,230,40]
    def response(x):
        s=chain(2*np.pi*f,[l2,L1,LR],[x[0],101.492048434,x[5]],x[1:5])[:,1,0]
        return s*np.exp(-2j*np.pi*f*x[6]*1e-12)
    def residual(x):
        d=response(x)-observed
        return np.r_[d.real,d.imag]
    starts=[[100,c,56.7,8.2865,8.2865,cr,5] for c in (135,459,900) for cr in (185.7,188.45)]
    results=[least_squares(residual,x,bounds=(lower,upper),x_scale='jac',max_nfev=800,ftol=1e-11) for x in starts]
    r=min(results,key=lambda r:r.fun@r.fun);x=r.x
    sv=np.linalg.svd(r.jac*np.maximum(abs(x),1)[None,:],compute_uv=False)
    result={'x':x.tolist(),'L_t2_nH':l2,'C0_t1_fixed_fF':101.492048434,
            'C0_t2_fF':float(x[0]),'Cg_t2_fF':float(x[1]),'Cg_t1_fF':float(x[2]),'Cg_fF':float(x[3]),
            'Cg_output_fF':float(x[4]),'C0_regular_fF':float(x[5]),'delay_ps':float(x[6]),
            'objective':'complex S21 in passband; measured taper-1 C0 fixed, new taper-2 C0 and all couplings free',
            's21_complex_rmse':float(np.sqrt(np.mean(abs(response(x)-observed)**2))),
            'fit_range_GHz':[float(f[0]/1e9),float(f[-1]/1e9)],'success':bool(r.success),
            'active_bounds':r.active_mask.tolist(),'relative_jacobian_condition':float(sv[0]/sv[-1])}
    Path(output).write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result
