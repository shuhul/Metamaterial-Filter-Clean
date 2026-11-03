"""Offline two-C fit with independently measured IDC fixture parasitics.

Hypothesis only: standalone IDC surroundings differ from the actual fixture.
Its three low-frequency samples determine a fixed C matrix and series L matrix,
not additional fit parameters. The node-side static ground C is subtracted
before adding C0_t2, so the latter still means TOTAL node-to-ground capacitance.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import numpy as np
from scipy.optimize import least_squares
from circuit_common import read_snp
from geometry.taper.locked_taper2 import ROOT,FIXED
from geometry.taper.diagnose_locked_response import downstream_abcd
from geometry.taper.check_reference_transfer import abcd_to_s


def fixture_matrices():
    net=read_snp(ROOT/'sonnet_files/taper_exploration/idc_check_pairs30_span172_t4_d4/generated_ncell_1.snp')
    eye=np.tile(np.eye(2),(len(net.f),1,1))
    z=50*(eye+net.s)@np.linalg.inv(eye-net.s)
    w=2*np.pi*net.f
    # Scale the columns to avoid a conditioning problem in SI units.
    basis=np.c_[-1/(w/1e9),w/1e9]
    coefficients=np.linalg.lstsq(basis,z.imag.reshape(len(w),4),rcond=None)[0]
    invc=coefficients[0].reshape(2,2)*1e9
    inductance=coefficients[1].reshape(2,2)*1e-9
    capacitance=np.linalg.inv(invc)
    reconstructed=-invc[None,:,:]/w[:,None,None]+w[:,None,None]*inductance
    return capacitance,inductance,float(np.max(abs(z.imag-reconstructed)))


def diagnose(name):
    out=ROOT/'output_data/taper_exploration'/name
    trial=json.loads((out/'trial_locked.json').read_text());n=read_snp(out/'response.snp')
    mask=(n.f>=4.5e9)&(n.f<=5.5e9);f=n.f[mask];w=2*np.pi*f;obs=n.s[mask,1,0]
    l2=trial['fit']['L_t2_nH'];cap,inductance,calerr=fixture_matrices()
    ground_in=float(cap[0].sum());ground_out=float(cap[1].sum())
    suffix=downstream_abcd(f);results={}
    for left_delay in (0,FIXED['fixture_delay_left_ps']):
        def predict(x):
            mutual=x[1]*1e-15
            c=np.array([[ground_in+mutual,-mutual],[-mutual,ground_out+mutual]])
            z=np.linalg.inv(c)[None,:,:]/(1j*w[:,None,None])+1j*w[:,None,None]*inductance
            a=np.empty_like(z);a[:,0,0]=z[:,0,0]/z[:,1,0]
            a[:,0,1]=(z[:,0,0]*z[:,1,1]-z[:,0,1]*z[:,1,0])/z[:,1,0]
            a[:,1,0]=1/z[:,1,0];a[:,1,1]=z[:,1,1]/z[:,1,0]
            shunt=np.tile(np.eye(2,dtype=complex),(len(f),1,1))
            shunt[:,1,0]=1j*w*(x[0]*1e-15-ground_out)+1/(1j*w*l2*1e-9)
            s=abcd_to_s(a@shunt@suffix)
            return s[:,1,0]*np.exp(-1j*w*(left_delay+FIXED['fixture_delay_right_ps'])*1e-12)
        def residual(x):
            d=predict(x)-obs;return np.r_[d.real,d.imag]
        grid=[np.array([c,g]) for c in np.linspace(0,400,21) for g in np.linspace(100,1000,21)]
        grid.sort(key=lambda x:float(residual(x)@residual(x)))
        fits=[least_squares(residual,x,bounds=([0,30],[500,1500]),x_scale='jac',max_nfev=300)
              for x in grid[:8]+[np.array([100,459])]]
        best=min(fits,key=lambda r:r.fun@r.fun);d=predict(best.x)-obs
        results[str(left_delay)]={'C0_t2_fF':float(best.x[0]),'Cg_t2_fF':float(best.x[1]),
                                 'active_bounds':best.active_mask.tolist(),
                                 'complex_rmse':float(np.sqrt(np.mean(abs(d)**2))),
                                 'max_complex_error':float(max(abs(d)))}
    result={'purpose':'Independent standalone IDC parasitic hypothesis only; NOT a final extraction',
            'case':name,'free_parameters':['C0_t2_fF','Cg_t2_fF'],
            'fixed_capacitance_matrix_fF':(cap*1e15).tolist(),
            'fixed_series_inductance_matrix_nH':(inductance*1e9).tolist(),
            'low_frequency_Z_fit_max_error_ohm':calerr,
            'input_ground_fF':ground_in*1e15,'node_ground_fF':ground_out*1e15,
            'fits_by_fixed_extra_left_delay_ps':results}
    (out/'idc_fixture_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result


if __name__=='__main__':
    diagnose('compact_idc30_span172_top148_w156_ground80_rise370_right536')
