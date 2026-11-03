"""Offline check of a measured input-fixture shunt; NOT a promoted extraction.

The 84.31 fF comes from the separate low-frequency IDC diagnostic. Its terminal
surroundings differ, so this tests a hypothesis rather than establishing a new
production fixture calibration. Exactly two capacitances vary in each fit.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import numpy as np
from scipy.optimize import least_squares
from geometry.taper.locked_taper2 import ROOT,FIXED,response
from geometry.taper.fitting import chain
from geometry.taper.check_reference_transfer import s_to_abcd,abcd_to_s
from circuit_common import read_snp


def diagnose(name):
    out=ROOT/'output_data/taper_exploration'/name
    trial=json.loads((out/'trial_locked.json').read_text())
    net=read_snp(out/'response.snp');mask=(net.f>=4.5e9)&(net.f<=5.5e9)
    f=net.f[mask];obs=net.s[mask];w=2*np.pi*f;l2=trial['fit']['L_t2_nH']
    idc=read_snp(ROOT/'sonnet_files/taper_exploration/idc_check_pairs30_span172_t4_d4/generated_ncell_1.snp')
    identity=np.tile(np.eye(2),(len(idc.f),1,1))
    y=(identity-idc.s)@np.linalg.inv(identity+idc.s)/50
    c_input=float(np.mean((y[:,0,0]+y[:,0,1]).imag/(2*np.pi*idc.f)*1e15))
    results={}
    for fixed_input in (0,c_input):
        def prediction(x):
            s=chain(w,[l2,FIXED['L_taper_nH'],FIXED['L_regular_nH']],
                    [x[0],FIXED['C_taper_fF'],FIXED['C_regular_fF']],
                    [x[1],FIXED['Cg_input_fF'],FIXED['Cg_junction_fF'],FIXED['Cg_output_fF']])
            a=np.tile(np.eye(2,dtype=complex),(len(f),1,1));a[:,1,0]=1j*w*fixed_input*1e-15
            s=abcd_to_s(a@s_to_abcd(s))
            left=FIXED['fixture_delay_left_ps']*1e-12;right=FIXED['fixture_delay_right_ps']*1e-12
            s[:,0,0]*=np.exp(-2j*w*left);s[:,1,1]*=np.exp(-2j*w*right)
            s[:,1,0]*=np.exp(-1j*w*(left+right));s[:,0,1]=s[:,1,0]
            return s
        if fixed_input==0 and not np.allclose(prediction([100,459]),response(f,[100,459],l2),atol=1e-8):
            raise ValueError('Zero-shunt model does not recover production response')
        for kind in ('S21','joint_S21_S11'):
            def residual(x):
                d=prediction(x)-obs
                z=d[:,1,0] if kind=='S21' else np.r_[d[:,1,0],d[:,0,0]]
                return np.r_[z.real,z.imag]
            grid=[np.array([c,g]) for c in np.linspace(10,300,21) for g in np.linspace(100,1000,21)]
            grid.sort(key=lambda x:float(residual(x)@residual(x)))
            starts=grid[:8]+[np.array([100,459]),np.array(trial['fit']['x'])]
            fits=[least_squares(residual,x,bounds=([0,30],[400,1500]),x_scale='jac',max_nfev=300) for x in starts]
            best=min(fits,key=lambda r:r.fun@r.fun);pred=prediction(best.x)
            results[f'input_C_{fixed_input:.3f}fF_{kind}']={
                'C0_t2_fF':float(best.x[0]),'Cg_t2_fF':float(best.x[1]),
                'active_bounds':best.active_mask.tolist(),
                'S21_complex_rmse':float(np.sqrt(np.mean(abs(pred[:,1,0]-obs[:,1,0])**2))),
                'S21_max_complex_error':float(np.max(abs(pred[:,1,0]-obs[:,1,0]))),
                'S11_complex_rmse':float(np.sqrt(np.mean(abs(pred[:,0,0]-obs[:,0,0])**2)))}
    result={'purpose':'Offline fixed-input-fixture hypothesis; NOT accepted design or final extraction',
            'case':name,'measured_standalone_input_ground_fF':c_input,
            'fixed_T1_and_regular':True,'free_parameters':['C0_t2_fF','Cg_t2_fF'],'fits':results}
    (out/'input_fixture_diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result


if __name__=='__main__':
    from geometry.taper.run_locked_taper2 import CASE_NAME
    diagnose(CASE_NAME)
