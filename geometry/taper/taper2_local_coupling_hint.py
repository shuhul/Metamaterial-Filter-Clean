"""Local coupling guidance from matched-body EM and extracted C values.

Purely offline. Neither interpolation nor its minimum accepts a design. Every
suggested geometry requires native Sonnet and the original fixed-target test.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import copy,json,itertools
import numpy as np
from scipy.optimize import least_squares,minimize_scalar
from geometry.taper import taper2_plain as w

CASES=[
 'idc_n24fg20_right6_dl256_dr42_pad148_w152_p594_i160_j40_r74tg2',
 'idc_n24fg20_right6_dl256_dr42_pad148_w152_p594_i100_j40_r74tg2',
 'idc_n24fg20_right6_dl256_dr42_pad148_w152_p594_i100_j32_r74tg2',
 'idc_n24fg20_right6_dl256_dr42_pad148_w152_p594_i160_j40_r78tg2']

def propose(cases=CASES):
    root=w.ROOT/'output_data/taper_exploration'/w.CAMPAIGN
    native=w.ROOT/'sonnet_files/taper_exploration'/w.CAMPAIGN
    dims=[];caps=[];observed=[];fixed=None;fits=[]
    for name in cases:
        settings=json.loads((native/name/'response/geometry.json').read_text())['trial_settings']
        settings=copy.deepcopy(settings);c=settings['coupling']
        dims.append([c.pop('d'),c.pop('d_t1'),c.pop('d_regular')])
        if fixed is not None and settings!=fixed:raise ValueError('Non-coupling geometry changed')
        fixed=settings
        fit=json.loads((root/name/'response/extraction/fit.json').read_text());fits.append(fit)
        if fits[0]['L_nH']!=fit['L_nH']:raise ValueError('Independent L changed')
        caps.append([*fit['C_fF'][:2],*fit['G_fF'][:3]])
        a=np.load(root/name/'response/comparison.npz')
        if observed and not np.array_equal(a['f'],f):raise ValueError('Frequency samples changed')
        f=a['f'];observed.append(a['observed']*np.exp(2j*np.pi*f*fit['delay_ps']*1e-12))
    def feature(p):
        p=np.asarray(p);return np.stack([p[...,0]/160,np.log(p[...,1]),np.log(p[...,2])],axis=-1)
    dims=np.asarray(dims,float);caps=np.asarray(caps);feat=feature(dims)
    df=feat[1:]-feat[0]
    if np.linalg.matrix_rank(df)!=3:raise ValueError('Three independent coupling directions required')
    slopes=np.linalg.lstsq(df,caps[1:]-caps[0],rcond=None)[0]
    inverse=np.linalg.lstsq(df,1/np.asarray(observed[1:])-1/observed[0],rcond=None)[0]
    ref=fits[0];ls=ref['L_nH'];full=np.arange(len(f));ix=full[::5]
    def circuit(cp,indices):
        return w.common.transmission(f[indices],ls,[*cp[:2],*ref['C_fF'][2:]],
                                     [*cp[2:],*ref['G_fF'][3:]],0)
    correction=1/observed[0]-1/circuit(caps[0],full)
    def predict(p,indices,mode):
        delta=feature(p)-feat[0]
        if mode=='native_inverse':return 1/(1/observed[0][indices]+delta@inverse[:,indices])
        cp=caps[0]+delta@slopes
        return 1/(1/circuit(cp,indices)+correction[indices])
    stages=[w.common.TAPER_2,w.common.TAPER_1,w.common.UNIT_CELL,w.common.UNIT_CELL]
    target=w.common.transmission(f,[q.L0_nH for q in stages],[q.C0_fF for q in stages],
        [q.Cg_input_fF for q in stages],0)
    report=dict(scope='Offline proposal only; no native acceptance or changed targets',sources=cases,
        dimensions=['input_gap','T2_T1_gap','T1_unit_gap'],
        capacitance_names=['C0_T2','C0_T1','Cg_T2','Cg_T1','Cg_unit'],
        feature_names=['input_gap/160','log(T2_T1_gap)','log(T1_unit_gap)'],
        capacitance_feature_slopes=slopes.tolist(),models={})
    lo=np.array([80,16,54,-100]);hi=np.array([260,120,98,100])
    for mode in ['capacitance_with_native_inverse_correction','native_inverse']:
        def residual(q):
            d=predict(q[:3],ix,mode)*np.exp(-2j*np.pi*f[ix]*q[3]*1e-12)-target[ix]
            return np.r_[d.real,d.imag]
        attempts=[least_squares(residual,np.r_[p,0],bounds=(lo,hi),
            x_scale=[100,40,74,20],max_nfev=300) for p in dims]
        best=min(attempts,key=lambda r:r.fun@r.fun)
        step=np.array([2,2,4]);offset=np.array([0,0,2])
        floor=np.floor((best.x[:3]-offset)/step)*step+offset;rows=[]
        for shift in itertools.product((0,1),repeat=3):
            p=floor+step*np.array(shift)
            if np.any(p<lo[:3]) or np.any(p>hi[:3]):continue
            curve=predict(p,full,mode)
            objective=lambda tau:float(np.mean(abs(curve*np.exp(-2j*np.pi*f*tau*1e-12)-target)**2))
            phase=minimize_scalar(objective,bounds=(-100,100),method='bounded')
            rows.append(dict(dimensions_um=p.tolist(),predicted_RMS=float(np.sqrt(phase.fun)),
                heuristic_C_fF=(caps[0]+(feature(p)-feat[0])@slopes).tolist()))
        report['models'][mode]=dict(continuous_dimensions_um=best.x[:3].tolist(),
            training_curve_RMS=[float(np.sqrt(np.mean(abs(predict(p,full,mode)-obs)**2)))
                                for p,obs in zip(dims,observed)],
            candidates=sorted(rows,key=lambda r:r['predicted_RMS']))
    (root/'sparse_coupling_proposals.json').write_text(json.dumps(report,indent=2)+'\n')
    for mode,data in report['models'].items():print(mode,data['candidates'][0],flush=True)
    return report

if __name__=='__main__':propose()
