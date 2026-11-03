"""Six independent native controls around the best body: proposal only.

The exact source curves provide local geometry sensitivities. Capacitance
derivatives are included as heuristics, never independent measurements.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import copy,itertools,json
import numpy as np
from scipy.optimize import least_squares,minimize_scalar
from geometry.taper import taper2_plain as w
from geometry.taper.taper2_cleared_tuning import peak_diagnostics

CASES=[
 'idc_rail8_right30_dl256_pad148_w148_p200_i42n70_j42_r74tg2',
 'idc_rail8_right30_dl256_pad148_w148_p200_i52n70_j42_r74tg2',
 'idc_rail8_right30_dl256_pad148_w148_p400_i42n70_j42_r74tg2',
 'idc_rail8_dl256_pad148_w148_p200_i42n70_j42_r74tg2',
 'idc_rail8_right30_dl256_pad148_w148_p200_i42n70_j42_r74tg6',
 'idc_rail8_right30_dl256_pad148_w148_p200_i42n70_j50_r74tg2',
 'idc_rail8_right30_dl256_pad148_w148_p200_i42n70_j42_r78tg2']

def propose(cases=CASES):
    out=w.ROOT/'output_data/taper_exploration'/w.CAMPAIGN
    native=w.ROOT/'sonnet_files/taper_exploration'/w.CAMPAIGN
    dims=[];curves=[];caps=[];peaks=[];fixed=None
    for name in cases:
        s=copy.deepcopy(json.loads((native/name/'response/geometry.json').read_text())['trial_settings'])
        dims.append([s['cap'].pop('t_right'),s['pad'].pop('l'),
            *[s['coupling'].pop(k) for k in ['d','d_t1','d_regular','t_regular_ground']]])
        if fixed is not None and s!=fixed:raise ValueError('Uncontrolled geometry changed')
        fixed=s
        a=np.load(out/name/'response/comparison.npz')
        if curves and not np.array_equal(a['f'],f):raise ValueError('Changed frequency samples')
        f=a['f'];tau=json.loads((out/name/'response/target_comparison.json').read_text())['fixture_delay_ps']
        curves.append(a['observed']*np.exp(2j*np.pi*f*tau*1e-12))
        p=peak_diagnostics(f,a['observed'])
        if len(p)!=2:raise ValueError('Two peaks required for the local peak diagnostic')
        peaks.append([v for q in p for v in [q['frequency_GHz']*1000,
            np.log(q['magnitude']),np.log(q['half_prominence_power_width_MHz'])]])
        fit=json.loads((out/name/'response/extraction/fit.json').read_text())
        caps.append([*fit['C_fF'][:2],*fit['G_fF'][:3]])
    def features(p):
        p=np.asarray(p)
        return np.stack([np.log(p[...,0]),p[...,1]/200,p[...,2]/42,
            np.log(p[...,3]),np.log(p[...,4]),p[...,5]/6],axis=-1)
    dims=np.array(dims,float);feat=features(dims);delta=feat[1:]-feat[0]
    if np.linalg.matrix_rank(delta)!=6:raise ValueError('Six independent controls required')
    curves=np.array(curves);inverse=1/curves;caps=np.array(caps);peaks=np.array(peaks)
    coef=np.linalg.lstsq(delta,inverse[1:]-inverse[0],rcond=None)[0]
    cpcoef=np.linalg.lstsq(delta,caps[1:]-caps[0],rcond=None)[0]
    pkcoef=np.linalg.lstsq(delta,peaks[1:]-peaks[0],rcond=None)[0]
    stages=[w.common.TAPER_2,w.common.TAPER_1,w.common.UNIT_CELL,w.common.UNIT_CELL]
    target=w.common.transmission(f,[q.L0_nH for q in stages],[q.C0_fF for q in stages],
        [q.Cg_input_fF for q in stages],0)
    target_peaks=peak_diagnostics(f,target)
    pkgoal=np.array([v for q in target_peaks for v in [q['frequency_GHz']*1000,
        np.log(q['magnitude']),np.log(q['half_prominence_power_width_MHz'])]])
    full=np.arange(len(f));ix=full[::5]
    def predict(p,indices):return 1/(inverse[0,indices]+(features(p)-feat[0])@coef[:,indices])
    lo=np.array([6,2,24,20,54,2,-100]);hi=np.array([62,594,80,120,98,10,100])
    reports={}
    for mode in ['native_complex_response','peak_features']:
        def residual(q):
            if mode=='peak_features':
                return np.r_[(peaks[0]+(features(q[:6])-feat[0])@pkcoef-pkgoal)/
                             [5,.15,.2,5,.15,.2],q[6]*1e-5]
            d=predict(q[:6],ix)*np.exp(-2j*np.pi*f[ix]*q[6]*1e-12)-target[ix]
            return np.r_[d.real,d.imag]
        fits=[least_squares(residual,np.r_[p,0],bounds=(lo,hi),
            x_scale=[30,200,42,42,74,6,20],max_nfev=300) for p in dims]
        best=min(fits,key=lambda q:q.fun@q.fun)
        step=np.array([4,2,2,2,2,2]);offset=np.array([2,0,0,0,0,0])
        floor=np.floor((best.x[:6]-offset)/step)*step+offset;rows=[]
        for shift in itertools.product((0,1),repeat=6):
            p=floor+step*np.array(shift)
            if np.any(p<lo[:6]) or np.any(p>hi[:6]) or (p[4]-p[5])%4:continue
            curve=predict(p,full)
            objective=lambda tau:float(np.mean(abs(curve*np.exp(-2j*np.pi*f*tau*1e-12)-target)**2))
            phase=minimize_scalar(objective,bounds=(-100,100),method='bounded')
            rows.append(dict(dimensions_um=p.tolist(),predicted_complex_RMS=float(np.sqrt(phase.fun)),
                heuristic_C_fF=(caps[0]+(features(p)-feat[0])@cpcoef).tolist(),
                predicted_peak_features=(peaks[0]+(features(p)-feat[0])@pkcoef).tolist()))
        reports[mode]=dict(continuous_dimensions_um=best.x[:6].tolist(),
            candidates=sorted(rows,key=lambda q:q['predicted_complex_RMS']))
    report=dict(scope='Local measured-control proposals only; no native acceptance',sources=cases,
        dimensions=['right_leg_t','pad_l','input_gap','T2_T1_gap','T1_unit_gap','ground_strip_t'],
        cap_names=['C0_T2','C0_T1','Cg_T2','Cg_T1','Cg_unit'],
        capacitance_feature_slopes=cpcoef.tolist(),models=reports)
    (out/'best_controlled_proposals.json').write_text(json.dumps(report,indent=2)+'\n')
    for mode,r in reports.items():print(mode,r['candidates'][:2],flush=True)
    return report

if __name__=='__main__':propose()
