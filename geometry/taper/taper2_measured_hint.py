"""Interpolate actual complex responses to propose nearby dimensions only.

No EM launches or acceptance decisions. A low predicted error must be checked
by a fresh native simulation. Responses must use the same frequency samples.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
import itertools
import copy
import numpy as np
from scipy.optimize import least_squares, minimize_scalar, lsq_linear
from geometry.taper import taper2_plain as w


def controlled_family(cases):
    """Phase-aligned actual-curve guidance for leg/pad/IDC dimensions only.

    Electrostatic proxies capture diminishing leg-width effects and the pad
    tip approaching the cap bar. This is a proposal model, never acceptance.
    """
    root=w.ROOT/'output_data/taper_exploration'/w.CAMPAIGN
    native=w.ROOT/'sonnet_files/taper_exploration'/w.CAMPAIGN
    dimensions=[];curves=[];fixed=None
    for case in cases:
        s=json.loads((native/case/'response/geometry.json').read_text())['trial_settings']
        unchanged=copy.deepcopy(s)
        if unchanged['coupling'].pop('w_input_stub',None) is not None:
            raise ValueError('Keep the original feed span for this controlled family')
        dimensions.append([s['cap']['t_right'],s['pad']['w'],s['pad']['l'],s['coupling']['d'],s['coupling']['d_t1'],s['coupling']['d_regular']])
        unchanged['cap'].pop('t_right');unchanged['pad'].pop('w');unchanged['pad'].pop('l')
        unchanged['coupling'].pop('d')
        unchanged['coupling'].pop('d_t1')
        unchanged['coupling'].pop('d_regular')
        if fixed is not None and unchanged!=fixed:raise ValueError('Unmodeled geometry changed')
        fixed=unchanged
        a=np.load(root/case/'response/comparison.npz')
        if curves and not np.array_equal(a['f'],f):raise ValueError('Frequency grid changed')
        f=a['f'];target=a['target']
        tau=json.loads((root/case/'response/target_comparison.json').read_text())['fixture_delay_ps']
        curves.append(a['observed']*np.exp(2j*np.pi*f*tau*1e-12))
    # All targets have identical magnitudes but separately matched delay. Use
    # the nominal zero-delay target for every aligned source and prediction.
    specs=[w.common.TAPER_2,w.common.TAPER_1,w.common.UNIT_CELL,w.common.UNIT_CELL]
    target=w.common.transmission(f,[q.L0_nH for q in specs],[q.C0_fF for q in specs],
        [q.Cg_input_fF for q in specs],0)
    def features(p):
        r,pw,pl,di,junction,regular=np.asarray(p).T
        return np.array([np.log(r),pl/596,pw*pl/(8*596),pw*np.log(596/(596-pl))/8,di/42,np.log(junction),np.log(regular)]).T
    dims=np.array(dimensions,float);feat=features(dims);delta=feat[1:]-feat[0]
    if np.linalg.matrix_rank(delta)<7:raise ValueError('Need seven independent measured feature directions')
    inverse=1/np.array(curves);coef=np.linalg.lstsq(delta,inverse[1:]-inverse[0],rcond=None)[0]
    def predict(p,indices):
        return 1/(inverse[0,indices]+(features(p)-feat[0])@coef[:,indices])
    full=np.arange(len(f));ix=full[::5]
    errors=[float(np.sqrt(np.mean(abs(predict(p,full)-actual)**2))) for p,actual in zip(dims,curves)]
    lower=[6,8,100,20,8,34,-100];upper=[94,16,594,70,100,110,100]
    def residual(q):
        d=predict(q[:6],ix)*np.exp(-2j*np.pi*f[ix]*q[6]*1e-12)-target[ix]
        return np.r_[d.real,d.imag]
    results=[least_squares(residual,np.r_[p,0],bounds=(lower,upper),x_scale=[30,8,400,42,42,74,20],
        max_nfev=300) for p in dims]
    best=min(results,key=lambda q:q.fun@q.fun)
    step=np.array([4,2,2,2,2,4]);offset=np.array([2,0,0,0,0,2])
    floor=np.floor((best.x[:6]-offset)/step)*step+offset;rows=[]
    for shift in itertools.product((0,1),repeat=6):
        p=floor+step*np.array(shift)
        if np.any(p<np.array(lower[:6])) or np.any(p>np.array(upper[:6])):continue
        curve=predict(p,full)
        objective=lambda t:float(np.mean(abs(curve*np.exp(-2j*np.pi*f*t*1e-12)-target)**2))
        phase=minimize_scalar(objective,bounds=(-100,100),method='bounded')
        rows.append(dict(dimensions_um=p.tolist(),predicted_RMS=float(np.sqrt(phase.fun))))
    report=dict(scope='Proposal only: phase-aligned actual native responses; fresh native EM required',
        cases=cases,dimensions=['right_leg_t','pad_w','pad_l','input_gap','T2_T1_gap','T1_unit_gap'],
        training_curve_RMS=errors,continuous_dimensions_um=best.x[:6].tolist(),
        candidates=sorted(rows,key=lambda q:q['predicted_RMS']))
    # Report the same complex RMS convention used by actual acceptance.
    report['predicted_continuous_RMS']=float(np.sqrt(np.mean(abs(
        predict(best.x[:6],full)*np.exp(-2j*np.pi*f*best.x[6]*1e-12)-target)**2)))
    (root/'controlled_family_proposals.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2));return report


def propose(cases, *, fixed_body=None):
    root=w.ROOT/'output_data/taper_exploration'/w.CAMPAIGN
    native=w.ROOT/'sonnet_files/taper_exploration'/w.CAMPAIGN
    arrays=[];dimensions=[];fixed=None
    for case in cases:
        data=np.load(root/case/'response/comparison.npz')
        if arrays and not np.array_equal(data['f'],f):raise ValueError('Frequency grid changed')
        f=data['f'];arrays.append(data['observed'])
        s=json.loads((native/case/'response/geometry.json').read_text())['trial_settings']
        unchanged=copy.deepcopy(s)
        unchanged['meander'].pop('w');unchanged['pad'].pop('l')
        for key in ('d','d_t1','d_regular','t_regular_ground'):unchanged['coupling'].pop(key)
        if fixed is not None and unchanged!=fixed:raise ValueError('An unmodeled geometry parameter changed')
        fixed=unchanged
        cp=s['coupling'];dimensions.append([s['meander']['w'],s['pad']['l'],cp['d'],cp['d_t1'],cp['d_regular'],cp['t_regular_ground']])
    dimensions=np.array(dimensions,dtype=float);inverse=1/np.array(arrays)
    delta=inverse[1:]-inverse[0];dd=dimensions[1:]-dimensions[0]
    projector=dd@np.linalg.pinv(dd)
    specs=[w.common.TAPER_2,w.common.TAPER_1,w.common.UNIT_CELL,w.common.UNIT_CELL]
    target=w.common.transmission(f,[s.L0_nH for s in specs],[s.C0_fF for s in specs],
        [s.Cg_input_fF for s in specs],0)
    idx=np.arange(0,len(f),5)
    def residual(q):
        weights=q[:-1]@projector;tau=q[-1]
        curve=1/(inverse[0,idx]+weights@delta[:,idx])*np.exp(-2j*np.pi*f[idx]*tau*1e-12)
        r=curve-target[idx]
        dims=dimensions[0]+weights@dd
        # Keep suggestions in a conservative dimensional neighborhood. These
        # penalties constrain a proposal, never the actual acceptance metric.
        limits=np.maximum([148,200,8,2,30,2]-dims,0)+np.maximum(dims-[152,586,56,140,98,12],0)
        body_error=np.array([]) if fixed_body is None else 100*(dims[:2]-fixed_body)
        # With more measured curves than dimensions, prohibit changing the
        # response through weights that leave every geometry dimension fixed.
        return np.r_[r.real,r.imag,limits,body_error,.01*(q[:-1]-weights)]
    best=None
    seeds=[np.zeros(len(cases)-1),np.ones(len(cases)-1)/(len(cases)-1),*np.eye(len(cases)-1)]
    for seed in seeds:
        result=least_squares(residual,np.r_[seed,0],bounds=(np.r_[np.full(len(seed),-2),-100],np.r_[np.full(len(seed),3),100]),max_nfev=300)
        if best is None or np.linalg.norm(result.fun)<np.linalg.norm(best.fun):best=result
    weights=best.x[:-1]@projector;tau=best.x[-1]
    curve=1/(inverse[0]+weights@delta)*np.exp(-2j*np.pi*f*tau*1e-12)
    report=dict(scope='Measured-curve local interpolation only; fresh native EM required',sources=cases,
        dimensions=['meander_w','pad_l','input_rail_gap','T2_T1_gap','T1_unit_gap','ground_strip_w'],
        predicted_dimensions_um=(dimensions[0]+weights@dd).tolist(),weights=weights.tolist(),
        predicted_RMS=float(np.sqrt(np.mean(abs(curve-target)**2))),delay_ps=float(tau),
        response_weight_nullspace_removed=True,fixed_body_dimensions_um=fixed_body)
    # Native-grid candidates are available once the measured directions span
    # all six dimensions. They remain predictions, never pass/fail results.
    rank=int(np.linalg.matrix_rank(dd));report['measured_dimension_rank']=rank
    if rank==6:
        continuous=dimensions[0]+weights@dd;steps=np.array([4,2,2,2,2,2]);offset=np.zeros(6)
        if fixed_body is not None:continuous[:2]=fixed_body
        floor=np.floor((continuous-offset)/steps)*steps+offset
        rows=[]
        for shift in itertools.product((0,1),repeat=6):
            dims=floor+steps*np.array(shift)
            if np.any(dims<[148,200,8,2,30,2]) or np.any(dims>[152,586,56,140,98,12]):continue
            if fixed_body is not None and not np.array_equal(dims[:2],fixed_body):continue
            if (dims[4]-dims[5])%4:continue
            wt=(dims-dimensions[0])@np.linalg.pinv(dd)
            predicted=1/(inverse[0]+wt@delta)
            objective=lambda delay:float(np.mean(abs(predicted*np.exp(-2j*np.pi*f*delay*1e-12)-target)**2))
            matched=minimize_scalar(objective,bounds=(-100,100),method='bounded')
            rows.append(dict(dimensions_um=dims.tolist(),predicted_RMS=float(np.sqrt(matched.fun)),delay_ps=float(matched.x)))
        report['grid_candidates']=sorted(rows,key=lambda r:r['predicted_RMS'])[:8]
    label='measured_curve_proposal' if fixed_body is None else 'measured_curve_proposal_fixed_body'
    (root/(label+'.json')).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2));return report


def peak_hint():
    """Local peak sensitivities from the measured, independent directions."""
    root=w.ROOT/'output_data/taper_exploration'/w.CAMPAIGN
    proposal=json.loads((root/'measured_curve_proposal.json').read_text())
    ledger=json.loads((root/'thin_rail_trials.json').read_text())
    lookup={r['case']:r for r in ledger['rows']};dims=[];peaks=[]
    for case in proposal['sources']:
        row=lookup[case];s=row['settings'];cp=s['coupling']
        dims.append([s['meander']['w'],s['pad']['l'],cp['d'],cp['d_t1'],cp['d_regular'],cp['t_regular_ground']])
        if len(row['EM_peaks'])!=2:raise ValueError('Two peaks required for this local diagnostic')
        peaks.append([v for p in row['EM_peaks'] for v in
                      (p['frequency_GHz']*1000,p['magnitude'],p['half_prominence_power_width_MHz'])])
    dims=np.array(dims);peaks=np.array(peaks);dd=dims[1:]-dims[0]
    if np.linalg.matrix_rank(dd)!=6:raise ValueError('Six independent measured dimensions required')
    jac=np.linalg.pinv(dd)@(peaks[1:]-peaks[0])
    target=np.array([v for p in row['target_peaks'] for v in
                     (p['frequency_GHz']*1000,p['magnitude'],p['half_prominence_power_width_MHz'])])
    scale=np.array([5,.1,5,5,.1,5])
    bounds=(np.array([148,16,8,2,30,2])-dims[0],np.array([152,586,80,140,98,12])-dims[0])
    fit=lsq_linear(jac.T/scale[:,None],(target-peaks[0])/scale,bounds=bounds)
    report=dict(scope='Local linear peak hint; not an EM prediction or acceptance',sources=proposal['sources'],
        dimensions=proposal['dimensions'],features=['lower_frequency_MHz','lower_height','lower_power_width_MHz',
        'upper_frequency_MHz','upper_height','upper_power_width_MHz'],
        feature_change_per_um=jac.tolist(),suggested_dimensions_um=(dims[0]+fit.x).tolist(),
        predicted_peaks=(peaks[0]+fit.x@jac).tolist(),target_peaks=target.tolist())
    (root/'peak_dimension_hint.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2));return report


if __name__=='__main__':
    propose(['idc_rail8_dl256_pad148_w152_i40n70_j36_r70tg2_preview',
             'idc_rail8_dl256_pad148_w152_i34n70_j32_r78tg2',
             'idc_rail8_dl256_pad148_w152_i32n70_j42_r78tg2'])
