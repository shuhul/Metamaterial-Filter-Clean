"""Suggest T2 parameters from completed EM curves; never accept a prediction.

Inverse transmission is the circuit's ABCD denominator. Interpolating it is
more useful near moving resonances than blending transmission magnitudes.
Descriptors are measured center-frequency shunt susceptance, inverse IDC
overlap, and the T2-to-T1 gap. Every suggestion needs fresh EM verification.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import json
import numpy as np
from scipy.optimize import least_squares
from geometry.taper.taper2_v3 import ROOT,CAMPAIGN,W,AP,transmission
from geometry.specifications import TAPER_1,TAPER_2,UNIT_CELL

CASES=['t2_w152_pad532_idc28_n72','t2_w152_pad24x586_d12_idc34_n72',
       't2_w152_pad532_idc28_n72_j96','t2_c264_top10_w152_pad532_idc34_n72_j112',
       't2_c264_top10_w152_pad532_idc34_n72_j112_rg64tg4',
       't2_c264_top10_w152_pad532_idc34_n72_j112_rg62tg6']


def describe(case,normalize_weak_coupling=False):
    base=ROOT/'output_data/taper_exploration'/CAMPAIGN/case
    meta=json.loads((ROOT/'geometry/taper/output'/CAMPAIGN/case/'response/geometry.json').read_text())
    net=W['read_network'](ROOT/'sonnet_files/taper_exploration'/CAMPAIGN/case/'Lcal/generated_ncell_1.snp',4)
    s,z=AP['reduce_pairs'](net.s);k,c,error=AP['extract_kc'](net.f,AP['s_to_y'](s,z))
    l=json.loads((base/'Lcal/results.json').read_text())['L0_nH']
    omega=2*np.pi*5.025e9
    effective_c=float(c[2].sum())-1/(omega**2*l*1e-9)*1e15
    coupling=meta['coupler'];area=coupling['n']*(coupling['d']-2*coupling['d_tip'])
    ground=coupling.get('t_regular_ground',6)
    gap=coupling.get('d_t1',56)
    if normalize_weak_coupling:
        # Measured baseline right-probe coupling; this is not a target value.
        gap*=6.68838984/float(-c[2,1])
    return [effective_c,1/area,gap,ground,(coupling.get('d_regular',66)-ground)/2]


def suggest(cases=CASES,ground_values=(4,6,8),clearance_values=(24,26,28,30,32),net_c_fixed=None,
            normalize_weak_coupling=False):
    base=ROOT/'output_data/taper_exploration'/CAMPAIGN
    observed=[];descriptors=[]
    for case in cases:
        data=np.load(base/case/'response/comparison.npz')
        if observed and not np.array_equal(f,data['f']):raise ValueError('Frequency grids differ')
        f=data['f'];observed.append(data['observed']);descriptors.append(describe(case,normalize_weak_coupling))
    features=np.array(descriptors);center=features.mean(axis=0);scale=np.ptp(features,axis=0)
    design=np.c_[np.ones(len(features)),(features-center)/scale]
    if np.linalg.matrix_rank(design)<6:raise ValueError('Need independent completed geometry directions')
    coefficients=np.linalg.lstsq(design,1/np.array(observed),rcond=None)[0]
    target=transmission(f,[TAPER_2.L0_nH,TAPER_1.L0_nH,UNIT_CELL.L0_nH,UNIT_CELL.L0_nH],
                        [TAPER_2.C0_fF,TAPER_1.C0_fF,UNIT_CELL.C0_fF,UNIT_CELL.C0_fF],
                        [TAPER_2.Cg_input_fF,TAPER_1.Cg_input_fF,UNIT_CELL.Cg_input_fF,UNIT_CELL.Cg_input_fF],0)
    def predict(x,ix,ground,clearance):
        descriptor=np.array([x[0],1/(72*(x[1]-4)),x[2],ground,clearance])
        row=np.r_[1,(descriptor-center)/scale]
        return 1/(row@coefficients[:,ix])*np.exp(-2j*np.pi*f[ix]*x[3]*1e-12)
    ix=np.arange(0,len(f),4)
    def residual(x,indices=ix,ground=4,clearance=30):
        delta=predict(x,indices,ground,clearance)-target[indices]
        return np.r_[delta.real,delta.imag]
    bounds=([-150,12,40,-100],[-40,70,220,100])
    seeds=[[v[0],4+1/v[1]/72,v[2],0] for v in features]+[[-100,34,100,0]]
    if net_c_fixed is not None:
        bounds[0][0]=net_c_fixed-1e-6;bounds[1][0]=net_c_fixed+1e-6
        for seed in seeds:seed[0]=net_c_fixed
    trials=[(g,c,least_squares(lambda x:residual(x,ground=g,clearance=c),s,bounds=bounds,
                              x_scale=[40,20,80,20],max_nfev=200))
            for g in ground_values for c in clearance_values for s in seeds]
    ground,clearance,best=min(trials,key=lambda item:item[2].fun@item[2].fun)
    all_indices=np.arange(len(f))
    best=least_squares(lambda x:residual(x,all_indices,ground,clearance),best.x,bounds=bounds,
                       x_scale=[40,20,80,20],max_nfev=300)
    descriptor=np.array([best.x[0],1/(72*(best.x[1]-4)),best.x[2],ground,clearance])
    weights=np.r_[1,(descriptor-center)/scale]@np.linalg.pinv(design)
    result={'method':'EM inverse-transmission local surrogate; not a simulation',
            'gap_normalized_by_measured_weak_Cg':normalize_weak_coupling,
            'fixed_body_net_C_fF':net_c_fixed,
            'cases':cases,'descriptors':descriptors,'design_condition':float(np.linalg.cond(design)),
            'center_frequency_GHz':5.025,'requested_net_C_fF':float(best.x[0]),
            'requested_weak_ground_C_fF_at_target_L':float(best.x[0]+1/((2*np.pi*5.025e9)**2*TAPER_2.L0_nH*1e-9)*1e15),
            'requested_IDC_gap_um_at_n72':float(best.x[1]),'requested_T2_T1_gap_um':float(best.x[2]),
            'requested_regular_ground_um':ground,'requested_regular_gap_um':ground+2*clearance,
            'requested_regular_ground_clearance_um':clearance,
            'native_mesh_manufacturable':bool(ground%2==0 and clearance%2==0),
            'delay_ps':float(best.x[3]),'source_weights':weights.tolist(),
            'predicted_target_RMS':float(np.sqrt(np.mean(abs(predict(best.x,all_indices,ground,clearance)-target)**2))),
            'discrete_ground_candidates':[{'ground_um':g,'clearance_um':c,
                 'predicted_RMS':min(float(np.sqrt(q.fun@q.fun/len(ix))) for gg,cc,q in trials if gg==g and cc==c)}
                 for g in ground_values for c in clearance_values],
            'requires_fresh_EM_validation':True}
    (base/'response_surrogate_proposal.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True);return result


def validate_prediction(held_case,cases=CASES,normalize_weak_coupling=False):
    """Check a genuinely held-out EM curve before trusting another proposal."""
    base=ROOT/'output_data/taper_exploration'/CAMPAIGN
    features=np.array([describe(case,normalize_weak_coupling) for case in cases])
    center=features.mean(axis=0);scale=np.ptp(features,axis=0)
    design=np.c_[np.ones(len(features)),(features-center)/scale]
    curves=[np.load(base/case/'response/comparison.npz')['observed'] for case in cases]
    coefficients=np.linalg.lstsq(design,1/np.array(curves),rcond=None)[0]
    point=np.array(describe(held_case,normalize_weak_coupling))
    predicted=1/(np.r_[1,(point-center)/scale]@coefficients)
    actual=np.load(base/held_case/'response/comparison.npz')['observed']
    result={'held_out_case':held_case,'training_cases':cases,
            'gap_normalized_by_measured_weak_Cg':normalize_weak_coupling,
            'prediction_complex_RMS':float(np.sqrt(np.mean(abs(predicted-actual)**2))),
            'maximum_prediction_error':float(abs(predicted-actual).max()),
            'scope':'Surrogate validation only, not target acceptance'}
    label='weak_normalized' if normalize_weak_coupling else 'raw'
    (base/f'surrogate_validation_{label}.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True);return result


if __name__=='__main__':suggest()
