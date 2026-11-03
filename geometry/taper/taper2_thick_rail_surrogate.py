"""Local proposals from identical-topology thick-rail EM experiments.

This is not an acceptance model. Every proposal needs a new full Sonnet sweep.
Unlike the earlier broad surrogate, this keeps all cap/coil/rail shapes fixed
and uses only input length, interstage spacing and existing ground-pad length.
"""
import json
import numpy as np
from scipy.optimize import least_squares
from geometry.taper.taper2_v3 import ROOT,CAMPAIGN,W,AP,transmission
from geometry.specifications import TAPER_1,TAPER_2,UNIT_CELL

CASES=['thick64_c344_w152_i30n70_j48_r54',
       'thick64_c344_w152_i12n70_j60_r54',
       'thick64_p200_w152_i16n70_j68_r54',
       'thick64_p200_w152_i28n70_j68_r54']


def describe(name):
    meta=json.loads((ROOT/'geometry/taper/output'/CAMPAIGN/name/'response/geometry.json').read_text())
    coupling=meta['coupler']
    if (coupling['n']!=70 or coupling.get('d_regular',66)!=54 or
        coupling.get('t_regular_ground',6)!=6):raise ValueError('Wrong surrogate geometry family')
    body=meta['unit']
    if (body['cap']['w']!=344 or body['cap']['t_side']!=64 or
        body['meander']['w']!=152 or body['pad']['w']!=8 or body['d_pad']!=20):
        raise ValueError('Wrong surrogate body')
    net=W['read_network'](ROOT/'sonnet_files/taper_exploration'/CAMPAIGN/name/'Lcal/generated_ncell_1.snp',4)
    s,z=AP['reduce_pairs'](net.s);k,c,error=AP['extract_kc'](net.f,AP['s_to_y'](s,z))
    measured=json.loads((ROOT/'output_data/taper_exploration'/CAMPAIGN/name/'Lcal/results.json').read_text())
    effective_c=float(c[2].sum())-1/((2*np.pi*5.025e9)**2*measured['L0_nH']*1e-9)*1e15
    # Measured 2-um plain proxy and 28-um IDC proxy give a screening scale.
    # This mapping is approximate; the input-only EM pair supplies the actual
    # local derivative. No proxy value is used for target acceptance.
    cg=90.17920078631728+(151.7868750439733-90.17920078631728)*(coupling['d']-2)/26
    return np.array([effective_c,1/cg,coupling['d_t1']])


def suggest():
    base=ROOT/'output_data/taper_exploration'/CAMPAIGN
    features=np.array([describe(name) for name in CASES])
    center=features.mean(axis=0);scale=np.ptp(features,axis=0)
    design=np.c_[np.ones(len(features)),(features-center)/scale]
    if np.linalg.matrix_rank(design)!=4:raise ValueError('Independent completed directions are required')
    data=[np.load(base/name/'response/comparison.npz') for name in CASES]
    f=data[0]['f']
    if any(not np.array_equal(f,q['f']) for q in data):raise ValueError('Source grids differ')
    observed=np.array([q['observed'] for q in data])
    coefficients=np.linalg.solve(design,1/observed)
    cells=[TAPER_2,TAPER_1,UNIT_CELL,UNIT_CELL]
    target=transmission(f,[v.L0_nH for v in cells],[v.C0_fF for v in cells],
                        [v.Cg_input_fF for v in cells],0)
    # Stay on the measured 200-um-pad body: no new L simulation is required.
    net_c=features[-1,0]
    def row(d,junction):
        cg=90.17920078631728+(151.7868750439733-90.17920078631728)*(d-2)/26
        return np.r_[1,(np.array([net_c,1/cg,junction])-center)/scale]
    ix=np.arange(0,len(f),5)
    def prediction(d,junction):return 1/(row(d,junction)@coefficients)
    def residual(x):
        curve=prediction(x[0],x[1])[ix]*np.exp(-2j*np.pi*f[ix]*x[2]*1e-12)
        delta=curve-target[ix]
        weights=row(x[0],x[1])@np.linalg.inv(design)
        penalty=np.maximum(abs(weights)-1.5,0)*.1
        return np.r_[delta.real,delta.imag,penalty]
    fits=[least_squares(residual,[d,j,0],bounds=([8,44,-100],[42,100,100]),
                       x_scale=[20,40,20],max_nfev=200)
          for d in (12,20,28,36) for j in (52,68,84)]
    best=min(fits,key=lambda q:q.fun@q.fun)
    rows=[]
    for d in range(max(8,2*round(best.x[0]/2)-4),min(42,2*round(best.x[0]/2)+4)+1,2):
        for j in range(max(44,2*round(best.x[1]/2)-4),min(100,2*round(best.x[1]/2)+4)+1,2):
            curve=prediction(d,j)
            def phase_residual(x):
                delta=curve*np.exp(-2j*np.pi*f*x[0]*1e-12)-target
                return np.r_[delta.real,delta.imag]
            phase=least_squares(phase_residual,[best.x[2]],bounds=([-100],[100]),max_nfev=100)
            rows.append(dict(d_input=d,d_t1=j,target_RMS=float(np.sqrt(np.mean(phase.fun**2)*2)),
                       phase_ps=float(phase.x[0]),weights=(row(d,j)@np.linalg.inv(design)).tolist(),
                       maximum_magnitude=float(abs(curve).max())))
    rows.sort(key=lambda q:q['target_RMS'])
    proposal=rows[0]
    report=dict(scope='Unvalidated local same-topology surrogate; fresh EM is mandatory',
                cases=CASES,features=features.tolist(),condition=float(np.linalg.cond(design)),
                pad_length_fixed_um=200,ground_gap_fixed_um=54,ground_width_fixed_um=6,
                best=proposal,discrete_candidates=rows)
    (base/'thick_rail_surrogate_proposal.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez(base/'thick_rail_surrogate_prediction.npz',f=f,
             prediction=prediction(proposal['d_input'],proposal['d_t1']),target=target)
    print(json.dumps({k:v for k,v in report.items() if k!='discrete_candidates'},indent=2))
    return report


if __name__=='__main__':suggest()
