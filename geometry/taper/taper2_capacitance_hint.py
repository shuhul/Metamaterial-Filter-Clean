"""Offline geometry guidance from extracted capacitances and measured S21.

An affine local capacitance model proposes dimensions; the last measured curve
supplies a fixed residual correction. Predictions never accept a design. Every
candidate still needs native Sonnet, independent L, and the fixed-target score.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json,itertools
import numpy as np
from scipy.optimize import least_squares,minimize_scalar
from geometry.taper import taper2_plain as w


def propose(cases,*,fixed_body=None):
    base=w.ROOT/'output_data/taper_exploration'/w.CAMPAIGN
    native=w.ROOT/'sonnet_files/taper_exploration'/w.CAMPAIGN
    dims=[];caps=[];fits=[]
    for case in cases:
        s=json.loads((native/case/'response/geometry.json').read_text())['trial_settings']
        c=s['coupling']
        if c['n']!=70 or c['d_finger']!=4:raise ValueError('Use the same IDC family')
        dims.append([s['meander']['w'],s['pad']['l'],c['d'],c['d_t1'],c['d_regular'],c['t_regular_ground']])
        fit=json.loads((base/case/'response/extraction/fit.json').read_text());fits.append(fit)
        caps.append([fit['C_fF'][0],fit['C_fF'][1],*fit['G_fF'][:3]])
    dims=np.array(dims,float);caps=np.array(caps,float)
    dx=dims[1:]-dims[0];dc=caps[1:]-caps[0]
    if np.linalg.matrix_rank(dx)<6:raise ValueError('Need six measured geometry directions')
    slopes=np.linalg.lstsq(dx,dc,rcond=None)[0]
    data=np.load(base/cases[0]/'response/comparison.npz');f=data['f'];obs=data['observed'];target=data['target']
    ref=fits[0];ls=ref['L_nH'];ix=np.arange(0,len(f),5)
    def circuit(point,indices):
        c0,c1,g0,g1,gu=caps[0]+(point-dims[0])@slopes
        return w.common.transmission(f[indices],ls,[c0,c1,*ref['C_fF'][2:]],
                                     [g0,g1,gu,*ref['G_fF'][3:]],0)
    old=circuit(dims[0],ix)
    lower=np.array([148,20,12,4,32,2],float);upper=np.array([152,586,56,120,120,10],float)
    if fixed_body is not None:
        # Freeze w and pad length exactly in the optimization.
        active=np.arange(2,6);anchor=np.array([*fixed_body,*dims[0,2:]])
    else:active=np.arange(1,6);anchor=dims[0].copy()
    def expand(q):
        p=anchor.copy();p[active]=q[:-1];return p
    def residual(q):
        p=expand(q);pred=obs[ix]*circuit(p,ix)/old*np.exp(-2j*np.pi*f[ix]*q[-1]*1e-12)
        delta=pred-target[ix]
        positivity=np.maximum(1-(caps[0]+(p-dims[0])@slopes),0)
        return np.r_[delta.real,delta.imag,positivity]
    best=None
    starts=[dims[0,active],*(d[active] for d in dims[1:])]
    for seed in starts:
        q=least_squares(residual,np.r_[np.clip(seed,lower[active],upper[active]),0],
            bounds=(np.r_[lower[active],-100],np.r_[upper[active],100]),max_nfev=150)
        if best is None or best.fun@best.fun>q.fun@q.fun:best=q
    continuous=expand(best.x);rows=[];full=np.arange(len(f));oldfull=circuit(dims[0],full)
    floor=np.floor(continuous/2)*2
    for shift in itertools.product((0,1),repeat=len(active)):
        p=continuous.copy();p[active]=floor[active]+2*np.array(shift)
        if (p[4]-p[5])%4:continue
        curve=obs*circuit(p,full)/oldfull
        objective=lambda t:float(np.mean(abs(curve*np.exp(-2j*np.pi*f*t*1e-12)-target)**2))
        phase=minimize_scalar(objective,bounds=(-100,100),method='bounded')
        rows.append(dict(dimensions_um=p.tolist(),predicted_RMS=float(np.sqrt(phase.fun)),
                         predicted_capacitances_fF=(caps[0]+(p-dims[0])@slopes).tolist()))
    report=dict(scope='Extracted-capacitance heuristic only; fresh native EM required',sources=cases,
        dimension_names=['meander_w','pad_l','input_gap','T2_T1_gap','T1_unit_gap','ground_strip_w'],
        capacitance_names=['C0_T2','C0_T1','Cg_T2','Cg_T1','Cg_unit'],
        capacitance_slopes_fF_per_um=slopes.tolist(),continuous_dimensions_um=continuous.tolist(),
        candidates=sorted(rows,key=lambda q:q['predicted_RMS'])[:8])
    label='capacitance_hint_fixed_body' if fixed_body else 'capacitance_hint'
    (base/(label+'.json')).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2));return report
