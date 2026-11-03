"""Offline checks of target transmission and local fit sensitivity.

These diagnostics do not release fixed parameters or change the primary fit.
Peak widths are estimates on the exported frequency grid, not separate EM data.
"""
import json
from pathlib import Path
import numpy as np
from scipy.signal import find_peaks, peak_widths
from scipy.optimize import minimize_scalar
from geometry.taper.progressive_fitting import transmission
from geometry.taper.corrected_iteration import REGULAR_C, REGULAR_G


def peaks(f,s):
    power=abs(s)**2
    ix,_=find_peaks(power,prominence=max(power.max()*.01,1e-8))
    ix=np.array(sorted(ix,key=lambda i:power[i],reverse=True)[:2])
    ix.sort()
    widths=peak_widths(power,ix,rel_height=.5)[0] if len(ix) else []
    return [dict(frequency_GHz=float(f[i]/1e9),magnitude=float(abs(s[i])),
                 power_FWHM_MHz=float(w*np.median(np.diff(f))/1e6)) for i,w in zip(ix,widths)]


def best_target_delay(f,base,observed):
    """Only match fixture phase; target circuit capacitors remain fixed."""
    w=2*np.pi*f;cross=base*np.conj(observed)
    def cost(tau):return -float(np.real(np.sum(cross*np.exp(-1j*w*tau*1e-12))))
    grid=np.linspace(-100,100,401);i=int(np.argmin([cost(q) for q in grid]))
    result=minimize_scalar(cost,bounds=(grid[max(0,i-1)],grid[min(len(grid)-1,i+1)]),
                           method='bounded',options={'xatol':1e-8})
    return float(result.x),base*np.exp(-1j*w*result.x*1e-12)


def assess(folder):
    folder=Path(folder);r=json.loads((folder/'fit.json').read_text())
    data=np.load(folder/'response.npz');f=data['f'];observed=data['actual'][:,1,0]
    ls=r['L_nH'];cs=r['C_fF'];gs=r['G_fF'];free=r['free_capacitances']
    if len(ls)!=4 or len(free)!=5:raise ValueError('Expected the corrected full IDC fixture')
    target=transmission(f,ls,[100,100,REGULAR_C,REGULAR_C],[459,56.7,8.2865,REGULAR_G],r['delay_ps'])
    target_base=transmission(f,ls,[100,100,REGULAR_C,REGULAR_C],[459,56.7,8.2865,REGULAR_G],0)
    target_delay,target_complex=best_target_delay(f,target_base,observed)
    target_error=abs(target_complex-observed)
    x=np.array([(cs if k=='C' else gs)[i] for k,i in free]+[r['delay_ps']])
    # Fixed target scales make sensitivity comparisons independent of a bound
    # optimum (e.g. fitted C0=0). Differentiate on a much finer parameter step
    # than the frequency grid's narrowest peak displacement.
    scales=np.array([100,459,100,56.7,8.2865,20.])
    def evaluate(q):
        cp=list(cs);gp=list(gs)
        for (k,i),v in zip(free,q): (cp if k=='C' else gp)[i]=v
        s=transmission(f,ls,cp,gp,q[-1]);return np.r_[s.real,s.imag]
    columns=[]
    for i,scale in enumerate(scales):
        delta=np.zeros(len(x));delta[i]=scale*1e-6
        columns.append((evaluate(x+delta)-evaluate(x-delta))/2e-6)
    _,sv,v=np.linalg.svd(np.array(columns).T,full_matrices=False)
    result=dict(scope='Diagnostics only; fixed L and unit parameters preserved',
        source=r['provenance']['source'],
        sample_step_kHz=float(np.median(np.diff(f))/1e3),physical_peaks=peaks(f,observed),
        fitted_peaks=peaks(f,data['model'][:,1,0]),target_capacitance_peaks_at_measured_L=peaks(f,target),
        target_capacitance_magnitude_RMS_at_measured_L=float(np.sqrt(np.mean((abs(observed)-abs(target))**2))),
        scaled_singular_values=sv.tolist(),scaled_condition=float(sv[0]/sv[-1]),
        weakest_scaled_parameter_direction=dict(zip(['C0_T2','Cg_input_T2','C0_T1','Cg_T1','Cg_unit','delay'],v[-1].tolist())),
        sensitivity_is_not_a_confidence_interval=True)
    result.update(target_C_fF=[100,100,REGULAR_C,REGULAR_C],target_G_fF=[459,56.7,8.2865,REGULAR_G],
        fixed_L_nH=ls,target_only_delay_ps=target_delay,
        target_complex_RMS=float(np.sqrt(np.mean(target_error**2))),target_complex_max=float(target_error.max()),
        target_response_curve_pass=bool(np.sqrt(np.mean(target_error**2))<.01 and target_error.max()<.05),
        target_comparison_scope='Taper target capacitances held fixed; measured inductances and accepted regular constants held fixed. Only total matched fixture delay is fitted. Does not replace the primary all-free-capacitance fit or certify unique capacitances.')
    np.savez(folder/'target_response.npz',f=f,actual=observed,target=target_complex)
    (folder/'transmission_assessment.json').write_text(json.dumps(result,indent=2)+'\n')
    return result
