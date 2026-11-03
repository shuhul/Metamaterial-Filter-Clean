"""Response-level diagnostics, independent of nonunique lumped-C extraction."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json
import numpy as np
from scipy.optimize import least_squares, lsq_linear
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from circuit_common import read_snp
from geometry.taper.target_iteration import target_response
from geometry.taper.run_bridged_stage import digest


def peak_features(f,s):
    power=abs(s)**2;indices,info=find_peaks(abs(s),prominence=.02)
    result=[]
    for i in indices:
        level=power[i]/2;left=i;right=i
        while left>0 and power[left]>level:left-=1
        while right<len(f)-1 and power[right]>level:right+=1
        width=None
        if power[left]<=level and power[right]<=level:
            fl=np.interp(level,power[left:left+2],f[left:left+2])
            fr=np.interp(level,power[right-1:right+1][::-1],f[right-1:right+1][::-1])
            width=float((fr-fl)/1e6)
        result.append(dict(index=int(i),frequency_GHz=float(f[i]/1e9),
            magnitude=float(abs(s[i])),power_FWHM_MHz=width))
    return result


def local_pole(f,channels,peak_GHz,width_MHz,half_window_MHz):
    """Variable projection: shared pole; separate complex residues/backgrounds.

    exp(+jwt) convention: denominator (f-f0-i*gamma), gamma>0.
    The reported linewidth is 2*gamma, not a separately extracted external Q.
    """
    mask=abs(f/1e9-peak_GHz)<=half_window_MHz/1000
    x=(f[mask]-peak_GHz*1e9)/1e6;y=channels[mask]
    def design(v):
        return np.column_stack([np.ones(len(x)),x/half_window_MHz,1/(x-v[0]-1j*v[1])])
    def residual(v):
        a=design(v);co=np.linalg.lstsq(a,y,rcond=None)[0];d=a@co-y
        return np.r_[d.real.ravel(),d.imag.ravel()]
    q=least_squares(residual,[0,width_MHz/2],
        bounds=([-min(5,half_window_MHz/2),.01],[min(5,half_window_MHz/2),half_window_MHz]),
        max_nfev=120,ftol=1e-11,xtol=1e-11,gtol=1e-11)
    center=peak_GHz+q.x[0]/1000;linewidth=2*q.x[1]
    co=np.linalg.lstsq(design(q.x),y,rcond=None)[0]
    return dict(frequency_GHz=float(center),loaded_linewidth_MHz=float(linewidth),
        loaded_Q=float(center*1000/linewidth),window_GHz=[float(f[mask][0]/1e9),float(f[mask][-1]/1e9)],
        local_complex_RMS=float(np.sqrt(np.mean(q.fun**2)*2)),optimizer_success=bool(q.success),
        active_bounds=q.active_mask.tolist(),residue_magnitudes=np.abs(co[-1]).tolist(),
        channels=['S11','S21','S22'],
        scope='Local common-pole fit; does not uniquely separate internal/external decay or lumped capacitors')


def analyze(source,folder):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    net=read_snp(source);f=net.f;s=net.s;goal=target_response(f)
    peaks=peak_features(f,s[:,1,0]);wanted=peak_features(f,goal)
    channels=np.column_stack([s[:,0,0],s[:,1,0],s[:,1,1]])
    poles=[]
    for p in peaks:
        separation=min([abs(p['frequency_GHz']-o['frequency_GHz'])*1000 for o in peaks if o!=p] or [250])
        window=min(.35*separation,max(12,6*(p['power_FWHM_MHz'] or 2)))
        if p['power_FWHM_MHz']:
            poles.append(local_pole(f,channels,p['frequency_GHz'],p['power_FWHM_MHz'],window))
    deficit=1-abs(s[:,0,0])**2-abs(s[:,1,0])**2
    result=dict(source=str(source),source_sha256=digest(source),peaks=peaks,target_peaks=wanted,
        local_poles=poles,reciprocity_max_error=float(np.max(abs(s[:,1,0]-s[:,0,1]))),
        missing_twoport_power=dict(minimum=float(deficit.min()),maximum=float(deficit.max()),
            values_at_peaks=[float(deficit[p['index']]) for p in peaks],
            scope='Incident-port-1 power not returned through S11/S21; can include radiation, material loss and calibration/numerical effects'),
        reflection_magnitudes_at_peaks=[dict(S11=float(abs(s[p['index'],0,0])),S22=float(abs(s[p['index'],1,1]))) for p in peaks])
    (folder/'response_diagnostics.json').write_text(json.dumps(result,indent=2)+'\n')
    fig,axs=plt.subplots(2,2,figsize=(15,10),sharex=True)
    for ax,(i,j) in zip(axs.flat[:3],[(1,0),(0,0),(1,1)]):
        ax.plot(f/1e9,abs(s[:,i,j]),lw=2,color='#2768aa',label='Sonnet')
        if (i,j)==(1,0):ax.plot(f/1e9,abs(goal),'--',lw=2,color='#db6d28',label='Target circuit')
        ax.set_ylabel(f'|S{i+1}{j+1}|',fontsize=18);ax.grid(alpha=.2);ax.legend()
    axs[1,1].plot(f/1e9,deficit,lw=2,color='#8561a8')
    axs[1,1].set_ylabel('Missing two-port power fraction',fontsize=14)
    axs[1,1].axhline(0,color='gray',lw=.7);axs[1,1].grid(alpha=.2)
    for ax in axs[1]:ax.set_xlabel('Frequency (GHz)',fontsize=16)
    for ax in axs.flat:ax.tick_params(labelsize=12)
    fig.suptitle('Plain-Gap Layout: Transmission, Reflection, and Power Balance',fontsize=21)
    fig.tight_layout();fig.savefig(folder/'s_parameter_diagnostics.png',dpi=170);plt.close(fig)
    print('RESPONSE DIAGNOSTICS',str(folder),json.dumps(dict(peaks=peaks,poles=poles)),flush=True)
    return result


def summarize_study():
    from geometry.taper.plain_gap_study import OUT,BASE_OUT,BASE_NATIVE,configurations
    OUT.mkdir(parents=True,exist_ok=True)
    base_path=OUT/'baseline/response_diagnostics.json'
    if not base_path.exists():analyze(BASE_NATIVE.with_suffix('.snp'),base_path.parent)
    base=json.loads(base_path.read_text());rows=[]
    for name in configurations():
        path=OUT/name/'idc/response_diagnostics.json'
        if not path.exists():continue
        d=json.loads(path.read_text());intent=json.loads((path.parent/'intent.json').read_text())
        score=json.loads((path.parent/'target_comparison.json').read_text())
        row=dict(case=name,parameter=intent['changed_parameter'],target_complex_RMS=score['complex_RMS'],
            peaks=d['peaks'],local_poles=d['local_poles'],changes_vs_baseline=[])
        if len(d['peaks'])==len(base['peaks'])==2:
            for new,old in zip(d['peaks'],base['peaks']):
                row['changes_vs_baseline'].append(dict(frequency_shift_MHz=1000*(new['frequency_GHz']-old['frequency_GHz']),
                    width_change_MHz=new['power_FWHM_MHz']-old['power_FWHM_MHz'],
                    height_change=new['magnitude']-old['magnitude']))
        rows.append(row)
    report=dict(baseline=base,cases=rows,scope='One-sided finite parameter changes; response trends within tested neighborhood, not proof of global feasibility')
    (OUT/'sensitivity_results.json').write_text(json.dumps(report,indent=2)+'\n')
    lines=['# Approved plain-gap parameter study','',
        'Straight upper ground; plain 4 um T2-to-T1 gap; equal 584 um facing legs; original T1/regular metal locked.',
        '','All variants change exactly one parameter. Acceptance compares complex S21 to the fixed target circuit over 4.65–5.40 GHz. Only a matched fixture delay is free.',
        '','| Case | Changed parameter | Target complex RMS | Lower shift (MHz) | Upper shift (MHz) | Width changes (MHz) |',
        '|---|---|---:|---:|---:|---|']
    for r in rows:
        delta=r['changes_vs_baseline'];description=', '.join(f"{k}: {v['baseline']} -> {v['variant']} um" for k,v in r['parameter'].items())
        if len(delta)==2:
            lines.append(f"| {r['case']} | {description} | {r['target_complex_RMS']:.5f} | {delta[0]['frequency_shift_MHz']:+.2f} | {delta[1]['frequency_shift_MHz']:+.2f} | {delta[0]['width_change_MHz']:+.3f}, {delta[1]['width_change_MHz']:+.3f} |")
    lines+=['','Local common-pole fits use complex S11/S21/S22 to estimate loaded linewidths. They do not assign unique capacitances or distinguish intrinsic and external losses by themselves.',
            '','Missing two-port power is a diagnostic including possible radiation/material loss and calibration/numerical effects. It is not automatically identified as conductor loss.',
            '','Current diagnostics preserve baseline geometry/stack/ports exactly and use three direct frequencies at the two baseline peaks and their midpoint. Their S-parameters are checked against the saved adaptive sweep.']
    lines+=['','## Absolute response features','',
        '| Response | Lower peak (GHz) | Upper peak (GHz) | Lower FWHM (MHz) | Upper FWHM (MHz) |',
        '|---|---:|---:|---:|---:|']
    features=[('Fixed target',base['target_peaks']),('Approved baseline',base['peaks'])]
    features.extend((r['case'],r['peaks']) for r in rows)
    for name,p in features:
        if len(p)==2:
            lines.append(f"| {name} | {p[0]['frequency_GHz']:.4f} | {p[1]['frequency_GHz']:.4f} | {p[0]['power_FWHM_MHz']:.3f} | {p[1]['power_FWHM_MHz']:.3f} |")
    lines+=['','## Files and reproducibility','',
        '`parameter_comparison.png` overlays all completed cases and the fixed target, with peak zooms.',
        '`baseline_currents/current_charge_maps.png` contains validated whole-layout current/charge maps.',
        '`baseline_currents/vertical_current_phase.png` shows the Jy component phase.',
        '`baseline_currents/field_validation.json` records coordinate checks and neck-current phase samples.',
        '`baseline_currents/exports/REJECTED_EXPORTS.md` identifies export formats that failed validation.',
        '`theory_sensitivity.json` contains circuit-only +10% electrical-parameter diagnostics; it does not redefine the target.',
        'Each case has source hashes, intent, solver log/status, fixed-target comparison and S11/S21/S22 diagnostics.',
        '', 'Run `.venv\\Scripts\\python.exe -B -m geometry.taper.plain_gap_study` for resumable sequential EM tests. '
        'Run `-m geometry.taper.report_plain_gap_study` to refresh overlays/status after all cases finish.',
        '', 'These three small one-sided perturbations establish local trends. They cannot establish that the topology is globally incapable of meeting the target.']
    (OUT/'README.md').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    return report


if __name__=='__main__':
    from geometry.taper.plain_gap_study import BASE_NATIVE,OUT
    analyze(BASE_NATIVE.with_suffix('.snp'),OUT/'baseline')
    summarize_study()
