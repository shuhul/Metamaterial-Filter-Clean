"""Offline proposals for the approved thin-rail IDC family, never acceptance.

Local dimension sensitivities are approximate. Measured corrections interpolate
the inverse complex transmission, and every chosen proposal requires fresh EM.
No solver, target changes or geometry mutations occur in this module.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import json
import numpy as np
from scipy.optimize import minimize_scalar
from scipy.signal import find_peaks, peak_widths
from geometry.taper import taper2_plain as workflow

BASE='idc_rail8_dl256_pad148_w152_i40n70_j36_r70tg2_preview'


def peak_diagnostics(f, curve):
    """Peak position, height and half-prominence power width (diagnostics only)."""
    power=abs(curve)**2
    peaks,_=find_peaks(power,prominence=.005)
    widths=peak_widths(power,peaks,rel_height=.5)[0]
    return [dict(frequency_GHz=float(f[i]/1e9),magnitude=float(abs(curve[i])),
                 half_prominence_power_width_MHz=float(width*(f[1]-f[0])/1e6))
            for i,width in zip(peaks,widths)]


def ledger(active_case=None):
    """Save actual completed trials for the user-approved cleared IDC family."""
    root=workflow.ROOT/'output_data/taper_exploration'/workflow.CAMPAIGN
    native=workflow.ROOT/'sonnet_files/taper_exploration'/workflow.CAMPAIGN
    selection_path=workflow.ROOT/'geometry/taper/selected_taper2_idc.json'
    selection=json.loads(selection_path.read_text()) if selection_path.exists() else None
    rows=[]
    for path in root.glob('*/response/target_comparison.json'):
        name=path.parent.parent.name
        case_native=native/name
        if selection and selection['case']==name:
            case_native=Path(selection['project']).parent.parent
        metadata=case_native/'response/geometry.json'
        if not metadata.exists():continue
        meta=json.loads(metadata.read_text());clear=meta.get('IDC_clearance_um',{})
        if (meta['unit']['cap'].get('t_left')!=8 or
                clear.get('meander',0)<264 or clear.get('nearest_internal_ground_pad',0)<156):continue
        score=json.loads(path.read_text());fit=json.loads((path.parent/'extraction/fit.json').read_text())
        calibration=path.parent.parent/'Lcal_explicit/results.json'
        if not calibration.exists():calibration=path.parent.parent/'Lcal/results.json'
        cal=json.loads(calibration.read_text())
        curves=np.load(path.parent/'comparison.npz')
        rows.append(dict(case=name,target_RMS=score['complex_RMS'],fit_RMS=fit['complex_RMS'],
            target_pass=score['pass'],L_nH=cal['L0_nH'],L_error_percent=cal['target_error_percent'],
            L_method=cal.get('method','Sparse four-port operator; pending explicit cross-check'),
            C_fF=fit['C_fF'],G_fF=fit['G_fF'],settings=meta['trial_settings'],
            snp_sha256=score['source_sha256'],project=str(case_native/'response/generated_ncell_3.sonx'),
            plot=str(path.parent/'s21_fit.png'),
            EM_peaks=peak_diagnostics(curves['f'],curves['observed']),
            target_peaks=peak_diagnostics(curves['f'],curves['target'])))
        # A changing mutual C can trade against fitted ground C while their
        # sum stays nearly fixed. These are fitted nodal matrix diagonals,
        # not independent physical capacitance measurements.
        cc,gg=fit['C_fF'],fit['G_fF']
        rows[-1]['effective_nodal_capacitance_fF']={
            'T2':cc[0]+gg[0]+gg[1],'T1':cc[1]+gg[1]+gg[2]}
    rows.sort(key=lambda q:q['target_RMS'])
    report=dict(scope='Actual completed Sonnet trials only; predictions excluded',active_case=active_case,rows=rows)
    (root/'thin_rail_trials.json').write_text(json.dumps(report,indent=2)+'\n')
    lines=['# Cleared thin-rail IDC trials','','Same fixed nominal target,7501 samples over4.65-5.40GHz; one matched delay only.',
           '','','|Case|Target RMS|Fit RMS|L (nH)|','|---|---:|---:|---:|']
    for row in rows:lines.append(f"|{row['case']}|{row['target_RMS']:.6f}|{row['fit_RMS']:.6f}|{row['L_nH']:.7f}|")
    (root/'thin_rail_trials.md').write_text('\n'.join(lines)+'\n')
    if rows:
        progress=workflow.ROOT/'geometry/taper/taper2_progress.json'
        q=json.loads(progress.read_text());best=rows[0]
        q.update(best_case=best['case'],best_target_RMS=best['target_RMS'],best_project=best['project'],
            best_plot=best['plot'],best_isolated_L_nH=best['L_nH'],best_fit_RMS=best['fit_RMS'],
            best_current_geometry_case=best['case'],best_current_target_RMS=best['target_RMS'],active_case=active_case)
        q['new_sims_running']=active_case is not None
        if selection and selection['specification']==workflow.common.specification():
            q.update(T2_accepted=True,accepted_case=selection['case'],state='accepted')
        if active_case:
            active_meta=native/active_case/'response/geometry.json'
            active_cal=root/active_case/'Lcal_explicit/results.json'
            q['active_project']=str(native/active_case/'response/generated_ncell_3.sonx')
            if active_meta.exists():
                am=json.loads(active_meta.read_text());q['active_geometry']=am['trial_settings']
                q['IDC_clearances_um']=am.get('IDC_clearance_um',{})
            if active_cal.exists():
                ac=json.loads(active_cal.read_text())
                q.update(active_isolated_L_nH=ac['L0_nH'],active_isolated_L_error_percent=ac['target_error_percent'])
            else:
                q.pop('active_isolated_L_nH',None);q.pop('active_isolated_L_error_percent',None)
            q['active_inductance_calibration_required']=not active_cal.exists()
        else:q['active_project']=None
        progress.write_text(json.dumps(q,indent=2)+'\n')
    print(json.dumps(dict(completed=len(rows),best=rows[0]['case'] if rows else None,
        best_RMS=rows[0]['target_RMS'] if rows else None,active_case=active_case),indent=2))
    return report


def propose(observed_cases=(), *, source_case=BASE, L2_nH=None, C2_delta_fF=0,
            input_gaps=range(22,55,2),junction_gaps=range(12,65,2),regular_gaps=range(66,99,4)):
    root=workflow.ROOT/'output_data/taper_exploration'/workflow.CAMPAIGN
    native=workflow.ROOT/'sonnet_files/taper_exploration'/workflow.CAMPAIGN
    def load(name):
        folder=root/name/'response'
        meta=json.loads((native/name/'response/geometry.json').read_text())
        return np.load(folder/'comparison.npz'),json.loads((folder/'extraction/fit.json').read_text()),meta
    data,fit,meta=load(source_case);f=data['f'];actual=data['observed']
    old=workflow.common.transmission(f,fit['L_nH'],fit['C_fF'],fit['G_fF'],fit['delay_ps'])
    specs=[workflow.common.TAPER_2,workflow.common.TAPER_1,workflow.common.UNIT_CELL,workflow.common.UNIT_CELL]
    target=workflow.common.transmission(f,[s.L0_nH for s in specs],[s.C0_fF for s in specs],
                                      [s.Cg_input_fF for s in specs],0)
    def dimensions(m):
        c=m['trial_settings']['coupling']
        if c['t_regular_ground']!=2:raise ValueError('This local model keeps the 2 um strip')
        return np.array([c['d'],c['d_t1'],c['d_regular']],dtype=float)
    origin=dimensions(meta);scales=np.array([2,2,4])
    def parameters(point):
        din,junction,regular=point-origin
        cs=list(fit['C_fF']);gs=list(fit['G_fF'])
        cs[0]+=.05*din+.049*junction
        cs[1]+=-.10889127018*regular+.049*junction
        gs[0]+=3.5*din
        gs[1]*=(origin[1]+37.2)/(point[1]+37.2)
        gs[2]+=-.05712171406*regular
        return cs,gs
    def uncorrected(point,indices):
        cs,gs=parameters(point)
        m=workflow.common.transmission(f[indices],fit['L_nH'],cs,gs,fit['delay_ps'])
        return actual[indices]*m/old[indices]
    full=np.arange(len(f));steps=[];errors=[]
    for name in observed_cases:
        observed,_,measured_meta=load(name)
        if measured_meta['unit']!=meta['unit'] or not np.array_equal(observed['f'],f):
            raise ValueError('Correction requires the same body and frequency grid')
        point=dimensions(measured_meta);step=(point-origin)/scales
        steps.append(step);errors.append(1/observed['observed']-1/uncorrected(point,full))
    coefficients=np.linalg.pinv(np.array(steps))@np.array(errors) if steps else None
    def predict(point,indices):
        curve=uncorrected(point,indices)
        if coefficients is not None:
            curve=1/(1/curve+((point-origin)/scales)@coefficients[:,indices])
        if L2_nH is not None or C2_delta_fF:
            cs,gs=parameters(point)
            previous=workflow.common.transmission(f[indices],fit['L_nH'],cs,gs,fit['delay_ps'])
            ls=list(fit['L_nH']);cs=list(cs)
            if L2_nH is not None:ls[0]=L2_nH
            cs[0]+=C2_delta_fF
            updated=workflow.common.transmission(f[indices],ls,cs,gs,fit['delay_ps'])
            curve*=updated/previous
        return curve
    indices=np.arange(0,len(f),5);rows=[]
    for din in input_gaps:
        for junction in junction_gaps:
            for regular in regular_gaps:
                point=np.array([din,junction,regular]);curve=predict(point,indices)
                seed=float(np.clip(np.angle(np.vdot(target[indices],curve))/(2*np.pi*np.mean(f[indices]))*1e12,-90,90))
                objective=lambda tau:float(np.mean(abs(curve*np.exp(-2j*np.pi*f[indices]*tau*1e-12)-target[indices])**2))
                phase=minimize_scalar(objective,bounds=(max(-100,seed-20),min(100,seed+20)),method='bounded')
                rows.append(dict(input_gap=din,junction=junction,regular_gap=regular,
                                 predicted_RMS=float(np.sqrt(phase.fun)),delay_ps=float(phase.x)))
    rows.sort(key=lambda q:q['predicted_RMS'])
    for row in rows[:20]:
        curve=predict(np.array([row['input_gap'],row['junction'],row['regular_gap']]),full)
        objective=lambda tau:float(np.mean(abs(curve*np.exp(-2j*np.pi*f*tau*1e-12)-target)**2))
        phase=minimize_scalar(objective,bounds=(max(-100,row['delay_ps']-10),min(100,row['delay_ps']+10)),method='bounded')
        row.update(predicted_RMS=float(np.sqrt(phase.fun)),delay_ps=float(phase.x))
    report=dict(scope='Offline local model only; predictions do not accept a design',source=source_case,
                observed_corrections=list(observed_cases),body_unchanged=L2_nH is None and not C2_delta_fF,
                assumed_body_L_nH=L2_nH,assumed_body_C_delta_fF=C2_delta_fF,
                dimension_ranges_um=dict(input=list(input_gaps),junction=list(junction_gaps),regular=list(regular_gaps)),
                rows=sorted(rows,key=lambda q:q['predicted_RMS'])[:40])
    label='corrected_parameter_proposals_'+str(len(observed_cases))
    if not report['body_unchanged']:label+='_body_override'
    (root/source_case/'response'/(label+'.json')).write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report['rows'][:6],indent=2));return report


if __name__=='__main__':
    propose(['idc_rail8_dl256_pad148_w152_i34n70_j32_r78tg2'])
