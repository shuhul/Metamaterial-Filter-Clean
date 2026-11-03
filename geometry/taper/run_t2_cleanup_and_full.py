"""Run the authorized cleanup comparisons, then one full mirrored filter ABS run."""
from pathlib import Path
from datetime import datetime,timezone
import json,os,time,shutil,traceback
import numpy as np
from geometry.taper import t2_cleanup as w

STATUS=w.OUT/'workflow_status.json'

def status(phase,**extra):
    r=w.read(STATUS) if STATUS.exists() else {}
    r.update(phase=phase,updated_utc=datetime.now(timezone.utc).isoformat(),launcher_pid=os.getpid(),**extra)
    w.save(STATUS,r);print('CLEANUP_STATUS',json.dumps(r),flush=True)

def archive_choice(choice):
    if choice['case']=='baseline_target5':return choice
    src=w.NATIVE/choice['case']
    dst=w.ROOT/'sonnet_files/taper_exploration/final_designs'/('T3_T2_cleanup_'+choice['case']+'_final')
    if not dst.exists():
        assert src.resolve().is_relative_to(w.NATIVE.resolve())
        assert dst.resolve().is_relative_to((w.ROOT/'sonnet_files/taper_exploration/final_designs').resolve())
        manifest=[dict(original=str(p),archive=str(dst/p.relative_to(src)),sha256=w.t3.digest(p)) for p in src.rglob('*') if p.is_file() and 'sondata' not in p.parts]
        shutil.move(str(src),str(dst));w.save(dst/'archive_manifest.json',manifest)
        shutil.copytree(w.OUT/choice['case'],dst/'evidence',ignore=shutil.ignore_patterns('*.resources.json'))
    else:assert not src.exists(),'Ambiguous archive and trial both exist'
    chosen=dict(choice,project=str(dst/'response/generated_ncell_4.sonx'),geometry=str(dst/'response/T2_closeup.png'),archive=str(dst))
    for ext in ('png','pdf'):shutil.copy2(Path(choice['plot']).with_suffix('.'+ext),dst/('final_s21_fit.'+ext))
    chosen['plot']=str(dst/'final_s21_fit.png');w.save(dst/'selection.json',chosen)
    return chosen

def full_run(chosen):
    from geometry.taper import full_filter as full
    case='full_filter_10_regular_3_tapers_each_side_T2_cleanup_2026-11-01'
    full.NATIVE=w.ROOT/'sonnet_files/taper_exploration/full_filter'/case
    full.OUT=w.ROOT/'output_data/taper_exploration/full_filter'/case
    full.PROJECT_NAME=case+'.sonx'
    full.source_default=lambda:Path(chosen['project'])
    project,meta=full.prepare()
    selection=dict(w.SELECTION,project=chosen['project'],analysis_fit_path=chosen['fit_path'],cleanup_choice=chosen)
    w.save(full.OUT/'selected_taper3_source.json',selection)
    w.save(w.ROOT/'geometry/taper/full_filter_cleanup_selection.json',dict(project=str(project),output=str(full.OUT),selection=selection,meta=meta))
    w.configure(project,4.65,5.4,7500)
    full.preflight(project)
    status('full_filter_EM_running',full_project=str(project),full_output=str(full.OUT),full_sweep='Native adaptive 4.65-5.4 GHz; allow healthy automatic refinement to finish',full_completed=False)
    # Fresh physical geometry starts its own cache; no earlier invalid interpolant.
    marker=full.OUT/'run_status.json'
    if marker.exists() and w.read(marker).get('status')=='complete':
        r=w.read(marker);assert w.t3.digest(project)==r['project_sha256'] and w.t3.digest(project.with_suffix('.snp'))==r['snp_sha256']
    else:w.simulate(project,full.OUT)
    w.check_geometry(project,full.NATIVE/'geometry.gds',meta,geometry_validator=full.validate)
    analysis=full.analyze(project)
    status('complete',full_completed=True,full_analysis=str(full.OUT/'full_filter_analysis.json'),full_plot=str(full.OUT/'full_filter_s21_magnitude.png'))
    return analysis

def main():
    w.OUT.mkdir(parents=True,exist_ok=True)
    status('baseline_target_bounds',state='running',started_utc=datetime.now(timezone.utc).isoformat(),original_selection_unchanged=True)
    # Correct the baseline extraction policy without changing its metal or old archive.
    baseline_project=w.SOURCE
    base_out=w.OUT/'baseline_target5'
    fit=w.fit_response(baseline_project,base_out,w.ANCHOR['L_nH'])
    assert fit['optimizer_success'] and fit['complex_RMS']<=.025,'Baseline cannot satisfy target bounds and fit cutoff'
    baseline=dict(case='baseline_target5',project=str(baseline_project),fit_path=str(base_out/'fit_all_caps_target5/fit.json'),fit_RMS=fit['complex_RMS'],fit_pass=True,L_nH=fit['L_nH'],plot=str(base_out/'final_s21_fit.png'))
    results=[baseline];w.save(w.OUT/'comparison.json',results)
    trim_l3=None
    for name,trim,center in [('trimmed_IDC',True,False),('centered_T2',False,True),('trimmed_IDC_centered_T2',True,True)]:
        try:
            ls=list(w.ANCHOR['L_nH'])
            if trim:
                if trim_l3 is None:
                    status('inductance_calibration',case=name,calibration_stage='T3')
                    trim_l3=w.calibrate('trimmed_IDC','T3',True,False)
                ls[0]=trim_l3['L0_nH']
            status('inductance_calibration',case=name,calibration_stage='T2')
            l2=w.calibrate(name,'T2',trim,center);ls[1]=l2['L0_nH']
            status('staged_EM_running',case=name,L_nH=ls)
            result=w.simulate_case(name,trim,center,ls)
            result['inductance_calibrations']=dict(T2=l2,T3=trim_l3 if trim else w.SELECTION['inductance'])
            result['L_pass']=bool(abs(ls[0]/2.77-1)<=.02 and abs(ls[1]/4.33-1)<=.02)
            result['accepted']=bool(result['fit_pass'] and result['L_pass'])
            results.append(result)
        except Exception as e:
            traceback.print_exc();results.append(dict(case=name,accepted=False,error=str(e)))
        w.save(w.OUT/'comparison.json',results)
        status('trial_complete',case=name,last_result=results[-1])
    accepted={r['case']:r for r in results if r.get('accepted',r['case']=='baseline_target5')}
    # Prefer both requested cleanups if they pass; otherwise retain a validated option.
    chosen=next(accepted[k] for k in ['trimmed_IDC_centered_T2','trimmed_IDC','centered_T2','baseline_target5'] if k in accepted)
    chosen=archive_choice(chosen);w.save(w.OUT/'selected_cleanup.json',chosen)
    status('cleanup_selected',selected=chosen)
    full_run(chosen)

if __name__=='__main__':
    try:main()
    except BaseException as e:
        traceback.print_exc();status('failed',state='failed',reason=str(e));raise
