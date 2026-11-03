"""Restore the accepted T2 input IDC and hold its extracted Cg fixed.

Only T3 C0, T3 input Cg, and the existing fixture delay are fitted.
No change to older capacitances or independently calibrated inductances.
"""
import json
import numpy as np
import gdstk
from geometry.taper import taper3 as work
from geometry.taper.sonnet import configure,simulate
from geometry.taper.progressive_fitting import fit


def run(case,p,l0):
    actual=(p.d_t2,p.n_t2,p.t_finger,p.d_finger,p.d_tip,p.t_input_rail,p.body.cap.t_right)
    if actual!=(22,70,4,4,2,8,8):
        raise ValueError('Restore the accepted T2 IDC: gap22,70fingers,trace4/gap4,tip2,rails8')
    project,out,meta=work.prepare(case,p,stage='response')
    verify_junction(case)
    status=out/'run_status.json'
    if status.exists() and json.loads(status.read_text()).get('status')=='complete':
        saved=json.loads(status.read_text())
        if work.digest(project)!=saved['project_sha256'] or work.digest(project.with_suffix('.snp'))!=saved['snp_sha256']:
            raise ValueError('Completed response changed')
    else:
        configure(project,4.65,5.4,7500)
        simulate(project,out)
    result=fit_saved(case,l0)
    progress_file=work.ROOT/'geometry/taper/taper3_progress.json'
    progress=json.loads(progress_file.read_text())
    progress.update(active_case=case,accepted=False,state='fixed_junction_trial_completed',latest_fixed_junction_result=result)
    progress_file.write_text(json.dumps(progress,indent=2)+'\n')
    return result


def verify_junction(case):
    stage=work.NATIVE/case/'response'
    meta=json.loads((stage/'geometry.json').read_text())
    old_left=work.SOURCE_META['cell_lefts'][1]
    shift=meta['downstream_shift_um']
    roi=gdstk.rectangle((old_left-30,-598),(old_left+8,-8))
    old=gdstk.read_gds(str(work.SOURCE_GDS)).top_level()[0].get_polygons()
    new=gdstk.read_gds(str(stage/'geometry.gds')).top_level()[0].get_polygons()
    a=gdstk.boolean(old,[roi],'and',precision=.001)
    b=gdstk.boolean([p.copy().translate(-shift,0) for p in new],[roi],'and',precision=.001)
    diff=gdstk.boolean(a,b,'xor',precision=.001)
    proof=dict(source_project=str(work.SOURCE),source_gds_sha256=work.digest(work.SOURCE_GDS),
        new_gds_sha256=work.digest(stage/'geometry.gds'),downstream_translation_um=shift,
        selected_junction_roi_um=[[old_left-30,-598],[old_left+8,-8]],
        rail_and_finger_XOR_area_um2=sum(p.area() for p in diff),junction_metal_identical=not diff)
    (work.OUT/case/'selected_IDC_metal_check.json').write_text(json.dumps(proof,indent=2)+'\n')
    if diff:raise ValueError('Selected T2 junction rails/fingers changed')
    return proof


def fit_saved(case,l0):
    work.require_reference()
    verify_junction(case)
    calibrated=json.loads((work.OUT/case/'inductance.json').read_text())['L0_nH']
    if abs(l0-calibrated)>1e-10:raise ValueError('Use independently calibrated L, never a fitted L')
    project=work.NATIVE/case/'response/generated_ncell_4.sonx'
    out=work.OUT/case/'response/fit_fixed_T2junction'
    refs=work.REFERENCE_FIT
    ls=[l0]+refs['L_nH'];cs=[251.]+refs['C_fF'];gs=[261.]+refs['G_fF']
    result=fit(project.with_suffix('.snp'),out,ls,cs,gs,[('C',0),('G',0)],
        capacitance_seeds=[[100,100],[251,261],[500,261],[800,500],[251,800]],
        title='Third Taper: Accepted T2 Junction Held Fixed')
    if result['G_fF'][1:]!=refs['G_fF'] or result['C_fF'][1:]!=refs['C_fF'] or result['L_nH']!=ls:
        raise ValueError('A fixed parameter changed')
    values=[result['C_fF'][0],result['G_fF'][0]]
    errors=(100*(np.array(values)/[251.,261.]-1)).tolist()
    summary=dict(case=case,fit_folder='fit_fixed_T2junction',values_fF=values,
        targets_fF=[251.,261.],errors_percent=errors,fit_RMS=result['complex_RMS'],
        components_pass=max(abs(v) for v in errors)<=5,
        fit_pass=result['curve_pass'] and result['interior_free_parameters'],
        measured_L_nH=l0,fixed_T3_T2_Cg_fF=gs[1],all_older_parameters_fixed=True,
        all_inductances_fixed=True,reference_file=str(work.REFERENCE_FILE),
        reference_sha256=work.REFERENCE_HASH,project=str(project))
    (out/'result.json').write_text(json.dumps(summary,indent=2)+'\n')
    rows=[json.loads(q.read_text()) for q in work.OUT.glob('*/response/fit_fixed_T2junction/result.json')]
    rows.sort(key=lambda q:(not q['fit_pass'],max(abs(v) for v in q['errors_percent'])))
    (work.OUT/'fixed_junction_ledger.json').write_text(json.dumps(rows,indent=2)+'\n')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    data=np.load(out/'response.npz')
    fig,ax=plt.subplots(figsize=(14,8))
    ax.plot(data['f']/1e9,abs(data['actual'][:,1,0]),lw=2.7,label='Sonnet EM (2 um mesh)')
    ax.plot(data['f']/1e9,abs(data['model'][:,1,0]),'k--',lw=2.2,label='Complex circuit fit')
    text=(f'T3 L = {l0:.4f} nH (target2.77; {100*(l0/2.77-1):+.2f}%; fixed)\n'
          f'T3 C0 = {values[0]:.2f} fF (target251; {errors[0]:+.2f}%)\n'
          f'T3 input Cg = {values[1]:.2f} fF (target261; {errors[1]:+.2f}%)\n'
          f'T3-T2 Cg = {gs[1]:.4f} fF (target153; {100*(gs[1]/153-1):+.2f}%; fixed)\n'
          f'Complex S21 RMS = {result["complex_RMS"]:.5f}\n'
          'All older Cs and all L fixed; fixture delay fitted')
    ax.text(.025,.97,text,transform=ax.transAxes,va='top',fontsize=14,
        bbox=dict(boxstyle='round,pad=.5',fc='white',ec='gray'))
    ax.set(xlabel='Frequency (GHz)',ylabel=r'$|S_{21}|$',xlim=(4.65,5.4),
        ylim=(0,max(abs(data['actual'][:,1,0]).max(),abs(data['model'][:,1,0]).max())*1.7))
    ax.xaxis.label.set_size(20);ax.yaxis.label.set_size(24);ax.tick_params(labelsize=15)
    ax.grid(alpha=.2);ax.legend(loc='upper right',fontsize=13)
    ax.set_title('Third Taper Cell: Transmission with T3-T2 Coupling Fixed',fontsize=21)
    fig.tight_layout()
    for ext in ('png','pdf'):fig.savefig(out/f's21_fit_parameters.{ext}',dpi=180)
    plt.close(fig)
    from geometry.taper.taper3_fit_sensitivity import analyze
    analyze(case,fit_folder='fit_fixed_T2junction')
    print('FIXED JUNCTION RESULT',json.dumps(summary),flush=True)
    return summary
