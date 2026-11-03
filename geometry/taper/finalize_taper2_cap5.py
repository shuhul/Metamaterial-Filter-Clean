"""Explicitly calibrate and archive a verified5%-capacitance T2 design.

No theoretical-S21 acceptance gate. Fit-to-Sonnet quality, explicit L, geometry,
source identity, materials and native completion still have to pass.
"""
from pathlib import Path
from hashlib import sha256
import json,shutil
import xml.etree.ElementTree as ET
import numpy as np
import gdstk
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from geometry.taper.taper2_cap_target_iteration import ROOT,OUT,OLD_OUT,NATIVE,workflow,fit_case,OLDER_REFERENCE
from geometry.specifications import TAPER_2,TAPER_1,UNIT_CELL
from geometry.prepare_sonnet import tag


def calibrate(case):
    native=NATIVE/case
    meta=json.loads((native/'response/geometry.json').read_text());s=workflow.from_dict(meta['trial_settings'])
    measured=workflow.measure(case,s)
    result=fit_case(case,l0=measured['L0_nH'],provisional=False)
    result['explicit_inductance']=measured
    (OUT/case/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


def finalize(case):
    native=NATIVE/case;out=OUT/case
    result=json.loads((out/'result.json').read_text());measured=result['explicit_inductance']
    if not result['capacitance_pass'] or not result['fit_quality_pass'] or result['provisional_L']:
        raise ValueError('New capacitances and fit quality must pass with calibrated L')
    if measured['Python_frequency_approximation'] is not False or abs(measured['target_error_percent'])>2:
        raise ValueError('Explicit T2 L must be within2%')
    meta=json.loads((native/'response/geometry.json').read_text());s=workflow.from_dict(meta['trial_settings'])
    workflow.check_geometry(native/'response/generated_ncell_3.sonx',native/'response/geometry.gds',meta,
                            geometry_validator=workflow.validate)
    workflow.common.require_locked(gdstk.read_gds(str(native/'response/geometry.gds')).top_level()[0].get_polygons(),meta)
    status=json.loads((out/'run_status.json').read_text())
    for suffix,key in [('.sonx','project_sha256'),('.snp','snp_sha256')]:
        if sha256((native/'response'/('generated_ncell_3'+suffix)).read_bytes()).hexdigest()!=status[key]:
            raise ValueError('Completed EM data changed')
    log=(out/'em.log').read_text()
    if not all(x in log for x in ['Em errors detected: 0','Em warnings detected: 0','Em simulation completed']):
        raise ValueError('Native completion/warnings check failed')
    if list(ET.parse(native/'response/generated_ncell_3.sonx').iter(tag('Component'))):raise ValueError('Ideal component in production')
    cal_meta=json.loads((native/'Lcal/geometry.json').read_text())
    if cal_meta['unit']!=meta['unit']:raise ValueError('Calibration and production body differ')
    workflow.check_geometry(native/'Lcal_explicit/generated_ncell_1.sonx',native/'Lcal/geometry.gds',cal_meta,
                            geometry_validator=workflow.validate)
    calibration=OLD_OUT/case/'Lcal_explicit'
    proof=json.loads((calibration/'completed.json').read_text())
    for suffix,key in [('.sonx','project_sha256'),('.snp','snp_sha256')]:
        if sha256((native/'Lcal_explicit'/('generated_ncell_1'+suffix)).read_bytes()).hexdigest()!=proof[key]:
            raise ValueError('Explicit calibration data changed')
    if sha256((calibration/'results.json').read_bytes()).hexdigest()!=proof['result_sha256']:
        raise ValueError('Explicit L result changed')
    fitted=json.loads((out/'fit/fit.json').read_text());data=np.load(out/'fit/response.npz')
    if fitted['provenance']['sha256']!=status['snp_sha256']:
        raise ValueError('Fit source differs from completed response')
    if result['fixed_older_reference_sha256']!=sha256(OLDER_REFERENCE.read_bytes()).hexdigest():
        raise ValueError('Older reference changed')
    if not np.isclose(fitted['L_nH'][0],measured['L0_nH'],rtol=1e-12):raise ValueError('Final fit L differs from calibration')
    ts=[TAPER_2,TAPER_1,UNIT_CELL,UNIT_CELL];rows=[]
    labels=[['T2 L0','T2 C0','Input-to-T2 Cg'],['T1 L0','T1 C0','T2-to-T1 Cg'],
            ['Regular L0','Regular C0','T1-to-regular Cg'],['Output-end L0','Output-end C0','Regular-to-output Cg']]
    for i,t in enumerate(ts):
        for j,(value,target,unit) in enumerate(zip([fitted['L_nH'][i],fitted['C_fF'][i],fitted['G_fF'][i]],
                [t.L0_nH,t.C0_fF,t.Cg_input_fF],['nH','fF','fF'])):
            if j>0 and abs(100*(value/target-1))>5:
                raise ValueError(f'Capacitance outside5%: {labels[i][j]}')
            free=(j==1 and i==0) or (j==2 and i in (0,1))
            rows.append([labels[i][j],f'{value:.3f} {unit}',f'{target:g} {unit}',f'{100*(value/target-1):+.2f}%',
                         'Fitted' if free else 'Fixed'])
    fig,(ax,tab)=plt.subplots(1,2,figsize=(20,9),gridspec_kw={'width_ratios':[1.2,1]})
    ax.plot(data['f']/1e9,abs(data['actual'][:,1,0]),lw=2.8,color='#2768aa',label='Sonnet EM (2 µm mesh)')
    ax.plot(data['f']/1e9,abs(data['model'][:,1,0]),'k--',lw=2.2,label='Three-capacitance complex fit')
    ax.text(.025,.97,f'All three new capacitances within 5%\nFit-to-Sonnet RMS: {fitted["complex_RMS"]:.6f}\nExplicit L T2: {measured["L0_nH"]:.5f} nH',
        transform=ax.transAxes,va='top',fontsize=14,bbox=dict(boxstyle='round,pad=.5',fc='white',ec='gray'))
    ax.set(xlim=(4.65,5.4),ylim=(0,max(abs(data['actual'][:,1,0]).max(),abs(data['model'][:,1,0]).max())*1.5),xlabel='Frequency (GHz)',ylabel=r'$|S_{21}|$')
    ax.xaxis.label.set_size(20);ax.yaxis.label.set_size(23);ax.tick_params(labelsize=15);ax.grid(alpha=.2);ax.legend(loc='upper right',fontsize=11)
    tab.axis('off');tab.set_title('Parameter Comparison',fontsize=21)
    table=tab.table(cellText=rows,colLabels=['Parameter','Used value','Target','Difference','Status'],
        colWidths=[.30,.22,.20,.16,.12],cellLoc='left',colLoc='left',bbox=[0,.18,1,.80])
    table.auto_set_font_size(False);table.set_fontsize(11.5)
    for (r,c),cell in table.get_celld().items():
        cell.set_edgecolor('#d6dce1');cell.set_facecolor('#e0e8ef' if r==0 else '#f4f7fa' if ((r-1)//3)%2==0 else 'white')
        if r==0 or rows[r-1][-1]=='Fitted':cell.set_text_props(weight='bold')
    tab.text(0,.12,'Older capacitances use the adopted 1% reference and stay fixed.\nAll inductances fixed; one fixture delay fitted.\nCapacitance targets are the current acceptance objective.',transform=tab.transAxes,fontsize=11.5,va='top',linespacing=1.5)
    fig.suptitle('Second Taper Cell: Transmission and Parameter Comparison',fontsize=24,y=.98)
    fig.tight_layout(rect=[0,0,1,.95],w_pad=2)
    for ext in ('png','pdf'):fig.savefig(out/f'final_s21_fit.{ext}',dpi=180)
    plt.close(fig)
    destination=ROOT/'sonnet_files/taper_exploration/final_designs'/f'new_T2_cap5_{case}_final'
    base=(ROOT/'sonnet_files/taper_exploration').resolve()
    if not native.resolve().is_relative_to(base) or not destination.resolve().is_relative_to(base) or destination.exists():
        raise ValueError('Unsafe or existing archive destination')
    destination.parent.mkdir(exist_ok=True);native.rename(destination)
    shutil.copytree(out,destination/'evidence');shutil.copytree(calibration,destination/'Lcal_explicit/evidence')
    for name in ['final_s21_fit.png','final_s21_fit.pdf']:shutil.copy2(out/name,destination/name)
    shutil.copy2(OLD_OUT/case/'response/geometry.png',destination/'geometry.png')
    shutil.copy2(OLDER_REFERENCE,destination/'fixed_older_reference.json')
    top,_=workflow.build(s);body=next(r.cell for r in top.references if r.cell.name=='ASYMMETRIC_T2')
    lib=gdstk.Library(unit=1e-6,precision=1e-9);unit=lib.new_cell('TOP');unit.add(gdstk.Reference(body));lib.add(body,*body.dependencies(True));lib.write_gds(str(destination/'t2_unit.gds'))
    selection=dict(case=case,acceptance_basis='Three new fitted Cs within5%, fit-to-Sonnet RMS<=0.025, explicit L within2%; theoretical target S21 is not an acceptance gate.',
        T2_capacitance_goal_met=True,all_circuit_capacitances_within_5_percent=True,
        settings=meta['trial_settings'],result=result,
        project=str(destination/'response/generated_ncell_3.sonx'),plot=str(destination/'final_s21_fit.png'),
        geometry_plot=str(destination/'geometry.png'),unit_gds=str(destination/'t2_unit.gds'),
        fixed_older_reference=str(OLDER_REFERENCE),native_material_mesh_geometry_audit_pass=True,
        T1_regular_output_component_shapes_unchanged=True)
    for p in [destination/'selection.json',ROOT/'geometry/taper/selected_taper2_cap5.json']:
        p.write_text(json.dumps(selection,indent=2)+'\n')
    (destination/'README.txt').write_text('Accepted T2 under the5%-capacitance objective.\n\n'
        'Design: response/generated_ncell_3.sonx\nPlot: final_s21_fit.png\nBody only: t2_unit.gds\n'
        'The production project contains no ideal components. Only the three new taper capacitances and fixture delay were fitted.\n'
        'The adopted older circuit reference, original T1/regular/output shapes, materials and2um mesh are preserved.\n'
        'Explicit Sonnet capacitor sweeps establish L; evidence is in Lcal_explicit/evidence.\n'
        'This acceptance is not a claim that the full theoretical filter S21 has been achieved.\n')
    progress=ROOT/'geometry/taper/taper2_progress.json';d=json.loads(progress.read_text())
    d.update(state='T2_capacitance_goal_met',active_case=None,T2_capacitance_goal_met=True,
             T2_accepted=True,acceptance_basis=selection['acceptance_basis'],
             capacitance_selection=str(destination/'selection.json'),capacitance_project=selection['project'],capacitance_plot=selection['plot'])
    progress.write_text(json.dumps(d,indent=2)+'\n')
    print(json.dumps(selection,indent=2),flush=True);return selection
