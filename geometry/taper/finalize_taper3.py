"""Audit and archive an actually simulated T3 meeting the component policy."""
from pathlib import Path
import json
import re
import shutil
import xml.etree.ElementTree as ET
import numpy as np
import gdstk
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from geometry.taper import taper3 as work
from geometry.specifications import TAPER_3, TAPER_2, TAPER_1, UNIT_CELL


def anchored_policy(folder):
    match=re.fullmatch(r'fit_bounded_(1|1\.5|1\.6|2|3)pct(?:_j(fixed|1pct|1\.5pct|1\.6pct|2pct))?',folder)
    if not match:return None
    percent=float(match.group(1))
    if match.group(2) and match.group(2)!='fixed' and float(match.group(2).removesuffix('pct'))>percent:
        raise ValueError('Junction bound exceeds overall old-C allowance')
    return int(percent) if percent.is_integer() else percent


def plot_case(case, *, stage='response', fit_folder='fit'):
    native = work.NATIVE/case/stage; out = work.OUT/case/stage
    policy=anchored_policy(fit_folder)
    if fit_folder not in ('fit','fit_older_1percent','fit_fixed_T2junction') and policy is None:
        raise ValueError('Unknown fit variant')
    anchored_all_older = policy is not None
    bound_percent = policy if anchored_all_older else 1
    bounded = fit_folder == 'fit_older_1percent' or anchored_all_older
    fixed_junction = fit_folder == 'fit_fixed_T2junction'
    alternate = fit_folder != 'fit'
    if anchored_all_older:
        from geometry.taper.taper3_pad_tuning import plot_result
        result = json.loads((out/fit_folder/'result.json').read_text())
        plot_result(case,bound_percent,result,fit_folder=fit_folder)
        return out/fit_folder/'s21_fit_parameters.png'
    fitted = json.loads((out/fit_folder/'fit.json').read_text())
    result = json.loads((out/fit_folder/'result.json' if alternate else out/'result.json').read_text())
    data = np.load(out/fit_folder/'response.npz')
    electrical_policy = 'All L fixed; older Cs bounded to Â±1%' if bounded else 'All inductances and older Cs fixed'
    targets = [TAPER_3,TAPER_2,TAPER_1,UNIT_CELL,UNIT_CELL]
    names = ['T3','T2','T1','Regular','Output-end']
    rows=[]
    for i,t in enumerate(targets):
        for kind,value,target,units in [
            ('L0',fitted['L_nH'][i],t.L0_nH,'nH'),
            ('C0',fitted['C_fF'][i],t.C0_fF,'fF'),
            ('input Cg',fitted['G_fF'][i],t.Cg_input_fF,'fF')]:
            if value is None:
                rows.append([names[i]+' '+kind,'Direct feed','â€”','â€”','Fixed'])
                continue
            key = 'C' if kind=='C0' else 'G'
            free=(kind!='L0' and [key,i] in fitted['free_capacitances'] and
                  ((kind=='C0' and i==0) or (kind=='input Cg' and i in (0,1))))
            status='Fitted' if free else 'Â±1% fit' if bounded and kind!='L0' else 'Fixed'
            rows.append([names[i]+' '+kind,f'{value:.3f} {units}',f'{target:g} {units}',
                         f'{100*(value/target-1):+.2f}%',status])
    fig,(ax,tab)=plt.subplots(1,2,figsize=(21,10),gridspec_kw={'width_ratios':[1.2,1]})
    ax.plot(data['f']/1e9,abs(data['actual'][:,1,0]),color='#2768aa',lw=2.7,label='Sonnet EM (2 Âµm mesh)')
    ax.plot(data['f']/1e9,abs(data['model'][:,1,0]),'k--',lw=2.2,label='Complex circuit fit')
    ax.text(.025,.97,f'Fit-to-Sonnet RMS: {fitted["complex_RMS"]:.6f}\nNew capacitance tolerance: Â±5%\n{electrical_policy}',
        transform=ax.transAxes,va='top',fontsize=14,bbox=dict(boxstyle='round,pad=.5',fc='white',ec='gray'))
    ax.set(xlabel='Frequency (GHz)',ylabel=r'$|S_{21}|$',xlim=(4.65,5.4),
           ylim=(0,max(abs(data['actual'][:,1,0]).max(),abs(data['model'][:,1,0]).max())*1.4))
    ax.xaxis.label.set_size(20);ax.yaxis.label.set_size(24);ax.tick_params(labelsize=15)
    ax.grid(alpha=.2);ax.legend(loc='upper right',fontsize=12)
    if not (result['components_pass'] and result['fit_pass']):
        ax.text(.025,.74,'Exploratory run: component/fit checks not passed',transform=ax.transAxes,
                fontsize=12,color='#a32929',bbox=dict(fc='white',ec='#a32929',boxstyle='round'))
    tab.axis('off');tab.set_title('Parameter Comparison',fontsize=21)
    table=tab.table(cellText=rows,colLabels=['Parameter','Extracted','Target','Difference','Status'],
        colWidths=[.28,.23,.20,.17,.12],cellLoc='left',colLoc='left',bbox=[0,.15,1,.83])
    table.auto_set_font_size(False);table.set_fontsize(11.5)
    for (r,c),cell in table.get_celld().items():
        cell.set_edgecolor('#d6dce1');cell.set_facecolor('#e0e8ef' if r==0 else '#f4f7fa' if ((r-1)//3)%2==0 else 'white')
        if r==0 or rows[r-1][-1]=='Fitted':cell.set_text_props(weight='bold')
    footer = ('New Cs and delay fitted; older Cs may move Â±1% from selected T2.\nAll inductances remain fixed.\n' if bounded else
              'Only the new stage capacitances and total fixture delay are fitted.\nOlder values use the selected T2 model.\n')
    tab.text(0,.10,footer+'Nominal-target S21 error is not the current acceptance gate.',
             transform=tab.transAxes,fontsize=11.5,va='top',linespacing=1.5)
    fig.suptitle('Third Taper Cell: '+('Direct-Feed Verification' if stage=='direct' else 'Transmission and Parameter Comparison'),fontsize=24,y=.98)
    fig.tight_layout(rect=[0,0,1,.95],w_pad=2)
    destination = out/fit_folder if alternate else out
    for ext in ('png','pdf'):fig.savefig(destination/f's21_fit_parameters.{ext}',dpi=180)
    plt.close(fig)
    return destination/'s21_fit_parameters.png'


def finalize(case, *, fit_folder='fit'):
    policy=anchored_policy(fit_folder)
    if fit_folder not in ('fit','fit_older_1percent','fit_fixed_T2junction') and policy is None:
        raise ValueError('Unknown fit variant')
    anchored_all_older = policy is not None
    bound_percent = policy if anchored_all_older else 1
    bounded = fit_folder == 'fit_older_1percent' or anchored_all_older
    fixed_junction = fit_folder == 'fit_fixed_T2junction'
    alternate = fit_folder != 'fit'
    work.require_reference()
    native=work.NATIVE/case;out=work.OUT/case
    result_path=out/'response'/fit_folder/'result.json' if alternate else out/'response/result.json'
    result=json.loads(result_path.read_text())
    measured=json.loads((out/'inductance.json').read_text())
    if not result['components_pass'] or not result['fit_pass'] or abs(measured['target_error_percent'])>2:
        raise ValueError('Measured L, extracted Cs and fit quality must all pass')
    if measured.get('Python_frequency_approximation') is not False:
        raise ValueError('Explicit native capacitor calibration required')
    meta=json.loads((native/'response/geometry.json').read_text())
    work.check_geometry(native/'response/generated_ncell_4.sonx',native/'response/geometry.gds',meta,geometry_validator=work.validate)
    actual=[gdstk.Polygon(work.points(q)) for q in ET.parse(native/'response/generated_ncell_4.sonx').iter(work.tag('PlanarPolygon'))]
    work.require_downstream(actual,meta,native=True)
    if list(ET.parse(native/'response/generated_ncell_4.sonx').iter(work.tag('Component'))):
        raise ValueError('No ideal components allowed in selected response')
    run=json.loads((out/'response/run_status.json').read_text())
    for suffix,key in [('.sonx','project_sha256'),('.snp','snp_sha256')]:
        if work.digest(native/'response'/('generated_ncell_4'+suffix))!=run[key]:
            raise ValueError('Native response identity changed')
    log=(out/'response/em.log').read_text()
    if not all(x in log for x in ['Em errors detected: 0','Em warnings detected: 0','Em simulation completed']):
        raise ValueError('Native completion/warnings check failed')
    fitted=json.loads((out/'response'/fit_folder/'fit.json').read_text())
    if fitted['provenance']['sha256']!=run['snp_sha256'] or abs(fitted['L_nH'][0]-measured['L0_nH'])>1e-10:
        raise ValueError('Wrong fit source or T3 L')
    if fitted['L_nH'][1:]!=work.REFERENCE_FIT['L_nH']:
        raise ValueError('Older inductances changed')
    older_actual=np.array(fitted['C_fF'][1:]+fitted['G_fF'][1 if anchored_all_older else 2:])
    older_reference=np.array(work.REFERENCE_FIT['C_fF']+work.REFERENCE_FIT['G_fF'][0 if anchored_all_older else 1:])
    older_change=100*(older_actual/older_reference-1)
    if max(abs(older_change))>(bound_percent+1e-7 if bounded else 1e-10):
        raise ValueError('Older capacitances exceed selected fit policy')
    values=np.array([fitted['C_fF'][0],fitted['G_fF'][0],fitted['G_fF'][1]])
    if max(abs(100*(values/[251,261,153]-1)))>5 or fitted['complex_RMS']>.025:
        raise ValueError('Actual fitted values fail component/curve limits')
    if anchored_all_older:
        expected_free=[['C',0],['G',0]]+[[k,i] for i in range(1,5) for k in ('C','G')]
        junction_bound=result.get('junction_bound_percent',bound_percent)
        if junction_bound==0:
            expected_free.remove(['G',1])
            if fitted['G_fF'][1]!=work.REFERENCE_FIT['G_fF'][0]:raise ValueError('Fixed reference junction changed')
        if abs(100*(fitted['G_fF'][1]/work.REFERENCE_FIT['G_fF'][0]-1))>junction_bound+1e-7:
            raise ValueError('Junction prior exceeded')
        if fitted['free_capacitances'] != expected_free:
            raise ValueError('Unexpected anchored fit parameters')
        if fitted['provenance']['initial_C_fF'] != [251.]+work.REFERENCE_FIT['C_fF'] or fitted['provenance']['initial_G_fF'] != [261.]+work.REFERENCE_FIT['G_fF']:
            raise ValueError('Anchored references changed')
        if bound_percent >= 2:
            previous=json.loads((out/'response/fit_bounded_1pct/fit.json').read_text())
            if previous['complex_RMS']<=.025 and not result.get('refinement_authorized',False):
                raise ValueError('Wider bounds need failed 1% RMS or explicit refinement authorization')
    if not fitted['optimizer_success'] or any(fitted['active_bounds'][:2 if anchored_all_older else 3]+fitted['active_bounds'][-1:]):
        raise ValueError('New values or delay are at a search boundary')
    if fixed_junction:
        from geometry.taper.taper3_fixed_junction import verify_junction
        verify_junction(case)
        if fitted['G_fF'][1]!=work.REFERENCE_FIT['G_fF'][0] or fitted['free_capacitances']!=[['C',0],['G',0]]:
            raise ValueError('Fixed junction or two-capacitance fit policy changed')
    source_case=measured.get('reused_from_case',case)
    proof_path=work.OUT/source_case/'Lcal_explicit/completed.json'
    proof=json.loads(proof_path.read_text())
    for path,hashed in proof['files_sha256'].items():
        if work.digest(path)!=hashed:raise ValueError('Explicit L source changed')
    source_meta=json.loads((work.NATIVE/source_case/'Lcal/geometry.json').read_text())
    if source_meta['unit']!=meta['unit']:raise ValueError('Calibrated body differs from production')
    plot_case(case,fit_folder=fit_folder)
    destination=work.ROOT/'sonnet_files/taper_exploration/final_designs'/f'new_T3_{case}_final'
    base=(work.ROOT/'sonnet_files/taper_exploration').resolve()
    if not native.resolve().is_relative_to(base) or not destination.resolve().is_relative_to(base) or destination.exists():
        raise ValueError('Unsafe or existing final destination')
    destination.parent.mkdir(exist_ok=True);native.rename(destination)
    shutil.copytree(out,destination/'evidence')
    plot_source=out/'response'/fit_folder if alternate else out/'response'
    for ext in ('png','pdf'):shutil.copy2(plot_source/f's21_fit_parameters.{ext}',destination/f'final_s21_fit.{ext}')
    shutil.copy2(out/'response/geometry.png',destination/'geometry.png')
    source_cal=work.NATIVE/source_case/'Lcal_explicit'
    if source_case==case:source_cal=destination/'Lcal_explicit'
    if source_case!=case:
        shutil.copytree(source_cal,destination/'Lcal_explicit',ignore=shutil.ignore_patterns('sondata','*.cache'))
        shutil.copytree(work.NATIVE/source_case/'Lcal',destination/'calibration_source/Lcal')
    shutil.copytree(work.OUT/source_case/'Lcal_explicit',destination/'Lcal_explicit/evidence')
    # Preserve original provenance documents, and separately map their old
    # absolute paths to verified archive copies after moving the native case.
    archived=[]
    for original,hashed in proof['files_sha256'].items():
        old=Path(original)
        if old.is_relative_to(native):
            target=destination/old.relative_to(native)
        elif old.is_relative_to(work.NATIVE/source_case/'Lcal'):
            target=destination/'calibration_source/Lcal'/old.relative_to(work.NATIVE/source_case/'Lcal')
        elif old.is_relative_to(work.NATIVE/source_case/'Lcal_explicit'):
            target=destination/'Lcal_explicit'/old.relative_to(work.NATIVE/source_case/'Lcal_explicit')
        elif old.is_relative_to(out):
            target=destination/'evidence'/old.relative_to(out)
        elif old.is_relative_to(work.OUT/source_case/'Lcal_explicit'):
            target=destination/'Lcal_explicit/evidence'/old.relative_to(work.OUT/source_case/'Lcal_explicit')
        else:
            raise ValueError('Unknown calibration archive source')
        if work.digest(target)!=hashed:
            raise ValueError('Archived calibration hash mismatch')
        archived.append(dict(original=str(old),archive=str(target),sha256=hashed))
    (destination/'calibration_archive_manifest.json').write_text(json.dumps(archived,indent=2)+'\n')
    result = dict(result,project=str(destination/'response/generated_ncell_4.sonx'))
    p=work.from_dict(meta['t3_settings']);_,body,_=work.isolated_body(p.body,d_top=p.d_top,pad_right=p.pad_right)
    lib=gdstk.Library(unit=1e-6,precision=1e-9);top=lib.new_cell('TOP');top.add(gdstk.Reference(body));lib.add(body,*body.dependencies(True));lib.write_gds(str(destination/'t3_unit.gds'))
    selection=dict(stage='T3',case=case,T3_accepted=True,settings=meta['t3_settings'],result=result,
        inductance=measured,T2_reference=str(work.REFERENCE_FILE),T2_reference_sha256=work.REFERENCE_HASH,
        older_geometry_preserved=True,older_electrical_values_fixed=not bounded,
        fit_folder=fit_folder,electrical_T3_T2_junction_fixed=bool(anchored_all_older and result.get('junction_bound_percent')==0),older_capacitance_bound_percent=bound_percent if bounded else 0,fixed_T3_T2_junction=fixed_junction,older_capacitance_change_percent=older_change.tolist(),
        calibration_archive_manifest=str(destination/'calibration_archive_manifest.json'),
        acceptance='New Cs within5%, fit-to-Sonnet RMS<=.025, explicit L within2%; not a nominal-target S21 claim.',
        project=str(destination/'response/generated_ncell_4.sonx'),plot=str(destination/'final_s21_fit.png'),
        geometry_plot=str(destination/'geometry.png'),unit_gds=str(destination/'t3_unit.gds'))
    for path in [destination/'selection.json',work.ROOT/'geometry/taper/selected_taper3.json']:
        path.write_text(json.dumps(selection,indent=2)+'\n')
    progress=work.ROOT/'geometry/taper/taper3_progress.json';d=json.loads(progress.read_text());d.update(state='accepted',accepted=True,selection=selection);progress.write_text(json.dumps(d,indent=2)+'\n')
    policy=f'Older Cs bounded to +/-{bound_percent}% of selected T2; all L fixed.' if bounded else 'Only new Cs and one fixture delay fitted; older selected T2 model fixed.'
    (destination/'README.txt').write_text('Selected T3 under the component-matching policy.\nProject: response/generated_ncell_4.sonx\nFinal plot: final_s21_fit.png\nGeometry preview: geometry.png\nT3 body only: t3_unit.gds\n'+policy+'\nIndependent L uses explicit native Sonnet capacitors.\nNative2um mesh and400um reference planes preserved.\nArchive calibration file paths/hashes: calibration_archive_manifest.json\nNot a claim of full-filter or fixed-target S21 agreement.\n')
    print(json.dumps(selection,indent=2),flush=True)
    return selection
