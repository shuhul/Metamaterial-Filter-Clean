"""Parameter-driven physical T2 trials, independent L first, then full assembly.

No hidden optimizer releases downstream constants. Each case gets its own GDS,
native projects, calibration, complex fit, target overlay and recorded errors.
"""
from dataclasses import replace
import json
import xml.etree.ElementTree as ET
import gdstk
import numpy as np
from geometry.taper.progressive_fixtures import ROOT,OUT,NATIVE,GDS,prepare as prepare_output,run as run_output
from geometry.taper.compact_stages import build,prepare as prepare_compact
from geometry.taper.second_layout import from_dict
from geometry.taper.sonnet import add_capacitor
from geometry.prepare_sonnet import tag,write_xml
from geometry.taper.run_bridged_stage import digest
from geometry.taper.measure_compact_inductance import calibrate
from geometry.taper.progressive_fitting import fit
from geometry.taper.staged_fitting import circuit,LT1,LR


def parameters(frame=700,wire=68,pairs=36,span=266,gap=2,left=600,right=584,junction=6):
    source=ROOT/'geometry/taper/output/compact_centered/plain_w700_wire68_idc26_r584_t24_gap6/idc/geometry.json'
    p=from_dict(json.loads(source.read_text())['second'])
    m=replace(p.taper.meander,t=wire,w=max(4*(3*wire//4+1),160))
    return replace(p,taper=replace(p.taper,w=frame,l_left=left,l_right=right,meander=m),
                   idc_pairs=pairs,idc_span=span,idc_d=gap,d_t2_t1=junction)


def source_gds(name,p):
    folder=GDS/name/'source';folder.mkdir(parents=True,exist_ok=True)
    path=folder/'geometry.gds';polys,meta,_=build(p,'idc')
    if path.exists():
        if json.loads(path.with_suffix('.json').read_text())!=meta:raise ValueError('Changed existing trial')
    else:
        lib=gdstk.Library(unit=1e-6,precision=1e-9);lib.new_cell('TOP').add(*polys);lib.write_gds(str(path))
        path.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
    return path


def prepare(name,p):
    source=source_gds(name,p)
    project,out=prepare_output(name,str(source),3,302)
    return project,out


def measure(name,p,*,calibration_stage='idc_Lcal'):
    project,out,meta=prepare_compact('progressive_'+name,p,calibration_stage)
    marker=out/'operator_prepared.json'
    if not marker.exists():
        add_capacitor(project,meta);tree=ET.parse(project);comp=next(tree.iter(tag('Component')))
        comp.set('Type','PORTS');comp.attrib.pop('Units',None)
        for child in list(comp):
            if child.tag==tag('Value'):comp.remove(child)
        sweeps=tree.find(tag('Sweeps'))
        for child in list(sweeps):sweeps.remove(child)
        block=ET.SubElement(sweeps,tag('Set'),On='TRUE');fs=ET.SubElement(block,tag('Frequencies'))
        for f in (4,5,5.5,6,7):ET.SubElement(fs,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(f),Stop=str(f),Step='1')
        for q in tree.iter(tag('File')):
            if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','FALSE')
        write_xml(tree,project);marker.write_text(json.dumps(dict(project_sha256=digest(project)),indent=2)+'\n')
    elif json.loads(marker.read_text())['project_sha256']!=digest(project):raise ValueError('L operator changed')
    result=out/'results.json'
    if result.exists():return json.loads(result.read_text())
    return calibrate(project,out,name,auto_loads=True,desired_GHz=(5.,4.7,4.4,4.1))


def review(name,l,*,junction_cg=None):
    baseline=json.loads((OUT/'01_T1_regular/twoC/fit.json').read_text())
    junction=baseline['G_fF'][0] if junction_cg is None else junction_cg
    project=NATIVE/name/'generated_ncell_3.sonx';out=OUT/name
    r=fit(project.with_suffix('.snp'),out/'twoC',[l,LT1,LR],
          [100,baseline['C_fF'][0],190.4],[459,junction,8.2865],
          [('C',0),('G',0)],title='Taper 2 Transmission: Full Assembly')
    errors=dict(L_T2=100*(l/.74-1),C0_T2=100*(r['C_fF'][0]/100-1),Cg_input_T2=100*(r['G_fF'][0]/459-1),
                Cg_T1=100*(junction/56.7-1))
    r['target_errors_percent']=errors
    r['passes_parameter_targets']=all(abs(v)<=2 for v in errors.values())
    data=np.load(out/'twoC/response.npz');f=data['f'];actual=data['actual'];model=data['model']
    target=circuit(f,[.74,6.06,4.9194],[100,100,190.4],[459,56.7,8.2865],left_ps=r['delay_ps'])
    r['target_S21_magnitude_RMS']=float(np.sqrt(np.mean((abs(target[:,1,0])-abs(actual[:,1,0]))**2)))
    r['target_S21_magnitude_max']=float(np.max(abs(abs(target[:,1,0])-abs(actual[:,1,0]))))
    r['passes_target_response']=r['target_S21_magnitude_RMS']<.02 and r['target_S21_magnitude_max']<.05
    r['final_acceptance_pending']=True
    r['independent_L_verified_for_this_geometry']=not (out/'provisional_inductance.json').exists()
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(12,7))
    ax.plot(f/1e9,abs(actual[:,1,0]),color='#b22222',lw=3,label='Sonnet simulation')
    ax.plot(f/1e9,abs(model[:,1,0]),'k--',lw=2.3,label='Complex S21 fit')
    ax.plot(f/1e9,abs(target[:,1,0]),color='#1565c0',ls=':',lw=2.3,label='Target circuit')
    text=(f"L T2 = {l:.5f} nH ({errors['L_T2']:+.2f}%)\n"
          f"C0 T2 = {r['C_fF'][0]:.2f} fF ({errors['C0_T2']:+.2f}%)\n"
          f"Cg T2 = {r['G_fF'][0]:.2f} fF ({errors['Cg_input_T2']:+.2f}%)\n"
          f"Cg T1 = {junction:.3f} fF ({errors['Cg_T1']:+.2f}%; fixed)\nComplex S21 RMS = {r['complex_RMS']:.5f}")
    ax.text(.03,.97,text,transform=ax.transAxes,va='top',fontsize=15,bbox=dict(boxstyle='round,pad=.5',fc='white',ec='gray'))
    ax.set_ylim(0,max(abs(actual[:,1,0]).max(),abs(target[:,1,0]).max())*1.6)
    ax.set_title('Exploratory Taper 2 Transmission',fontsize=24,pad=15);ax.set_xlabel('Frequency (GHz)',fontsize=19);ax.set_ylabel(r'$|S_{21}|$',fontsize=24)
    ax.tick_params(labelsize=16);ax.legend(fontsize=13);ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(out/'s21_fit_linear.png',dpi=200);plt.close(fig)
    (out/'target_review.json').write_text(json.dumps(r,indent=2)+'\n')
    print('TARGET REVIEW',name,errors,r['complex_RMS'],r['target_S21_magnitude_RMS'],flush=True)
    return r


def run(name,p):
    project,out=prepare(name,p)
    cal=measure(name,p)
    if cal['R_squared']<.9999:raise ValueError('L calibration is not linear enough')
    # Early C-geometry trials may miss L slightly; fit using measured L and
    # keep the strict 2% final target gate in review(). Avoid recalibrating a
    # wire repeatedly before learning which IDC/frame geometry is needed.
    if abs(cal['L0_nH']/.74-1)>.10:raise ValueError('Independent L misses exploratory 10% guard')
    source=source_gds(name,p)
    run_output(name,str(source),3,302)
    return review(name,cal['L0_nH'])


def screen_with_reference_L(name,p,l,reference):
    """One controlled coupling-spacing screen; L is provisional, never final."""
    project,out=prepare(name,p)
    note=dict(L_estimate_nH=l,reference=str(reference),reference_sha256=digest(reference),
              independent_L_verified_for_this_geometry=False,
              reason='Only junction separation changed; use prior L for a geometry sensitivity screen, then recalibrate selected final geometry.')
    (out/'provisional_inductance.json').write_text(json.dumps(note,indent=2)+'\n')
    run_output(name,str(source_gds(name,p)),3,302)
    result=review(name,l);result['independent_L_verified_for_this_geometry']=False
    (out/'target_review.json').write_text(json.dumps(result,indent=2)+'\n')
    return result


if __name__=='__main__':
    import sys
    name='trial_w700_wire68_idc36_span266_gap2_left600'
    p=parameters()
    if '--prepare-only' in sys.argv:prepare(name,p)
    else:run(name,p)
