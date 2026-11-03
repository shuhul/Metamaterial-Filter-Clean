"""Capacitance heuristics from saved full two-port data, then parameter-only EM.

Electrical fitting is a guide, never the fixed-target acceptance test. All
inductances and regular-cell parameters are fixed to documented measurements.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
from dataclasses import replace
import json
import numpy as np
from scipy.optimize import least_squares
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from circuit_common import read_snp
from geometry.taper.plain_gap_study import ROOT, BASE_NATIVE, BASE_OUT, baseline
from geometry.taper.run_bridged_stage import digest
from geometry.taper.staged_fitting import LT1, LR
from geometry.taper.corrected_iteration import REGULAR_C, REGULAR_G
from geometry.taper.progressive_fitting import transmission

CAMPAIGN='capacitance_guided_2026-10-25'
OUT=ROOT/'output_data/taper_exploration'/CAMPAIGN
NATIVE=ROOT/'sonnet_files/taper_exploration'/CAMPAIGN
LABELS=['C0_T2_fF','Cg_input_T2_fF','C0_T1_fF','Cg_T2_T1_fF','Cg_T1_unit_fF']
L_BASE=.7342278563325565


def network(f,values,l=L_BASE):
    c0,g0,c1,g1,g2,left,right=values
    ls=[l,LT1,LR,LR];cs=[c0,c1,REGULAR_C,REGULAR_C];gs=[g0,g1,g2,REGULAR_G]
    w=2*np.pi*f;a=np.ones(len(f),complex);b=np.zeros(len(f),complex);c=b.copy();d=a.copy()
    for ll,cc,gg in zip(ls,cs,gs):
        z=1/(1j*w*gg*1e-15);b=b+a*z;d=d+c*z
        y=1j*w*cc*1e-15+1/(1j*w*ll*1e-9);a=a+b*y;c=c+d*y
    den=a+b/50+50*c+d
    return np.column_stack([(a+b/50-50*c-d)/den*np.exp(-2j*w*left*1e-12),
        2/den*np.exp(-1j*w*(left+right)*1e-12),
        (-a+b/50-50*c+d)/den*np.exp(-2j*w*right*1e-12)])


def extract(source,folder,l=L_BASE):
    folder=Path(folder);folder.mkdir(parents=True,exist_ok=True)
    provenance=dict(source=str(source),source_sha256=digest(source),code_sha256=digest(__file__),
        fixed_L_nH=[l,LT1,LR,LR],fixed_regular_C_fF=REGULAR_C,fixed_output_Cg_fF=REGULAR_G,
        free=LABELS+['left_fixture_delay_ps','right_fixture_delay_ps'],
        fitted_channels=['complex S11','complex S21','complex S22'])
    marker=folder/'extraction.json'
    if marker.exists():
        saved=json.loads(marker.read_text())
        if saved['provenance']==provenance:return saved
    net=read_snp(source);f=net.f;y=np.column_stack([net.s[:,0,0],net.s[:,1,0],net.s[:,1,1]])
    ix=np.unique(np.linspace(0,len(f)-1,2001).astype(int))
    lower=np.r_[np.log([.01,10,20,1,2]),[-50,-50]]
    upper=np.r_[np.log([1500,6000,400,300,20]),[50,50]]
    def expand(x):return np.r_[np.exp(x[:5]),x[5:]]
    def residual(x,indices=ix):
        delta=network(f[indices],expand(x),l)-y[indices]
        return np.r_[delta.real.ravel(),delta.imag.ravel()]
    def smooth(x):
        m=network(f[ix],expand(x),l)[:,1];s=y[ix,1]
        delta=(1/m-1/s)*np.maximum(abs(s),.03)**2
        return np.r_[delta.real,delta.imag]
    seeds=[[100,459,100,56.7,8.2865,0,0],
           [10,1200,150,15,8.2,-10,-10],
           [200,400,120,30,8.2,0,0],
           [400,800,100,56.7,8.2,0,0]]
    starts=[np.r_[np.log(v[:5]),v[5:]] for v in seeds]
    smooths=[least_squares(smooth,x,bounds=(lower,upper),max_nfev=180,x_scale=np.r_[np.ones(5),20,20]).x for x in starts]
    trials=[]
    for x in starts+smooths:
        q=least_squares(residual,x,bounds=(lower,upper),max_nfev=300,
                        x_scale=np.r_[np.ones(5),20,20],ftol=1e-10,xtol=1e-10,gtol=1e-10)
        trials.append(q)
    q=min(trials,key=lambda q:q.fun@q.fun)
    q=least_squares(lambda x:residual(x,np.arange(len(f))),q.x,bounds=(lower,upper),
                    max_nfev=200,x_scale=np.r_[np.ones(5),20,20],ftol=1e-11,xtol=1e-11,gtol=1e-11)
    values=expand(q.x);model=network(f,values,l);rms=np.sqrt(np.mean(abs(model-y)**2,axis=0))
    singular=np.linalg.svd(q.jac*np.r_[np.ones(5),20,20],compute_uv=False)
    near=[expand(t.x).tolist() for t in trials if
          np.linalg.norm(t.fun)/np.sqrt(len(t.fun))<=np.linalg.norm(q.fun)/np.sqrt(len(q.fun))*1.05]
    result=dict(provenance=provenance,values=dict(zip(LABELS,values[:5].tolist())),
        fixture_delay_ps=values[5:].tolist(),complex_RMS=dict(zip(['S11','S21','S22'],rms.tolist())),
        all_channel_complex_RMS=float(np.sqrt(np.mean(abs(model-y)**2))),
        optimizer_success=bool(q.success),active_bounds=q.active_mask.tolist(),
        scaled_jacobian_singular_values=singular.tolist(),similar_residual_start_solutions=near,
        scope='Full two-port circuit heuristics, all L and regular cells fixed. '
              'Imperfect circuit structure or fixture model can bias fitted capacitances; no uniqueness or physical certification implied.')
    marker.write_text(json.dumps(result,indent=2)+'\n')
    np.savez(folder/'extraction.npz',f=f,actual=y,model=model)
    fig,axs=plt.subplots(1,3,figsize=(18,5))
    for i,ax in enumerate(axs):
        ax.plot(f/1e9,abs(y[:,i]),lw=2,label='Sonnet');ax.plot(f/1e9,abs(model[:,i]),'--',lw=2,label='Circuit heuristic')
        ax.set_title(['S11','S21','S22'][i]+f' · RMS {rms[i]:.4f}',fontsize=16)
        ax.set_xlabel('Frequency (GHz)',fontsize=13);ax.set_ylabel('Magnitude',fontsize=13);ax.grid(alpha=.2);ax.legend()
    fig.suptitle('Capacitance Extraction: Independent Reflection Information Included',fontsize=20)
    fig.tight_layout();fig.savefig(folder/'extraction.png',dpi=160);plt.close(fig)
    print('EXTRACTION',str(folder),json.dumps(result),flush=True)
    return result


def saved_study():
    sources={'baseline':BASE_NATIVE.with_suffix('.snp')}
    for name in ('input_span156','facing_thickness24','frame_width416'):
        sources[name]=ROOT/'sonnet_files/taper_exploration/plain_gap_sensitivity_2026-10-25'/name/'idc/generated_ncell_3.snp'
    results={name:extract(source,OUT/name/'extraction') for name,source in sources.items()}
    (OUT/'saved_geometry_extractions.json').write_text(json.dumps(results,indent=2)+'\n')
    return results


def candidate(name,p,*,calibration_case=None,record_currents=False):
    """One sequential guarded EM job; preserve original T1/downstream metal."""
    import xml.etree.ElementTree as ET
    import gdstk
    from geometry.prepare_sonnet import tag,write_xml
    from geometry.taper.clearance_extraction import prepare
    from geometry.taper.sonnet import configure,simulate
    from geometry.taper.target_iteration import START,STOP,score
    from geometry.taper.response_diagnostics import analyze
    cal=OUT/(calibration_case or name)/'idc_Lcal/results.json'
    if calibration_case:
        from geometry.taper.second_layout import from_dict
        data=ROOT/'geometry/taper/output'/CAMPAIGN/calibration_case/'idc_Lcal/geometry.json'
        cp=from_dict(json.loads(data.read_text())['second'])
        frame_wire=lambda q:(q.taper.w,q.taper.t_top,q.taper.meander.t,q.taper.meander.l,q.straight_inductor,q.x_meander)
        if frame_wire(cp)!=frame_wire(p):
            raise ValueError('Shared calibration requires exactly the same frame/bar and straight wire')
    inductance=json.loads(cal.read_text())['L0_nH'] if cal.exists() else L_BASE
    if (p.taper.w,p.taper.meander.t)!=(400,60) and not cal.exists():
        raise ValueError('Changed frame/wire needs its own inductance calibration first')
    if abs(inductance/.74-1)>.02:raise ValueError('Retune wire width before the full response: L outside 2%')
    if not p.equal_idc_sides or not p.straight_inductor or p.x_meander:
        raise ValueError('Preserve uniform IDC and centered straight inductor')
    if p.taper1!=baseline().taper1:raise ValueError('T1 shape must stay original')
    if (p.t_coupling_finger or p.w_upper_plate or p.w_t1_ground_slot or p.w_ground_pad or
        p.t_inner_ground or p.l_ground_rise or p.t_ground_tail is not None):
        raise ValueError('No added structures or ground-plane cuts')
    project,out,meta=prepare(name,p,'idc',campaign=CAMPAIGN)
    gd=ROOT/'geometry/taper/output'/CAMPAIGN/name/'idc/geometry.gds'
    metal=gdstk.read_gds(str(gd)).top_level()[0].get_polygons()
    if gdstk.boolean([gdstk.rectangle((0,0),(meta['w_sim'],448))],metal,'not',precision=.001):
        raise ValueError('Top ground must be the complete straight rectangle')
    from dataclasses import asdict
    from geometry.taper.plain_gap_study import leaves
    old=leaves(asdict(baseline()));new=leaves(asdict(p))
    changes={k:dict(baseline=old[k],candidate=v) for k,v in new.items() if old[k]!=v}
    intent=dict(topology='Approved plain-gap layout; no added structures',
        T1_downstream_native_lock='passed',changed_parameters=changes,
        rationale='Parameter-only iteration guided by extracted loading and measured response sensitivities; see campaign README.',
        IDC_span_um=p.idc_span,T2_T1_gap_um=p.d_t2_t1,inductance_for_heuristic_nH=inductance,
        inductance_scope=('Same frame/wire reference; recheck full final geometry' if calibration_case else 'Own isolated calibration') if cal.exists() else 'Baseline isolated calibration; recheck final candidate',
        calibration_reference=str(cal) if cal.exists() else None,
        existing_junction_ground_offset_um=p.x_junction_ground,
        native_lock_exception='Existing junction gap only; cell shapes still unchanged' if p.x_junction_ground else None,
        record_currents=record_currents)
    (out/'design_intent.json').write_text(json.dumps(intent,indent=2)+'\n')
    marker=out/'run_status.json'
    if marker.exists() and json.loads(marker.read_text()).get('status')=='complete':
        record=json.loads(marker.read_text())
        if record['project_sha256']!=digest(project) or record['snp_sha256']!=digest(project.with_suffix('.snp')):
            raise ValueError('Completed candidate changed')
        source=project.with_suffix('.snp')
    else:
        tree=ET.parse(project);root=tree.getroot();sweeps=root.find(tag('Sweeps'))
        if record_currents:
            control=root.find(tag('Control'));current=control.find(tag('Currents'))
            if current is None:current=ET.SubElement(control,tag('Currents'))
            current.set('On','TRUE')
        for child in list(sweeps):sweeps.remove(child)
        fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
        ET.SubElement(fs,tag('Sweep'),Type='ABS_ENTRY',On='TRUE',Adaptive='TRUE',Start=str(START),Stop=str(STOP),Target='7500')
        for q in root.iter(tag('File')):
            if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','TRUE')
        write_xml(tree,project);configure(project,START,STOP,7500);source=simulate(project,out)
    result=score(source,out);analyze(source,out)
    intent['fixed_target_complex_RMS']=result['complex_RMS']
    (out/'design_intent.json').write_text(json.dumps(intent,indent=2)+'\n')
    return result


if __name__=='__main__':
    import sys
    if '--candidate' in sys.argv:
        name=sys.argv[sys.argv.index('--candidate')+1]
        span,gap=map(int,name.removeprefix('span').split('_gap'))
        candidate(name,replace(baseline(),idc_span=span,d_t2_t1=gap))
    else:saved_study()
