"""New three-stage targets; tune T1 using the regular cell's modular parts.

All settings are Python parameters. No older accepted taper geometry is used.
Each candidate has immutable GDS, native Sonnet, calibration and response files.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from dataclasses import asdict, replace
from pathlib import Path
from types import SimpleNamespace
import json
import shutil
import xml.etree.ElementTree as ET
import gdstk
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from scipy.optimize import minimize_scalar, least_squares
from geometry.capacitor import CapParams
from geometry.meander import MeanderParams
from geometry.pad import PadParams
from geometry.unit_cell import UnitCellParams, make_unit_cell
from geometry.specifications import specification, TAPER_1, UNIT_CELL, S21_COMPLEX_RMS_CUTOFF
from geometry.prepare_sonnet import main as import_native, tag, write_xml
from geometry.taper.sonnet import add_capacitor, configure, simulate
from geometry.taper.measure_compact_inductance import calibrate, W
from geometry.taper.progressive_fitting import transmission
from circuit_common import read_snp

ROOT=Path(__file__).resolve().parents[2]
CAMPAIGN='three_tapers_2026-10-25'
BASE=UnitCellParams(cap=CapParams(w=328),meander=MeanderParams(w=224),
                    pad=PadParams(w=8),d_pad=20)
REGULAR=UnitCellParams(meander=MeanderParams(w=186))


def from_dict(v):
    return UnitCellParams(cap=CapParams(**v['cap']),meander=MeanderParams(**v['meander']),
                          pad=PadParams(**v['pad']),d_top=v['d_top'],d_pad=v['d_pad'])


def build(p,*,isolation=False,d_input=40,d_junction=56):
    # Shared ground height and the original 2 um Sonnet mesh are preserved.
    if p.d_top+p.cap.t_top+p.meander.l!=618:
        raise ValueError('T1 must meet the flat bottom ground at -618 um')
    for v in (d_input,d_junction):
        if v<8 or v%4:raise ValueError('Grounded gaps must be multiples of 4 um and >=8 um')
    units=[p,p,p] if isolation else [p,p,REGULAR,REGULAR]
    gaps=[d_input,d_input] if isolation else [d_input,d_junction,56]
    lefts=[426]
    for u,gap in zip(units,gaps):lefts.append(lefts[-1]+u.cap.w+gap)
    centers=[x+u.cap.w/2 for x,u in zip(lefts,units)]
    width=lefts[-1]+units[-1].cap.w+426
    top=gdstk.Cell('TOP')
    for i,(u,x) in enumerate(zip(units,centers)):
        top.add(gdstk.Reference(make_unit_cell(u,name=f'CELL_{i}'),origin=(x,0)))
    parts=[gdstk.rectangle((0,0),(width,448)),
           gdstk.rectangle((0,-1052),(width,-618)),
           gdstk.rectangle((0,-1052),(412,-24)),
           gdstk.rectangle((width-412,-1052),(width,-24)),
           gdstk.rectangle((0,-18),(centers[0],-8)),
           gdstk.rectangle((centers[-1],-18),(width,-8))]
    for i,gap in enumerate(gaps):
        x=lefts[i]+units[i].cap.w+gap/2
        parts.append(gdstk.rectangle((x-2,-618),(x+2,0)))
    top.add(*gdstk.boolean(parts,[],'or',precision=.001,layer=1))
    meta={'layout_kind':'new_t1_from_unit_cell','unit':asdict(p),'isolation':isolation,
          'd_input':d_input,'d_junction':d_junction,'w_sim':width,'cell_centers':centers,
          'parameters':{'n_cells':1 if isolation else 2,'l_sim':1500,'l_top_ground':448,
                        'w_side':412,'t_feed':10,'unit':{'d_top':8,'cap':{'t_top':14}}}}
    return top,meta


def validate(path,meta):
    top,rebuilt=build(from_dict(meta['unit']),isolation=meta['isolation'],
                      d_input=meta['d_input'],d_junction=meta['d_junction'])
    if rebuilt!=meta:raise ValueError('Geometry metadata changed')
    lib=gdstk.read_gds(str(path));tops=lib.top_level()
    if len(tops)!=1 or tops[0].name!='TOP' or not np.isclose(lib.unit,1e-6):raise ValueError('Invalid GDS')
    polys=tops[0].get_polygons()
    if any(q.layer!=1 or q.datatype!=0 for q in polys):raise ValueError('Invalid metal layer')
    if gdstk.boolean(polys,top.get_polygons(),'xor',precision=.001):raise ValueError('GDS differs from parameters')
    # Ground strips join upper/lower ground; every resonator joins lower
    # ground through its inductor, as in the original regular n-cell fixture.
    if len(gdstk.boolean(polys,[],'or',precision=.001))!=1:raise ValueError('Wrong ground connectivity')
    return SimpleNamespace(n_cells=meta['parameters']['n_cells'],l_sim=1500,l_top_ground=448),polys


def prepare(name,p,*,isolation=False,d_input=40,d_junction=56,
            campaign=CAMPAIGN,builder=build,validator=validate):
    stage='Lcal' if isolation else 'response'
    gd=ROOT/'geometry/taper/output'/campaign/name/stage
    native=ROOT/'sonnet_files/taper_exploration'/campaign/name/stage
    out=ROOT/'output_data/taper_exploration'/campaign/name/stage
    for folder in (gd,native,out):folder.mkdir(parents=True,exist_ok=True)
    top,meta=builder(p,isolation=isolation,d_input=d_input,d_junction=d_junction)
    path=gd/'geometry.gds'
    if path.exists():
        if json.loads(path.with_suffix('.json').read_text())!=meta:raise ValueError('Case geometry changed')
    else:
        lib=gdstk.Library(unit=1e-6,precision=1e-9);lib.add(top,*top.dependencies(True));lib.write_gds(str(path))
        path.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
    validator(path,meta)
    fig,ax=plt.subplots(figsize=(15,6))
    for poly in top.get_polygons():ax.add_patch(Polygon(poly.points,fc='#297b54',ec='none'))
    ax.set(xlim=(360,meta['w_sim']-360),ylim=(-650,45),aspect='equal',xlabel='x (um)',ylabel='y (um)',
           title=f'New T1: {name} — '+('Isolated L Calibration' if isolation else 'T1 + Regular Cell'))
    fig.tight_layout();fig.savefig(out/'geometry.png',dpi=190);plt.close(fig)
    project=native/f"generated_ncell_{meta['parameters']['n_cells']}.sonx"
    if not project.exists():import_native(meta['parameters']['n_cells'],gds_path=path,output_dir=native,geometry_validator=validator)
    (out/'specification.json').write_text(json.dumps(specification(),indent=2)+'\n')
    return project,out,meta


def measure(name,p,*,d_input=40,prepare_case=prepare):
    project,out,meta=prepare_case(name,p,isolation=True,d_input=d_input)
    tree=ET.parse(project)
    if not list(tree.iter(tag('Component'))):
        add_capacitor(project,meta);tree=ET.parse(project)
        comp=next(tree.iter(tag('Component')));comp.set('Type','PORTS');comp.attrib.pop('Units',None)
        for q in list(comp):
            if q.tag==tag('Value'):comp.remove(q)
        sweeps=tree.find(tag('Sweeps'))
        for q in list(sweeps):sweeps.remove(q)
        fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
        for f in (1,2,2.5,3,4):
            ET.SubElement(fs,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(f),Stop=str(f),Step='1')
        for q in tree.findall(f"{tag('OutputFiles')}//{tag('File')}"):
            if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','FALSE')
        write_xml(tree,project)
    result=calibrate(project,out,name,auto_loads=True,target_L_nH=TAPER_1.L0_nH,
                     held_out_GHz=2.5,frequency_range_GHz=(1,4),desired_GHz=[3.8,3.4,3,2.6])
    return result


def response(name,p,l,*,d_input=40,d_junction=56,prepare_case=prepare,input_end_cell=True):
    project,out,meta=prepare_case(name,p,d_input=d_input,d_junction=d_junction)
    marker=out/'run_status.json'
    if marker.exists() and json.loads(marker.read_text()).get('status')=='complete':
        old=json.loads(marker.read_text())
        if old['project_sha256']!=W['digest'](project) or old['snp_sha256']!=W['digest'](project.with_suffix('.snp')):
            raise ValueError('Completed simulation changed')
        snp=project.with_suffix('.snp')
    else:
        configure(project,4.65,5.4,7500);snp=simulate(project,out)
    # All L and unit parameters are fixed. T1 C0, input Cg and junction Cg
    # may vary in the heuristic complex fit, along with fixture phase delay.
    net=read_snp(snp);f=net.f;obs=net.s[:,1,0]
    # V1 had an input-end LC; the corrected V2 electrode has none. The
    # directly fed regular output end stays in the model in both cases.
    def predict(x,frequencies,nominal=False):
        lr=UNIT_CELL.L0_nH if nominal else 4.922723945068658
        cr=UNIT_CELL.C0_fF if nominal else 188.8781163861302
        gout=UNIT_CELL.Cg_input_fF if nominal else 8.241620220231225
        lt=TAPER_1.L0_nH if nominal else l
        ls=[lt,lt,lr,lr] if input_end_cell else [lt,lr,lr]
        cs=[x[0],x[0],cr,cr] if input_end_cell else [x[0],cr,cr]
        gs=[None,x[1],x[2],gout] if input_end_cell else [x[1],x[2],gout]
        return transmission(frequencies,ls,cs,gs,x[3])
    ix=np.unique(np.linspace(0,len(f)-1,min(len(f),2001)).astype(int))
    def residual(x,indices=ix,smooth=False):
        pred=predict(x,f[indices]);actual=obs[indices]
        delta=(1/pred-1/actual)*np.maximum(abs(actual),.03)**2 if smooth else pred-actual
        return np.r_[delta.real,delta.imag]
    seed=np.array([TAPER_1.C0_fF,TAPER_1.Cg_input_fF,UNIT_CELL.Cg_input_fF,5.])
    bounds=([20,.1,.1,-100],[400,100,30,100]);scale=[150,12,8,20]
    starts=[seed]+[seed*np.r_[np.exp(np.random.default_rng(i).normal(0,.15,3)),1] for i in range(4)]
    coarse=[least_squares(lambda x:residual(x,smooth=True),x,bounds=bounds,x_scale=scale,max_nfev=250) for x in starts]
    fits=[least_squares(residual,q.x,bounds=bounds,x_scale=scale,max_nfev=400) for q in coarse]
    best=min(fits,key=lambda q:q.fun@q.fun)
    best=least_squares(lambda x:residual(x,np.arange(len(f))),best.x,bounds=bounds,x_scale=scale,max_nfev=400)
    model=predict(best.x,f)
    fitted={'C_fF':[float(best.x[0])],'G_fF':best.x[1:3].tolist(),'L_T1_nH':l,
            'delay_ps':float(best.x[3]),'complex_RMS':float(np.sqrt(np.mean(abs(model-obs)**2))),
            'success':bool(best.success),'fixed_unit_values':True,'input_end_LC_present':input_end_cell,
            'free_parameters':['C0_T1','Cg_input_T1','Cg_T1_unit','fixture_delay_ps'],
            'active_bounds':best.active_mask.tolist(),'source_sha256':W['digest'](snp)}
    singular=np.linalg.svd(best.jac*np.maximum(abs(best.x),[1,1,1,1])[None,:],compute_uv=False)
    fitted['relative_jacobian_condition']=float(singular[0]/singular[-1])
    fitted['multistart_branches']=[{'parameters':q.x.tolist(),
                                   'sampled_complex_RMS':float(np.sqrt(np.mean(abs(predict(q.x,f[ix])-obs[ix])**2)))} for q in fits]
    (out/'extraction').mkdir(exist_ok=True)
    (out/'extraction/fit.json').write_text(json.dumps(fitted,indent=2)+'\n')
    nominal=predict([TAPER_1.C0_fF,TAPER_1.Cg_input_fF,UNIT_CELL.Cg_input_fF,0],f,nominal=True)
    objective=lambda tau:float(np.mean(abs(nominal*np.exp(-2j*np.pi*f*tau*1e-12)-obs)**2))
    grid=np.linspace(-100,100,801);i=int(np.argmin([objective(t) for t in grid]))
    opt=minimize_scalar(objective,bounds=(max(-100,grid[i]-.5),min(100,grid[i]+.5)),method='bounded')
    target=nominal*np.exp(-2j*np.pi*f*opt.x*1e-12)
    score={'complex_RMS':float(np.sqrt(opt.fun)),'pass':bool(opt.success and np.sqrt(opt.fun)<=S21_COMPLEX_RMS_CUTOFF),
           'cutoff':S21_COMPLEX_RMS_CUTOFF,
           'fixture_delay_ps':float(opt.x),'all_electrical_targets_fixed':True,'specification':specification(),
           'source_sha256':W['digest'](snp),'frequency_range_GHz':[4.65,5.4]}
    (out/'target_comparison.json').write_text(json.dumps(score,indent=2)+'\n')
    np.savez(out/'extraction/response.npz',f=f,actual=net.s,model_S21=model,target_S21=target)
    values=[('L T1',l,TAPER_1.L0_nH,'nH'),('C0 T1',fitted['C_fF'][0],TAPER_1.C0_fF,'fF'),
            ('Cg T1',fitted['G_fF'][0],TAPER_1.Cg_input_fF,'fF'),
            ('Cg unit',fitted['G_fF'][1],UNIT_CELL.Cg_input_fF,'fF')]
    text='\n'.join(f'{label} = {v:.4g} {unit} ({100*(v/t-1):+.2f}%)' for label,v,t,unit in values)
    fig,ax=plt.subplots(figsize=(14,8))
    ax.plot(f/1e9,abs(obs),lw=2.8,label='Sonnet EM',color='#2768aa')
    ax.plot(f/1e9,abs(model),'k--',lw=2.2,label='Complex circuit fit')
    ax.plot(f/1e9,abs(target),':',lw=2.6,label='Specified target circuit',color='#db6d28')
    ax.text(.025,.97,text+f"\nTarget RMS = {score['complex_RMS']:.5f} (limit {S21_COMPLEX_RMS_CUTOFF:g})\nFit RMS = {fitted['complex_RMS']:.5f}",
            transform=ax.transAxes,va='top',fontsize=16,bbox=dict(boxstyle='round,pad=.5',fc='white',ec='gray'))
    ax.set(xlabel='Frequency (GHz)',ylabel=r'$|S_{21}|$',title='First Taper Cell: Transmission',ylim=(0,1.65),xlim=(4.65,5.4))
    ax.xaxis.label.set_size(20);ax.yaxis.label.set_size(23);ax.title.set_size(23);ax.tick_params(labelsize=16)
    ax.grid(alpha=.2);ax.legend(fontsize=15,loc='upper right');fig.tight_layout()
    for suffix in ('png','pdf'):fig.savefig(out/f's21_fit.{suffix}',dpi=190)
    plt.close(fig)
    print('NEW T1',name,'TARGET RMS',score['complex_RMS'],'FIT',fitted['C_fF'],fitted['G_fF'],flush=True)
    return score,fitted


def finalize(name,*,campaign=CAMPAIGN,geometry_validator=validate):
    """Promote only a measured, passing T1 and keep the complete native case."""
    from geometry.taper.audit import check_geometry
    out=ROOT/'output_data/taper_exploration'/campaign/name
    score=json.loads((out/'response/target_comparison.json').read_text())
    fitted=json.loads((out/'response/extraction/fit.json').read_text())
    measured=json.loads((out/'Lcal/results.json').read_text())
    if not score['pass'] or score['complex_RMS']>S21_COMPLEX_RMS_CUTOFF or not fitted['success'] or fitted['complex_RMS']>S21_COMPLEX_RMS_CUTOFF:
        raise ValueError('T1 does not pass target and fitted-curve checks')
    if abs(100*(measured['L0_nH']/TAPER_1.L0_nH-1))>2:
        raise ValueError('T1 isolated inductance is outside tolerance')
    if not np.isclose(fitted['L_T1_nH'],measured['L0_nH'],rtol=1e-12):
        raise ValueError('The fit used a different isolated inductance')
    if score['specification']!=specification():raise ValueError('Results use a superseded target set')
    source=(ROOT/'sonnet_files/taper_exploration'/campaign/name).resolve()
    destination=(ROOT/'sonnet_files/taper_exploration/final_designs'/f'new_T1_{name}_final').resolve()
    if not source.is_relative_to(ROOT) or not destination.is_relative_to(ROOT):
        raise ValueError('Finalization paths must stay inside the workspace')
    if destination.exists():raise ValueError('Final directory already exists')
    for stage,n in (('Lcal',1),('response',2)):
        gds=ROOT/'geometry/taper/output'/campaign/name/stage/'geometry.gds'
        meta=json.loads(gds.with_suffix('.json').read_text())
        check_geometry(source/stage/f'generated_ncell_{n}.sonx',gds,meta,geometry_validator=geometry_validator)
    calibration_meta=json.loads((ROOT/'geometry/taper/output'/campaign/name/'Lcal/geometry.json').read_text())
    response_meta=json.loads((ROOT/'geometry/taper/output'/campaign/name/'response/geometry.json').read_text())
    if calibration_meta['unit']!=response_meta['unit']:
        raise ValueError('Calibration core differs from the response geometry')
    calibration_marker=json.loads((out/'Lcal/operator_completed.json').read_text())
    calibration_project=source/'Lcal/generated_ncell_1.sonx'
    if (W['digest'](calibration_project)!=calibration_marker['project_sha256'] or
        W['digest'](calibration_project.with_suffix('.snp'))!=calibration_marker['snp_sha256']):
        raise ValueError('Calibration source identity changed')
    for log in (source/'Lcal/em.log',out/'response/em.log'):
        text=log.read_text()
        if 'Em errors detected: 0' not in text or 'Em warnings detected: 0' not in text or 'Em simulation completed' not in text:
            raise ValueError('Incomplete or flagged native simulation')
    project=source/'response/generated_ncell_2.sonx';snp=project.with_suffix('.snp')
    run=json.loads((out/'response/run_status.json').read_text())
    if run['status']!='complete' or W['digest'](project)!=run['project_sha256'] or W['digest'](snp)!=score['source_sha256']:
        raise ValueError('Simulation identity check failed')
    if list(ET.parse(project).iter(tag('Component'))):raise ValueError('Production geometry contains ideal components')
    destination.parent.mkdir(exist_ok=True)
    source.rename(destination)
    for file in ('s21_fit.png','s21_fit.pdf','geometry.png','target_comparison.json'):
        shutil.copy2(out/'response'/file,destination/file)
    gds=ROOT/'geometry/taper/output'/campaign/name/'response/geometry.gds'
    for suffix in ('.gds','.json'):shutil.copy2(gds.with_suffix(suffix),destination/('geometry'+suffix))
    unit=make_unit_cell(from_dict(response_meta['unit']),name='NEW_T1')
    library=gdstk.Library(unit=1e-6,precision=1e-9)
    library.new_cell('TOP').add(gdstk.Reference(unit));library.add(unit,*unit.dependencies(True))
    library.write_gds(str(destination/'t1_unit.gds'))
    selected={'specification':specification(),'stage':'T1','case':name,'T1_accepted':True,
              'T2_accepted':False,'T3_accepted':False,'geometry':json.loads(gds.with_suffix('.json').read_text()),
              'isolated_L':measured,'fixed_target_comparison':score,'heuristic_fit':fitted,
              'fixed_regular_extraction_model':{'L0_nH':4.922723945068658,
                  'C0_fF':188.8781163861302,'output_Cg_fF':8.241620220231225},
              'project':str(destination/'response/generated_ncell_2.sonx'),
              'plot':str(destination/'s21_fit.png'),'gds':str(destination/'geometry.gds'),
              'unit_gds':str(destination/'t1_unit.gds'),
              'audit':{'GDS_native_XOR_empty':True,'mesh_um':[2,2],'materials_preserved':True,
                       'source_hashes_checked':True,'production_ideal_components':0,
                       'calibration_core_matches':True,'native_errors':0,'native_warnings':0}}
    selected['code_sha256']={q:W['digest'](ROOT/q) for q in (
        'geometry/specifications.py','geometry/taper/three_tapers.py','geometry/taper/three_tapers_v2.py',
        'geometry/taper/three_tapers_v3.py',
        'geometry/taper/measure_compact_inductance.py','geometry/unit_cell.py','geometry/meander.py',
        'geometry/prepare_sonnet.py','geometry/taper/sonnet.py','geometry/taper/solver_guard.py')}
    (ROOT/'geometry/taper/selected_three_tapers.json').write_text(json.dumps(selected,indent=2)+'\n')
    (destination/'selection.json').write_text(json.dumps(selected,indent=2)+'\n')
    (ROOT/'geometry/taper/three_tapers_status.json').write_text(json.dumps(selected,indent=2)+'\n')
    (destination/'README.md').write_text(
        '# New T1: accepted coarse design\n\n'
        'Open `response/generated_ncell_2.sonx` for the simulated T1 + regular fixture.\n'
        'Open `s21_fit.png` for Sonnet, fitted circuit, fixed target, and parameter errors.\n'
        f'The nominal targets are {TAPER_1.L0_nH} nH, {TAPER_1.C0_fF} fF C0, and {TAPER_1.Cg_input_fF} fF input Cg.\n'
        f"Measured L = {measured['L0_nH']:.8f} nH. Fixed-target complex RMS = {score['complex_RMS']:.6f}.\n"
        'The fixture GDS is `geometry.gds`; `t1_unit.gds` contains just the modular T1.\n'
        'The isolated L measurement is in `Lcal/`.\n'
        'This verifies T1 in this test fixture. T2, T3 and the complete filter are not yet verified.\n'
        'See `selection.json` for all parameters, extraction, specification and audit details.\n',encoding='utf-8')
    print('FINAL T1',destination,flush=True)
    return selected


if __name__=='__main__':
    measure('w224_pad8_d20',BASE)
