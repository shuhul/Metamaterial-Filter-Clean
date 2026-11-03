"""T2 experiments with verified T1/right-hand metal and exactly two fit variables.

Historical second_experiment fits are preserved separately. They are never
called to fit these trials. No T1, regular-cell, coupling or delay parameter
is released by this module.
"""
from pathlib import Path
from dataclasses import asdict
from hashlib import sha256
import json
import xml.etree.ElementTree as ET
import gdstk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon,Rectangle,Patch
from scipy.optimize import least_squares
from circuit_common import read_snp
from geometry.meander import MeanderParams
from geometry.taper.layout import TaperParams
from geometry.taper.second_layout import SecondParams,build,write,validate,from_dict,require_same_isolation_geometry
from geometry.taper.second_experiment import measure
from geometry.taper.fitting import chain
from geometry.taper.sonnet import configure,simulate
from geometry.prepare_sonnet import main as prepare,points,tag
from geometry.specifications import TAPER_1,TAPER_2,specification

ROOT=Path(__file__).resolve().parents[2]
SOURCE=ROOT/'geometry/taper/selected_design.json'
VERIFIED=json.loads(SOURCE.read_text())
FIXED=VERIFIED['primary_fit']
g=dict(VERIFIED['geometry']);g['meander']=MeanderParams(**g['meander'])
T1=TaperParams(**g)


def require_locked(p,*,allow_t2_junction=False):
    if p.taper1!=T1:raise ValueError('Verified T1 geometry is locked')
    if p.d_t1_regular!=VERIFIED['d_junction'] or p.t_junction_ground!=4:
        raise ValueError('Verified T1-to-regular junction is locked')
    if any((p.w_t1_ground_slot,p.w_upper_plate,p.t_coupling_finger,p.n_coupling_teeth)):
        raise ValueError('No extra fingers, plates or ground slots around T1')
    coupler=VERIFIED['input_coupler']
    if p.taper.t_top!=T1.t_top or (not allow_t2_junction and
        (p.taper.t_right,p.d_t2_t1)!=(coupler['t'],coupler['d_input'])):
        raise ValueError('T2 coupling-leg width, cap top and separation are locked; its length is a T2 control')


def locked_validate(gds,meta,*,allow_t2_junction=False):
    p=from_dict(meta['second']);require_locked(p,allow_t2_junction=allow_t2_junction)
    result=validate(gds,meta)
    if meta['isolation']:return result
    # Prove that every metal polygon at/after T1 is the verified layout,
    # translated as a rigid block. This includes both weak gaps and grounds.
    ref_gds=Path(VERIFIED['artifacts']['fixture_gds'])
    ref_meta=json.loads(ref_gds.with_suffix('.json').read_text())
    old_left=ref_meta['cell_centers'][1]-T1.w/2
    new_left=meta['cell_centers'][2]-T1.w/2;shift=new_left-old_left
    ref=gdstk.read_gds(str(ref_gds)).top_level()[0].get_polygons()
    ref=[poly.copy().translate(shift,0) for poly in ref]
    box=gdstk.rectangle((new_left,-1052),(meta['w_sim'],448))
    expected=gdstk.boolean(ref,[box],'and',precision=.001)
    actual=gdstk.boolean(result[1],[box],'and',precision=.001)
    if gdstk.boolean(actual,expected,'xor',precision=.001):
        raise ValueError('Metal at or after T1 differs from the verified design')
    return result


def require_native_locked(project,meta,*,allowed_regions=None):
    """Also compare imported Sonnet polygons directly to the verified .sonx."""
    old_meta=json.loads(Path(VERIFIED['artifacts']['fixture_gds']).with_suffix('.json').read_text())
    old_left=old_meta['cell_centers'][1]-T1.w/2
    new_left=meta['cell_centers'][2]-T1.w/2;shift=new_left-old_left
    ref=[gdstk.Polygon(points(q)).translate(shift,0) for q in ET.parse(VERIFIED['artifacts']['sonnet_project']).iter(tag('PlanarPolygon'))]
    actual=[gdstk.Polygon(points(q)) for q in ET.parse(project).iter(tag('PlanarPolygon'))]
    box=gdstk.rectangle((new_left,0),(meta['w_sim'],1500))
    expected=gdstk.boolean(ref,[box],'and',precision=.001)
    clipped=gdstk.boolean(actual,[box],'and',precision=.001)
    difference=gdstk.boolean(clipped,expected,'xor',precision=.001)
    if allowed_regions:
        difference=gdstk.boolean(difference,allowed_regions,'not',precision=.001)
    if difference:
        raise ValueError('Imported Sonnet metal at/after T1 differs from verified project')
    return True


def preview(gds,path):
    meta=json.loads(Path(gds).with_suffix('.json').read_text());p=from_dict(meta['second'])
    top,_=build(p,isolation=meta['isolation'])
    fig,ax=plt.subplots(figsize=(14,6))
    for poly in top.get_polygons():ax.add_patch(Polygon(poly.points,facecolor='#388563',edgecolor='none'))
    for ref,color in zip(top.references,('#246cb5','#db8b2b','#388563','#388563')):
        for poly in ref.get_polygons():ax.add_patch(Polygon(poly.points,facecolor=color,edgecolor='none'))
    if not meta['isolation']:
        # Color is only a review aid; the GDS contains one metal layer.
        for x0,y0,x1,y1 in ((0,0,meta['w_sim'],448),(0,-1052,meta['w_sim'],-618),
                             (0,-1052,412,-24),(meta['w_sim']-412,-1052,meta['w_sim'],-24)):
            ax.add_patch(Rectangle((x0,y0),x1-x0,y1-y0,facecolor='#b7bdc3',edgecolor='none'))
        if p.l_ground_rise:
            x=meta['cell_centers'][1]+p.x_meander
            ax.add_patch(Rectangle((x-p.w_ground_riser/2,-618),p.w_ground_riser,p.l_ground_rise,
                                   facecolor='#b7bdc3',edgecolor='none'))
            ax.text(x,-618+p.l_ground_rise/2,'GND',ha='center',va='center',fontsize=8,rotation=90)
        trace_x=meta['cell_centers'][1]+p.x_meander
        ax.annotate(f'{p.taper.meander.t:g} um trace',xy=(trace_x,-200),
                    xytext=(trace_x+32,-240),fontsize=10,
                    arrowprops={'arrowstyle':'->','lw':.8})
    names=['Input IDC','T2','Verified T1','Regular','Output end']
    if meta['isolation']:names=['Input feed','T2','Calibration feed']
    for x,label in zip(meta['cell_centers'],names):ax.text(x,25,label,ha='center',fontsize=12)
    ax.set(xlim=(350,meta['w_sim']-350),ylim=(-660,65),xlabel='x (um)',ylabel='y (um)',
           title='T2 Candidate — Verified T1 and Right-Hand Geometry Locked')
    ax.set_aspect('equal');fig.tight_layout();fig.savefig(path,dpi=200);plt.close(fig)


def response(f,x,l2):
    """x contains exactly C0_t2 and Cg_t2, in fF. Everything else is fixed."""
    if len(x)!=2:raise ValueError('Exactly two free capacitances are allowed')
    s=chain(2*np.pi*f,[l2,FIXED['L_taper_nH'],FIXED['L_regular_nH']],
            [x[0],FIXED['C_taper_fF'],FIXED['C_regular_fF']],
            [x[1],FIXED['Cg_input_fF'],FIXED['Cg_junction_fF'],FIXED['Cg_output_fF']])
    left=FIXED['fixture_delay_left_ps']*1e-12;right=FIXED['fixture_delay_right_ps']*1e-12
    s[:,0,0]*=np.exp(-4j*np.pi*f*left);s[:,1,1]*=np.exp(-4j*np.pi*f*right)
    s[:,1,0]*=np.exp(-2j*np.pi*f*(left+right));s[:,0,1]=s[:,1,0]
    return s


def fit(snp,output,l2):
    output=Path(output);output.mkdir(parents=True,exist_ok=True);net=read_snp(snp)
    def residual(x):
        d=response(net.f,x,l2)[:,1,0]-net.s[:,1,0]
        return np.r_[d.real,d.imag]
    # Narrow resonances can trap local solvers even with two variables.
    # Screen a coarse two-dimensional grid using the same objective/data,
    # then refine the best seeds. This adds no EM simulations or free values.
    candidates=[np.array([c0,cg]) for c0 in np.linspace(10,300,21) for cg in np.linspace(100,1000,21)]
    ranked=sorted(candidates,key=lambda x:float(residual(x)@residual(x)))
    starts=[[100,459],[100,300],[100,700]]
    for x in ranked:
        if all(np.linalg.norm((x-np.array(old))/[290,900])>.08 for old in starts):starts.append(x)
        if len(starts)>=11:break
    fits=[least_squares(residual,x,bounds=([10,100],[300,1000]),x_scale='jac',max_nfev=250,ftol=1e-11)
          for x in starts]
    r=min(fits,key=lambda q:q.fun@q.fun);model=response(net.f,r.x,l2)
    errors={'L0_t2':100*(l2/TAPER_2.L0_nH-1),
            'C0_t2':100*(r.x[0]/TAPER_2.C0_fF-1),'Cg_t2':100*(r.x[1]/TAPER_2.Cg_input_fF-1)}
    result={'free_parameters':['C0_t2_fF','Cg_t2_fF'],'x':r.x.tolist(),'L_t2_nH':l2,
            'C0_t2_fF':float(r.x[0]),'Cg_t2_fF':float(r.x[1]),'fixed_model':FIXED,
            'fixed_reference_manifest':str(SOURCE),'fixed_reference_sha256':sha256(SOURCE.read_bytes()).hexdigest(),
            'objective':'Complex S21 real and imaginary parts; exactly two fitted parameters; fixed reference delays',
            'success':bool(r.success),'active_bounds':r.active_mask.tolist(),'errors_percent':errors,
            'passes_targets':bool(r.success and not any(r.active_mask) and all(abs(v)<=2 for v in errors.values())),
            's21_complex_rmse':float(np.sqrt(np.mean(abs(model[:,1,0]-net.s[:,1,0])**2))),
            's21_max_complex_error':float(np.max(abs(model[:,1,0]-net.s[:,1,0]))),
            'fit_range_GHz':[float(net.f[0]/1e9),float(net.f[-1]/1e9)]}
    result['optimizer_grid_points']=len(candidates);result['optimizer_starts']=len(starts)
    (output/'fit_locked.json').write_text(json.dumps(result,indent=2)+'\n')
    plot(snp,result,output/'s21_fit_linear.png')
    print('LOCKED FIT',json.dumps(result),flush=True)
    return result


def plot(snp,result,path):
    net=read_snp(snp);pred=response(net.f,result['x'],result['L_t2_nH'])
    fig,ax=plt.subplots(figsize=(10,6))
    ax.plot(net.f/1e9,abs(net.s[:,1,0]),'r-',lw=2.2,label='Sonnet simulation')
    ax.plot(net.f/1e9,abs(pred[:,1,0]),'k--',lw=2.2,label='Circuit fit')
    labels=[(r'$L_{t2}$',result['L_t2_nH'],'nH','L0_t2'),
            (r'$C_{0,t2}$',result['C0_t2_fF'],'fF','C0_t2'),
            (r'$C_{g,t2}$',result['Cg_t2_fF'],'fF','Cg_t2')]
    text='\n'.join(f'{name} = {value:.4f} {unit} ({result["errors_percent"][key]:+.2f}%)' for name,value,unit,key in labels)
    cg1=FIXED['Cg_input_fF']
    text+=f'\n$C_{{g,t1}}$ = {cg1:.3f} fF ({100*(cg1/TAPER_1.Cg_input_fF-1):+.2f}%, fixed)'
    text+='\n\nT1 and regular parameters fixed\nOnly two capacitances fitted'
    ax.text(.035,.95,text,transform=ax.transAxes,va='top',fontsize=13,
            bbox=dict(boxstyle='round,pad=.5',facecolor='white',edgecolor='gray',alpha=.95))
    peak=max(abs(net.s[:,1,0]).max(),abs(pred[:,1,0]).max());ax.set_ylim(0,peak*1.8)
    ax.set_xlabel('Frequency (GHz)',fontsize=15);ax.set_ylabel(r'$|S_{21}|$',fontsize=16)
    title='Second Taper Cell: Transmission and Circuit Fit'
    if not result.get('passes_targets',False):title+='\nCandidate — Targets Not Met'
    ax.set_title(title,fontsize=17,pad=12)
    ax.tick_params(labelsize=13);ax.legend(loc='upper right',fontsize=12);ax.grid(linestyle=':',alpha=.6)
    fig.tight_layout();fig.savefig(path,dpi=220);plt.close(fig)


def prepare_case(name,p,*,allow_t2_junction=False):
    require_locked(p,allow_t2_junction=allow_t2_junction)
    gds=ROOT/'geometry/taper/output'/f'{name}.gds'
    if gds.exists():
        meta=json.loads(gds.with_suffix('.json').read_text())
        if from_dict(meta['second'])!=p:raise ValueError('Changed geometry needs a new case name')
    else:
        gds=write(p,gds);meta=json.loads(gds.with_suffix('.json').read_text())
    locked_validate(gds,meta,allow_t2_junction=allow_t2_junction)
    out=ROOT/'output_data/taper_exploration'/name;out.mkdir(parents=True,exist_ok=True)
    preview(gds,out/'geometry_review.png')
    report={'case':name,'t1_and_right_hand_metal_xor_empty':True,'verified_T1':asdict(T1),
            'parameters':asdict(p),'free_fit_parameters':['C0_t2_fF','Cg_t2_fF'],
            'fixed_reference_sha256':sha256(SOURCE.read_bytes()).hexdigest(),'specification':specification()}
    (out/'geometry_check.json').write_text(json.dumps(report,indent=2)+'\n')
    return gds,out


def run(name,p,*,inductance_reference=None,sweep=None,sweep_target=500):
    cached=ROOT/'output_data/taper_exploration'/name/'trial_locked.json'
    if cached.exists():
        record=json.loads(cached.read_text())
        if from_dict(record['parameters'])!=p:raise ValueError('Changed geometry needs a new case name')
        if sweep is not None and not np.allclose(record['fit']['fit_range_GHz'],sweep,rtol=0,atol=1e-9):
            raise ValueError('Changed simulation band needs a new case name')
        for path,h in record['hashes'].items():
            if sha256(Path(path).read_bytes()).hexdigest()!=h:raise ValueError('Completed trial artifact changed')
        locked_validate(record['gds'],json.loads(Path(record['gds']).with_suffix('.json').read_text()))
        print('Reusing completed locked trial; no EM launched:',name,flush=True)
        return record
    gds,out=prepare_case(name,p)
    if inductance_reference is None:
        # Match the 5 GHz operating band rather than extracting a slope at
        # 8-11 GHz. Separate case names preserve the earlier calibration data.
        inductance_reference,lr=measure(name+'_designband',p,c_add_sweep_pf=(.8,1.,1.2,1.4),
                                       calibration_band=(4.,6.5),calibration_target=2500)
    else:
        inductance_reference=Path(inductance_reference)
        require_same_isolation_geometry(p,ROOT/'geometry/taper/output'/f'{inductance_reference.parent.name}.gds')
        lr=json.loads(inductance_reference.read_text())
    if abs(lr['L0_nH']/TAPER_2.L0_nH-1)>.02:
        print('L outside target; skipping full mixed simulation:',lr['L0_nH'],flush=True)
        screening={'status':'inductance_outside_target','L_t2_nH':lr['L0_nH'],
                   'inductance_reference':str(inductance_reference),'parameters':asdict(p)}
        (out/'screening.json').write_text(json.dumps(screening,indent=2)+'\n')
        return screening
    dest=ROOT/'sonnet_files/taper_exploration'/name
    project=dest/'generated_ncell_3.sonx'
    prepare(3,gds_path=gds,output_dir=dest,geometry_validator=locked_validate)
    require_native_locked(project,json.loads(gds.with_suffix('.json').read_text()))
    band=(4.5,5.5) if sweep is None else tuple(sweep)
    configure(project,*band,sweep_target)
    snp=simulate(project,out)
    result=fit(snp,out,lr['L0_nH'])
    record={'case':name,'status':'target_screening_pass' if result['passes_targets'] else 'targets_not_met',
            'gds':str(gds),'project':str(project),'snp':str(snp),'inductance_reference':str(inductance_reference),
            'parameters':asdict(p),'fit':result,'T1_geometry_unchanged':True,
            'simulation_band_GHz':list(band),'simulation_target':sweep_target,
            'hashes':{str(path):sha256(Path(path).read_bytes()).hexdigest() for path in (gds,project,snp,inductance_reference)}}
    (out/'trial_locked.json').write_text(json.dumps(record,indent=2)+'\n')
    return record
