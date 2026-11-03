"""T2 for the current specification; accepted T1 cell components are locked.

Edit the grouped parameters and flags below. Every changed layout needs a new
CASE name. Declared assembly gaps/centered ground widths are tunable separately.
All geometry uses microns, layer 1/datatype 0 and the 2 um Sonnet grid.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from dataclasses import dataclass,asdict,replace
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
from scipy.optimize import least_squares,minimize_scalar
from geometry.unit_cell import make_unit_cell
from geometry.specifications import TAPER_1,TAPER_2,UNIT_CELL,specification,S21_COMPLEX_RMS_CUTOFF
from geometry.taper import three_tapers as shared
from geometry.taper import three_tapers_v3 as verified
from geometry.taper.measure_compact_inductance import calibrate,W,AP
from geometry.taper.sonnet import add_capacitor,configure,simulate
from geometry.taper.progressive_fitting import transmission,fit as fit_circuit
from geometry.prepare_sonnet import main as import_native,tag,points,write_xml
from circuit_common import read_snp

ROOT=shared.ROOT
CAMPAIGN='three_tapers_v3_T2_2026-10-26'
REFERENCE_FILE=ROOT/'geometry/taper/taper2_T1_reference.json'
REFERENCE=json.loads(REFERENCE_FILE.read_text())
T1=shared.from_dict(REFERENCE['geometry']['unit'])
REGULAR=shared.REGULAR
REGULAR_MODEL=REFERENCE['fixed_regular_extraction_model']
T1_MODEL=REFERENCE['heuristic_fit']


@dataclass(frozen=True)
class InputCoupler:
    d: int=2                 # Edge-to-edge rail separation.
    n: int=0                 # Zero = simple parallel legs; otherwise IDC teeth.
    t: int=4                 # Finger thickness.
    d_finger: int=4          # Clear vertical distance between adjacent fingers.
    d_tip: int=2             # Finger tip clearance to the opposite rail.
    d_t1: int=56              # T2-to-T1 cap-edge separation (no ground strip).
    d_regular: int=66         # T1-to-regular cap-edge separation.
    t_regular_ground: int=6   # Centered assembly ground strip, not a cell part.
    t_input: int|None=None    # Optional feed-electrode width for isolated junction tests.
    w_input_stub: int|None=None  # Horizontal feed span before the input rail; legacy default = cap.w.


BASE=replace(T1,cap=replace(T1.cap,t_side=64),
             meander=replace(T1.meander,w=152),pad=replace(T1.pad,l=260))
COUPLER=InputCoupler(d=22,n=70,d_t1=126,d_regular=62,t_regular_ground=6)
CASE='thick64_p260_w152_i22n70_j126_r62'
RUN_INDUCTANCE=True
RUN_RESPONSE=True
FINALIZE_ON_PASS=True


def require_reference():
    if not REFERENCE['T1_accepted'] or not REFERENCE['fixed_target_comparison']['pass']:
        raise ValueError('A verified T1 reference is required')
    if REFERENCE['specification']!=specification():raise ValueError('T1 target specification changed')
    project=Path(REFERENCE['project'])
    if W['digest'](project.with_suffix('.snp'))!=REFERENCE['fixed_target_comparison']['source_sha256']:
        raise ValueError('Verified T1 S-parameters changed')
    run=ROOT/'output_data/taper_exploration'/REFERENCE['specification']['version']/REFERENCE['case']/'response/run_status.json'
    if W['digest'](project)!=json.loads(run.read_text())['project_sha256']:
        raise ValueError('Verified T1 native project changed')


def build(p,*,isolation=False,coupler=InputCoupler(),input_probe=False):
    require_reference()
    if p.d_top+p.cap.t_top+p.meander.l!=618:raise ValueError('Inductor must meet flat bottom ground')
    if p.cap.l!=T1.cap.l:
        raise ValueError('Keep facing cap legs at the same vertical length')
    if p.cap.t_side<2 or p.cap.t_side%2:raise ValueError('Cap legs must be even and >=2 um')
    if input_probe and not isolation:raise ValueError('Input probe requires the isolated fixture')
    if isolation and not input_probe:
        top,meta=verified.build(p,isolation=True,d_input=56,d_junction=66,t_junction_ground=6)
        meta['layout_kind']='T2_weak_probe_isolation_v3'
    else:
        if coupler.d<2 or coupler.d%2:raise ValueError('Input gap must be even and >=2 um')
        if p.cap.w%4:raise ValueError('T2 cap width must be divisible by 4')
        units=[p,p,p] if input_probe else [p,p,T1,REGULAR,REGULAR]
        if coupler.d_t1<2 or coupler.d_t1%2:raise ValueError('T2-to-T1 gap must be even and >=2 um')
        if (coupler.t_regular_ground<2 or coupler.t_regular_ground%2 or
            coupler.d_regular<coupler.t_regular_ground+4 or
            (coupler.d_regular-coupler.t_regular_ground)%4):
            raise ValueError('Centered regular-junction strip must lie on the 2 um grid')
        gaps=[coupler.d,98] if input_probe else [coupler.d,coupler.d_t1,
              coupler.d_regular,56]
        spans=[u.cap.w for u in units]
        if coupler.w_input_stub is not None:
            if coupler.w_input_stub<2 or coupler.w_input_stub%2:
                raise ValueError('Input feed span must be positive and even')
            spans[0]=coupler.w_input_stub
        lefts=[426]
        for span,gap in zip(spans,gaps):lefts.append(lefts[-1]+span+gap)
        centers=[x+span/2 for x,span in zip(lefts,spans)]
        width=lefts[-1]+units[-1].cap.w+426
        top=gdstk.Cell('TOP')
        for i in ((1,) if input_probe else (1,2,3,4)):
            top.add(gdstk.Reference(make_unit_cell(units[i],name=f'STAGE_{i}'),origin=(centers[i],0)))
        edge=lefts[0]+spans[0]
        rail_t=p.cap.t_side if coupler.t_input is None else coupler.t_input
        if rail_t<2 or rail_t%2 or rail_t>p.cap.w:raise ValueError('Input rail must be even and fit its scaffold')
        if spans[0]<rail_t:raise ValueError('Input rail must fit within the feed span')
        parts=[gdstk.rectangle((edge-rail_t,-p.d_top-p.cap.l),(edge,-p.d_top)),
               gdstk.rectangle((0,-18),(edge,-8)),
               gdstk.rectangle((0,0),(width,448)),
               gdstk.rectangle((0,-1052),(width,-618)),
               gdstk.rectangle((0,-1052),(412,-24)),
               gdstk.rectangle((width-412,-1052),(width,-24))]
        if input_probe:
            right=lefts[2]
            parts.extend([gdstk.rectangle((right,-p.d_top-p.cap.l),(right+p.cap.t_side,-p.d_top)),
                          gdstk.rectangle((right,-18),(width,-8))])
        else:parts.append(gdstk.rectangle((centers[-1],-18),(width,-8)))
        # No ground strip in either input->T2 or T2->T1 gap.
        for i,thickness in (((1,6),) if input_probe else ((2,coupler.t_regular_ground),(3,4))):
            mid=lefts[i]+units[i].cap.w+gaps[i]/2
            parts.append(gdstk.rectangle((mid-thickness/2,-618),(mid+thickness/2,0)))
        top.add(*gdstk.boolean(parts,[],'or',precision=.001,layer=1))
        if coupler.n:
            if (not isinstance(coupler.n,int) or coupler.n<2 or
                any(v<2 or v%2 for v in (coupler.t,coupler.d_finger,coupler.d_tip)) or
                coupler.d<=2*coupler.d_tip):raise ValueError('Invalid on-grid IDC')
            upper=-p.d_top-p.cap.t_top-2
            bottom=upper-(coupler.n-1)*(coupler.t+coupler.d_finger)-coupler.t
            if bottom<-p.d_top-p.cap.l+2:raise ValueError('IDC teeth extend below the equal-length rails')
            for i in range(coupler.n):
                y=upper-i*(coupler.t+coupler.d_finger)
                if i%2:x0,x1=edge+coupler.d_tip,edge+coupler.d
                else:x0,x1=edge,edge+coupler.d-coupler.d_tip
                top.add(gdstk.rectangle((x0,y-coupler.t),(x1,y),layer=1))
        coupling_meta=asdict(coupler)
        if coupler.w_input_stub is None:coupling_meta.pop('w_input_stub')
        if coupler.t_input is None:coupling_meta.pop('t_input')
        if coupler.d_t1==56:coupling_meta.pop('d_t1') # Preserve existing case identities.
        if coupler.d_regular==66:coupling_meta.pop('d_regular')
        if coupler.t_regular_ground==6:coupling_meta.pop('t_regular_ground')
        meta={'layout_kind':'T2_input_coupler_probe_v3' if input_probe else 'T2_T1_regular_v3',
              'unit':asdict(p),'isolation':input_probe,
              'taper1':asdict(T1),'regular':asdict(REGULAR),'coupler':coupling_meta,
              'w_sim':width,'cell_centers':centers,'cell_lefts':lefts,'gaps':gaps,
              'expected_metal_components':3 if input_probe else 2,
              'parameters':{'n_cells':1 if input_probe else 3,'l_sim':1500,'l_top_ground':448,'w_side':412,'t_feed':10,
                            'unit':{'d_top':p.d_top,'cap':{'t_top':p.cap.t_top}}}}
        if input_probe:meta['input_probe']=True
    meta['target_stage']='taper_2'
    meta['T1_reference_sha256']=W['digest'](REFERENCE_FILE)
    return top,meta


def validate(path,meta):
    p=shared.from_dict(meta['unit'])
    top,expected=build(p,isolation=meta['isolation'],coupler=InputCoupler(**meta.get('coupler',{})),
                       input_probe=meta.get('input_probe',False))
    if meta!=expected:raise ValueError('Geometry metadata changed; use a new CASE')
    lib=gdstk.read_gds(str(path));tops=lib.top_level()
    if len(tops)!=1 or tops[0].name!='TOP' or not np.isclose(lib.unit,1e-6):raise ValueError('Invalid GDS')
    polys=tops[0].get_polygons()
    if any(q.layer!=1 or q.datatype!=0 for q in polys):raise ValueError('Wrong metal layer')
    if gdstk.boolean(polys,top.get_polygons(),'xor',precision=.001):raise ValueError('GDS differs from parameters')
    if len(gdstk.boolean(polys,[],'or',precision=.001))!=meta['expected_metal_components']:
        raise ValueError('Wrong feed/ground connectivity')
    if not meta['isolation']:require_locked(polys,meta)
    return SimpleNamespace(n_cells=meta['parameters']['n_cells'],l_sim=1500,l_top_ground=448),polys


def require_locked(polys,meta,*,native=False):
    old=REFERENCE['geometry'];old_left=old['cell_centers'][1]-T1.cap.w/2
    new_left=meta['cell_lefts'][2];shift=new_left-old_left
    if shift%2:raise ValueError('Translation would change native snapping')
    if native:
        ref=[gdstk.Polygon(points(q)) for q in ET.parse(REFERENCE['project']).iter(tag('PlanarPolygon'))]
        y0,y1=0,1500
    else:
        ref=gdstk.read_gds(REFERENCE['gds']).top_level()[0].get_polygons();y0,y1=-1052,448
    ref=[q.copy().translate(shift,0) for q in ref]
    coupling=InputCoupler(**meta['coupler'])
    if coupling.d_regular!=66 or coupling.t_regular_ground!=6:
        # The cells are immutable; only this explicitly named assembly
        # junction may change. Translate the regular/output group together.
        edge=new_left+T1.cap.w;old_regular=edge+66
        delta=coupling.d_regular-66
        left=gdstk.rectangle((-1e6,y0),(old_regular,y1))
        right=gdstk.rectangle((old_regular,y0),(1e6,y1))
        a=gdstk.boolean(ref,[left],'and',precision=.001)
        b=[q.translate(delta,0) for q in gdstk.boolean(ref,[right],'and',precision=.001)]
        lo,hi=(448,1066) if native else (-618,0)
        a=gdstk.boolean(a,[gdstk.rectangle((edge,lo),(old_regular,hi))],'not',precision=.001)
        mid=edge+coupling.d_regular/2;half=coupling.t_regular_ground/2
        ref=a+b+[gdstk.rectangle((mid-half,lo),(mid+half,hi))]
        # Keep the two planes solid for either sign of the group translation.
        bands=((0,448),(1066,1500)) if native else ((0,448),(-1052,-618))
        ref += [gdstk.rectangle((new_left,low),(meta['w_sim'],high)) for low,high in bands]
    box=gdstk.rectangle((new_left,y0),(meta['w_sim'],y1))
    a=gdstk.boolean(polys,[box],'and',precision=.001)
    b=gdstk.boolean(ref,[box],'and',precision=.001)
    if gdstk.boolean(a,b,'xor',precision=.001):raise ValueError('Cell metal or undeclared assembly geometry changed')


def prepare(name,p,*,isolation=False,coupler=InputCoupler(),input_probe=False):
    stage='input_probe' if input_probe else ('Lcal' if isolation else 'response')
    gd=ROOT/'geometry/taper/output'/CAMPAIGN/name/stage
    native=ROOT/'sonnet_files/taper_exploration'/CAMPAIGN/name/stage
    out=ROOT/'output_data/taper_exploration'/CAMPAIGN/name/stage
    for folder in (gd,native,out):folder.mkdir(parents=True,exist_ok=True)
    top,meta=build(p,isolation=isolation,coupler=coupler,input_probe=input_probe)
    path=gd/'geometry.gds'
    if path.exists():
        if json.loads(path.with_suffix('.json').read_text())!=meta:raise ValueError('Case geometry changed')
    else:
        lib=gdstk.Library(unit=1e-6,precision=1e-9);lib.add(top,*top.dependencies(True));lib.write_gds(str(path))
        path.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
    validate(path,meta)
    fig,ax=plt.subplots(figsize=(17,6))
    for poly in top.get_polygons():ax.add_patch(Polygon(poly.points,fc='#297b54',ec='none'))
    purpose='Isolated Coupling Probe' if input_probe else ('Independent L Calibration' if isolation else 'T2 + Verified T1 + Regular')
    ax.set(xlim=(360,meta['w_sim']-360),ylim=(-650,45),aspect='equal',xlabel='x (um)',ylabel='y (um)',
           title=f'{name} - {purpose}')
    fig.tight_layout();fig.savefig(out/'geometry.png',dpi=190);plt.close(fig)
    project=native/f"generated_ncell_{meta['parameters']['n_cells']}.sonx"
    if not project.exists():import_native(meta['parameters']['n_cells'],gds_path=path,output_dir=native,geometry_validator=validate)
    if not isolation:
        require_locked([gdstk.Polygon(points(q)) for q in ET.parse(project).iter(tag('PlanarPolygon'))],meta,native=True)
    (out/'specification.json').write_text(json.dumps(specification(),indent=2)+'\n')
    return project,out,meta


def prepare_probe(project,meta,frequencies=(1,2,2.5,3,4)):
    """Turn the capacitor-load pins into an internal analysis port pair."""
    tree=ET.parse(project)
    if not list(tree.iter(tag('Component'))):
        add_capacitor(project,meta);tree=ET.parse(project)
        comp=next(tree.iter(tag('Component')));comp.set('Type','PORTS');comp.attrib.pop('Units',None)
        for q in list(comp):
            if q.tag==tag('Value'):comp.remove(q)
        sweeps=tree.find(tag('Sweeps'))
        for q in list(sweeps):sweeps.remove(q)
        fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
        for f in frequencies:
            ET.SubElement(fs,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(f),Stop=str(f),Step='1')
        for q in tree.findall(f"{tag('OutputFiles')}//{tag('File')}"):
            if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','FALSE')
        write_xml(tree,project)
    return project


def measure(name,p):
    # Weak calibration electrodes isolate intrinsic L from the strong input IDC.
    project,out,meta=prepare(name,p,isolation=True)
    prepare_probe(project,meta)
    return calibrate(project,out,name,auto_loads=True,target_L_nH=TAPER_2.L0_nH,
                     held_out_GHz=2.5,frequency_range_GHz=(1,4),desired_GHz=[3.8,3.4,3,2.6])


def reuse_calibration(source_name,name,p):
    """Reuse measured L only for the exact same body and weak-probe geometry."""
    if source_name==name:return
    _,expected=build(p,isolation=True)
    source_gd=ROOT/'geometry/taper/output'/CAMPAIGN/source_name/'Lcal'
    if json.loads((source_gd/'geometry.json').read_text())!=expected:raise ValueError('Calibration geometry differs')
    src=ROOT/'sonnet_files/taper_exploration'/CAMPAIGN/source_name/'Lcal'
    out=ROOT/'output_data/taper_exploration'/CAMPAIGN/source_name/'Lcal'
    marker=json.loads((out/'operator_completed.json').read_text())
    project=src/'generated_ncell_1.sonx'
    if W['digest'](project)!=marker['project_sha256'] or W['digest'](project.with_suffix('.snp'))!=marker['snp_sha256']:
        raise ValueError('Calibration source changed')
    for source,dest in ((source_gd,ROOT/'geometry/taper/output'/CAMPAIGN/name/'Lcal'),
                        (src,ROOT/'sonnet_files/taper_exploration'/CAMPAIGN/name/'Lcal'),
                        (out,ROOT/'output_data/taper_exploration'/CAMPAIGN/name/'Lcal')):
        if dest.exists():raise ValueError('Calibration destination already exists')
        shutil.copytree(source,dest,ignore=shutil.ignore_patterns('sondata'))
    result=ROOT/'output_data/taper_exploration'/CAMPAIGN/name/'Lcal/results.json'
    saved=json.loads(result.read_text());saved['reused_from_case']=source_name
    result.write_text(json.dumps(saved,indent=2)+'\n')


def measure_input(name,p,*,coupler=InputCoupler(),frequencies=(1,2,4),fixed_L_nH=None):
    project,out,meta=prepare(name,p,isolation=True,coupler=coupler,input_probe=True)
    # Three points suffice for the two-term K/omega + omega*C screening model.
    # Independent L calibration retains all five points and its held-out check.
    prepare_probe(project,meta,frequencies=frequencies)
    marker=out/'operator_completed.json'
    if marker.exists():
        old=json.loads(marker.read_text())
        if old['project_sha256']!=W['digest'](project) or old['snp_sha256']!=W['digest'](project.with_suffix('.snp')):
            raise ValueError('Input-probe source changed')
        net=W['read_network'](project.with_suffix('.snp'),4)
    else:
        net=W['em'](project,4)
        marker.write_text(json.dumps({'project_sha256':W['digest'](project),
                                     'snp_sha256':W['digest'](project.with_suffix('.snp'))},indent=2)+'\n')
    if not np.allclose(net.f,np.array(frequencies)*1e9,rtol=0,atol=1):
        raise ValueError('Probe frequencies changed; use a fresh case name')
    s,z=AP['reduce_pairs'](net.s)
    k,c,error=AP['extract_kc'](net.f,AP['s_to_y'](s,z))
    result={'scope':'Isolated terminal-admittance proxies; not full S21 acceptance',
            'coupler':asdict(coupler),'ground_C_proxy_fF':float(c[2].sum()),
            'input_Cg_proxy_fF':float(-c[2,0]),'dummy_Cg_proxy_fF':float(-c[2,1]),
            'low_frequency_inverse_L_per_nH':float(k[2,2]),
            'low_frequency_L_nH':float(1/k[2,2]),'KC_max_admittance_residual_S':float(error)}
    result['frequency_range_GHz']=[float(net.f[0]/1e9),float(net.f[-1]/1e9)]
    if min(frequencies)>=4:
        result['scope']='Operating-band terminal admittance; local K/C coefficients are not physical component acceptance'
        result['KC_local_fit_coefficients']={
            'ground_C_fF':result.pop('ground_C_proxy_fF'),
            'input_Cg_fF':result.pop('input_Cg_proxy_fF'),
            'output_Cg_fF':result.pop('dummy_Cg_proxy_fF'),
            'inverse_L_per_nH':result.pop('low_frequency_inverse_L_per_nH'),
            'L_nH':result.pop('low_frequency_L_nH')}
    if fixed_L_nH is not None:
        y=AP['s_to_y'](s,z);omega=2*np.pi*net.f
        result['fixed_L_nH']=fixed_L_nH
        result['terminal_admittance_at_each_frequency']=[
            {'frequency_GHz':float(f/1e9),
             'input_Cg_fF':float(-yy[2,0].imag/w*1e15),
             'output_Cg_fF':float(-yy[2,1].imag/w*1e15),
             'ground_C_at_fixed_L_fF':float((yy[2].sum().imag+1/(w*fixed_L_nH*1e-9))/w*1e15)}
            for f,w,yy in zip(net.f,omega,y)]
        result['fixed_L_scope']='Nodal ground susceptance with separately measured L held fixed; not an S21 fitted capacitance'
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result


def response(name,p,l,*,coupler=InputCoupler(),prepared_case=None,fixed_T1_L_nH=None):
    """Keep L fixed; fit taper C values as guides, score fixed nominal targets."""
    project,out,meta=prepare(name,p,coupler=coupler) if prepared_case is None else prepared_case
    marker=out/'run_status.json'
    if marker.exists() and json.loads(marker.read_text()).get('status')=='complete':
        old=json.loads(marker.read_text())
        if old['project_sha256']!=W['digest'](project) or old['snp_sha256']!=W['digest'](project.with_suffix('.snp')):
            raise ValueError('Completed simulation source changed')
        snp=project.with_suffix('.snp')
    else:
        configure(project,4.65,5.4,7500);snp=simulate(project,out)
    ls=[l,REFERENCE['isolated_L']['L0_nH'] if fixed_T1_L_nH is None else fixed_T1_L_nH,
        REGULAR_MODEL['L0_nH'],REGULAR_MODEL['L0_nH']]
    cs=[TAPER_2.C0_fF,T1_MODEL['C_fF'][0],REGULAR_MODEL['C0_fF'],REGULAR_MODEL['C0_fF']]
    gs=[TAPER_2.Cg_input_fF,T1_MODEL['G_fF'][0],T1_MODEL['G_fF'][1],REGULAR_MODEL['output_Cg_fF']]
    free=[('C',0),('G',0),('C',1),('G',1),('G',2)]
    fitted=fit_circuit(snp,out/'extraction',ls,cs,gs,free,title='Second Taper Cell: Transmission')
    net=read_snp(snp);f=net.f;obs=net.s[:,1,0]
    if len(f)!=7501 or not np.allclose([f[0],f[-1]],[4.65e9,5.4e9],rtol=0,atol=1):
        raise ValueError('Acceptance frequency range/grid changed')
    model=transmission(f,ls,fitted['C_fF'],fitted['G_fF'],fitted['delay_ps'])
    # Local conditioning reports whether extracted C values trade off; it
    # never changes the fixed-target curve or silently pins a capacitance.
    columns=[]
    for kind,index in free:
        cp=list(fitted['C_fF']);gp=list(fitted['G_fF'])
        values=cp if kind=='C' else gp
        step=max(values[index],1)*1e-4;values[index]+=step
        delta=(transmission(f,ls,cp,gp,fitted['delay_ps'])-model)/step*max(values[index],1)
        columns.append(np.r_[delta.real,delta.imag])
    delta=-2j*np.pi*f*1e-12*model*20
    columns.append(np.r_[delta.real,delta.imag]);jac=np.array(columns).T
    singular=np.linalg.svd(jac,compute_uv=False)
    cov=np.linalg.pinv(jac.T@jac);den=np.sqrt(np.outer(np.diag(cov),np.diag(cov)))
    diagnostic={'parameters':[f'{kind}[{index}]' for kind,index in free]+['delay_ps'],
                'relative_jacobian_condition':float(singular[0]/singular[-1]),
                'local_parameter_correlation':(cov/den).tolist(),
                'scope':'Local fit sensitivity; capacitances are geometry heuristics, not independently accepted component values'}
    (out/'extraction/conditioning.json').write_text(json.dumps(diagnostic,indent=2)+'\n')
    nominal=transmission(f,[TAPER_2.L0_nH,TAPER_1.L0_nH,UNIT_CELL.L0_nH,UNIT_CELL.L0_nH],
                         [TAPER_2.C0_fF,TAPER_1.C0_fF,UNIT_CELL.C0_fF,UNIT_CELL.C0_fF],
                         [TAPER_2.Cg_input_fF,TAPER_1.Cg_input_fF,UNIT_CELL.Cg_input_fF,UNIT_CELL.Cg_input_fF],0)
    objective=lambda tau:float(np.mean(abs(nominal*np.exp(-2j*np.pi*f*tau*1e-12)-obs)**2))
    grid=np.linspace(-100,100,801);i=int(np.argmin([objective(t) for t in grid]))
    opt=minimize_scalar(objective,bounds=(max(-100,grid[i]-.5),min(100,grid[i]+.5)),method='bounded')
    target=nominal*np.exp(-2j*np.pi*f*opt.x*1e-12)
    score={'complex_RMS':float(np.sqrt(opt.fun)),'pass':bool(opt.success and np.sqrt(opt.fun)<=S21_COMPLEX_RMS_CUTOFF),
           'cutoff':S21_COMPLEX_RMS_CUTOFF,'fixture_delay_ps':float(opt.x),'all_electrical_targets_fixed':True,
           'specification':specification(),'source_sha256':W['digest'](snp),
           'T1_reference_sha256':W['digest'](REFERENCE_FILE),'frequency_range_GHz':[4.65,5.4],'samples':len(f)}
    (out/'target_comparison.json').write_text(json.dumps(score,indent=2)+'\n')
    np.savez(out/'comparison.npz',f=f,observed=obs,fit=model,target=target)
    values=[('L T2',l,TAPER_2.L0_nH,'nH'),('C0 T2',fitted['C_fF'][0],TAPER_2.C0_fF,'fF'),
            ('Cg T2',fitted['G_fF'][0],TAPER_2.Cg_input_fF,'fF'),
            ('C0 T1',fitted['C_fF'][1],TAPER_1.C0_fF,'fF'),
            ('Cg T1',fitted['G_fF'][1],TAPER_1.Cg_input_fF,'fF'),
            ('Cg unit',fitted['G_fF'][2],UNIT_CELL.Cg_input_fF,'fF')]
    text='\n'.join(f'{label} = {v:.4g} {unit} ({100*(v/t-1):+.2f}%)' for label,v,t,unit in values)
    fig,ax=plt.subplots(figsize=(14,8))
    ax.plot(f/1e9,abs(obs),lw=2.8,label='Sonnet EM',color='#2768aa')
    ax.plot(f/1e9,abs(model),'k--',lw=2.2,label='Complex circuit fit')
    ax.plot(f/1e9,abs(target),':',lw=2.6,label='Specified target circuit',color='#db6d28')
    ax.text(.025,.97,text+f"\nTarget RMS = {score['complex_RMS']:.5f} (limit {S21_COMPLEX_RMS_CUTOFF:g})\nFit RMS = {fitted['complex_RMS']:.5f}",
            transform=ax.transAxes,va='top',fontsize=16,bbox=dict(boxstyle='round,pad=.5',fc='white',ec='gray'))
    ax.set(xlabel='Frequency (GHz)',ylabel=r'$|S_{21}|$',title='Second Taper Cell: Transmission',ylim=(0,1.5),xlim=(4.65,5.4))
    ax.xaxis.label.set_size(20);ax.yaxis.label.set_size(23);ax.title.set_size(23);ax.tick_params(labelsize=16)
    ax.grid(alpha=.2);ax.legend(fontsize=15,loc='upper right');fig.tight_layout()
    for suffix in ('png','pdf'):fig.savefig(out/f's21_fit.{suffix}',dpi=190)
    plt.close(fig)
    print('T2',name,'TARGET RMS',score['complex_RMS'],'FIT',fitted['C_fF'],fitted['G_fF'],flush=True)
    return score,fitted


def main():
    selected=ROOT/'geometry/taper/selected_taper2_v3.json'
    if selected.exists():
        saved=json.loads(selected.read_text())
        if saved['case']==CASE and saved['specification']==specification():
            _,expected=build(BASE,coupler=COUPLER)
            if saved['geometry']!=expected:raise ValueError('Use a fresh CASE for changed geometry')
            print('Already verified:',saved['project']);print('Plot:',saved['plot']);return
    if RUN_INDUCTANCE:measure(CASE,BASE)
    if RUN_RESPONSE:
        saved=ROOT/'geometry/taper/output'/CAMPAIGN/CASE/'Lcal/geometry.json'
        _,expected=build(BASE,isolation=True)
        if json.loads(saved.read_text())!=expected:raise ValueError('Cached calibration body differs')
        measured=json.loads((ROOT/'output_data/taper_exploration'/CAMPAIGN/CASE/'Lcal/results.json').read_text())
        if abs(measured['target_error_percent'])>2:raise ValueError('Inductance must pass before full response')
        score,fitted=response(CASE,BASE,measured['L0_nH'],coupler=COUPLER)
        if score['pass'] and fitted['curve_pass'] and FINALIZE_ON_PASS:finalize(CASE)


def finalize(name):
    """Audit and move the complete native case, never the verified T1 folder."""
    from geometry.taper.audit import check_geometry
    require_reference()
    out=ROOT/'output_data/taper_exploration'/CAMPAIGN/name
    score=json.loads((out/'response/target_comparison.json').read_text())
    fitted=json.loads((out/'response/extraction/fit.json').read_text())
    measured=json.loads((out/'Lcal/results.json').read_text())
    if (not score['pass'] or score['complex_RMS']>S21_COMPLEX_RMS_CUTOFF or
        not fitted['curve_pass'] or fitted['complex_RMS']>S21_COMPLEX_RMS_CUTOFF):
        raise ValueError('Fixed-target and fitted-curve checks must pass')
    if abs(100*(measured['L0_nH']/TAPER_2.L0_nH-1))>2:raise ValueError('Isolated L outside tolerance')
    if not np.isclose(fitted['L_nH'][0],measured['L0_nH'],rtol=1e-12):raise ValueError('Wrong fixed inductance')
    if score['specification']!=specification():raise ValueError('Superseded targets')
    source=(ROOT/'sonnet_files/taper_exploration'/CAMPAIGN/name).resolve()
    destination=(ROOT/'sonnet_files/taper_exploration/final_designs'/f'new_T2_{name}_final').resolve()
    if not source.is_relative_to(ROOT) or not destination.is_relative_to(ROOT):raise ValueError('Unsafe move path')
    if destination.exists():raise ValueError('Final directory already exists')
    metadata={}
    for stage,n in (('Lcal',1),('response',3)):
        gds=ROOT/'geometry/taper/output'/CAMPAIGN/name/stage/'geometry.gds'
        meta=json.loads(gds.with_suffix('.json').read_text());metadata[stage]=meta
        check_geometry(source/stage/f'generated_ncell_{n}.sonx',gds,meta,geometry_validator=validate)
    if metadata['Lcal']['unit']!=metadata['response']['unit']:raise ValueError('Calibration core changed')
    cal=source/'Lcal/generated_ncell_1.sonx'
    marker=json.loads((out/'Lcal/operator_completed.json').read_text())
    if W['digest'](cal)!=marker['project_sha256'] or W['digest'](cal.with_suffix('.snp'))!=marker['snp_sha256']:
        raise ValueError('Calibration source changed')
    for log in (source/'Lcal/em.log',out/'response/em.log'):
        text=log.read_text()
        if not all(s in text for s in ('Em errors detected: 0','Em warnings detected: 0','Em simulation completed')):
            raise ValueError('Native simulation incomplete or flagged')
    project=source/'response/generated_ncell_3.sonx';snp=project.with_suffix('.snp')
    run=json.loads((out/'response/run_status.json').read_text())
    if run['status']!='complete' or W['digest'](project)!=run['project_sha256'] or W['digest'](snp)!=score['source_sha256']:
        raise ValueError('Response source changed')
    if fitted['provenance']['sha256']!=score['source_sha256']:raise ValueError('Fit source differs')
    if list(ET.parse(project).iter(tag('Component'))):raise ValueError('Ideal components in production geometry')
    destination.parent.mkdir(exist_ok=True);source.rename(destination)
    for file in ('s21_fit.png','s21_fit.pdf','geometry.png','target_comparison.json'):
        shutil.copy2(out/'response'/file,destination/file)
    gds=ROOT/'geometry/taper/output'/CAMPAIGN/name/'response/geometry.gds'
    for suffix in ('.gds','.json'):shutil.copy2(gds.with_suffix(suffix),destination/('geometry'+suffix))
    unit=make_unit_cell(shared.from_dict(metadata['response']['unit']),name='NEW_T2')
    library=gdstk.Library(unit=1e-6,precision=1e-9);library.new_cell('TOP').add(gdstk.Reference(unit))
    library.add(unit,*unit.dependencies(True));library.write_gds(str(destination/'t2_unit.gds'))
    coupling=InputCoupler(**metadata['response']['coupler'])
    shared_junction_changed=coupling.d_regular!=66 or coupling.t_regular_ground!=6
    selected={'specification':specification(),'stage':'T2','campaign':CAMPAIGN,'case':name,
              'T1_accepted':True,'T2_accepted':True,'T3_accepted':False,
              'T1_reference':str(REFERENCE_FILE),'T1_reference_sha256':W['digest'](REFERENCE_FILE),
              'geometry':metadata['response'],'isolated_L':measured,'fixed_target_comparison':score,
              'heuristic_fit':fitted,'fixed_regular_extraction_model':REGULAR_MODEL,
              'project':str(destination/'response/generated_ncell_3.sonx'),'plot':str(destination/'s21_fit.png'),
              'gds':str(destination/'geometry.gds'),'unit_gds':str(destination/'t2_unit.gds'),
              'audit':{'GDS_native_XOR_empty':True,'T1_and_downstream_metal_unchanged':not shared_junction_changed,
                       'T1_regular_output_component_shapes_unchanged':True,
                       'shared_junction_adjusted':shared_junction_changed,
                       'shared_junction_gap_um':coupling.d_regular,'shared_junction_ground_um':coupling.t_regular_ground,
                       'mesh_um':[2,2],'materials_preserved':True,'source_hashes_checked':True,
                       'production_ideal_components':0,'calibration_core_matches':True,'native_errors':0,'native_warnings':0}}
    selected['code_sha256']={q:W['digest'](ROOT/q) for q in ('geometry/specifications.py',
        'geometry/taper/taper2_v3.py','geometry/taper/progressive_fitting.py','geometry/taper/measure_compact_inductance.py',
        'geometry/unit_cell.py','geometry/meander.py','geometry/prepare_sonnet.py','geometry/taper/sonnet.py')}
    for path in (ROOT/'geometry/taper/selected_taper2_v3.json',destination/'selection.json',
                 ROOT/'geometry/taper/three_tapers_status.json'):
        path.write_text(json.dumps(selected,indent=2)+'\n')
    (destination/'README.md').write_text(
        '# Second taper: accepted fixture\n\n'
        'Open `response/generated_ncell_3.sonx` for T2 + locked T1 + regular.\n'
        'Open `s21_fit.png` for Sonnet, complex fit, fixed nominal target and extracted-value errors.\n'
        f"Measured T2 L = {measured['L0_nH']:.8f} nH; fixed-target complex RMS = {score['complex_RMS']:.6f} (limit .025).\n"
        'The comparison covers 7501 frequencies over 4.65–5.40 GHz; only one matched fixture delay is adjusted for the target curve.\n'
        'The heuristic extraction allows all taper capacitances and fixture delay to vary, with every inductance and regular-cell value fixed.\n'
        'T1, regular and output component shapes match the accepted T1. Groups translate only by even-grid distances.\n'
        f'Assembly T1-regular junction: gap {coupling.d_regular} um, centered ground strip {coupling.t_regular_ground} um; '
        f'original accepted fixture used 66/6 um. Shared-junction adjustment: {shared_junction_changed}.\n'
        'No input–T2 or T2–T1 ground strip; flat upper ground. Production contains no ideal circuit elements.\n'
        'Independent weak-probe inductance calibration is in `Lcal/`. The added capacitors are evaluated offline from one four-port EM operator.\n'
        '`geometry.gds` contains the fixture; `t2_unit.gds` contains the isolated modular T2 body (input coupler is fixture geometry).\n'
        'See `selection.json` for complete parameters, measurements and provenance. T3 and the complete filter remain unverified.\n',encoding='utf-8')
    print('FINAL T2',destination,flush=True)
    return selected


if __name__=='__main__':main()
