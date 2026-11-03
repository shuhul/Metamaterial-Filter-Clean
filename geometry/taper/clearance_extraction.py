"""Separate the left grounded side block from the input electrode on a 2 um grid.

T1 and downstream metal translate rigidly. Their shapes and spacings do not
change. The centered straight wire still terminates at the original flat ground.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
from types import SimpleNamespace
import xml.etree.ElementTree as ET
import gdstk
import numpy as np
from geometry.taper.compact_stages import build as original_build, JunctionIDC
from geometry.taper.second_layout import from_dict
from geometry.taper.corrected_iteration import ROOT, review
from geometry.taper.bridged_fixtures import preview, require_native_locked
from geometry.prepare_sonnet import main as import_native, tag, write_xml
from geometry.taper.sonnet import add_capacitor, configure, simulate
from geometry.taper.measure_compact_inductance import calibrate
from geometry.taper.run_bridged_stage import digest

NAME='clearance_extraction_2026-10-25'
NATIVE=ROOT/'sonnet_files/taper_exploration'/NAME
GDS=ROOT/'geometry/taper/output'/NAME
OUT=ROOT/'output_data/taper_exploration'/NAME


def build(p,stage,extra,*,junction=None,recess_t1=0):
    if extra<0 or extra%2:raise ValueError('Additional left clearance must be nonnegative and even')
    polys,meta,overlays=original_build(p,stage,junction=junction)
    translated=[q.copy().translate(extra,0) for q in polys]
    side=412;bottom=-618
    opening=[gdstk.rectangle((side,bottom),(side+extra,-24))] if extra else []
    changed=gdstk.boolean(translated,opening,'not',precision=.001,layer=1)
    if extra:
        changed.extend([gdstk.rectangle((0,0),(extra,448),layer=1),
                        gdstk.rectangle((0,-1052),(extra,bottom),layer=1),
                        gdstk.rectangle((0,bottom),(extra,-24),layer=1),
                        gdstk.rectangle((0,-18),(extra,-8),layer=1)])
    changed=gdstk.boolean(changed,[],'or',precision=.001,layer=1)
    # No changes to any device, top ground, or bottom ground to the right of
    # the source electrode. The only cut is the left grounded side block.
    device_roi=[gdstk.rectangle((side+extra+14,-1052),(meta['w_sim']+extra,448))]
    if gdstk.boolean(gdstk.boolean(translated,device_roi,'and',precision=.001),
                     gdstk.boolean(changed,device_roi,'and',precision=.001),'xor',precision=.001):
        raise ValueError('Clearance change touched device metal')
    if recess_t1:
        if stage not in ('direct','idc') or recess_t1%2 or not 0<recess_t1<448:
            raise ValueError('T1 ground recess requires a response stage and an even depth below 448 um')
        t1=p.taper1
        left=meta['cell_centers'][2]+extra-t1.w/2
        opening=gdstk.rectangle((left,0),(left+t1.w,recess_t1))
        changed=gdstk.boolean(changed,[opening],'not',precision=.001,layer=1)
        meta['T1_top_ground_recess_um']=recess_t1
    meta['layout_kind']='left_clearance_extraction'
    meta['left_clearance_extra_um']=extra;meta['left_ground_clearance_um']=14+extra
    meta['w_sim']+=extra;meta['cell_centers']=[x+extra for x in meta['cell_centers']]
    meta['expected_metal_components']=len(changed)
    meta['title']=f"{'No IDC' if stage in ('direct','Lcal') else 'Uniform IDC'}: Left Ground Clearance {14+extra} um"
    meta['overlay_bounds']=[]
    return changed,json.loads(json.dumps(meta))


def validate(path,meta):
    junction=JunctionIDC(**meta['junction_idc']) if 'junction_idc' in meta else None
    expected,rebuilt=build(from_dict(meta['second']),meta['stage'],meta['left_clearance_extra_um'],junction=junction,
                           recess_t1=meta.get('T1_top_ground_recess_um',0))
    if rebuilt!=meta:raise ValueError('Stale clearance geometry metadata')
    lib=gdstk.read_gds(str(path));tops=lib.top_level()
    if len(tops)!=1 or tops[0].name!='TOP' or not np.isclose(lib.unit,1e-6):raise ValueError('Invalid GDS')
    actual=tops[0].get_polygons()
    if any(q.layer!=1 or q.datatype!=0 for q in actual):
        raise ValueError('Metal must be layer 1/datatype 0')
    # As in the original pipeline, regular-cell source GDS contains odd-um
    # coordinates. The native importer snaps them by at most 1 um and checks
    # every final vertex against the 2 um grid; native T1/downstream XOR is
    # then checked against the same already verified snapped reference.
    if gdstk.boolean(actual,expected,'xor',precision=.001):raise ValueError('GDS differs from parameters')
    return SimpleNamespace(n_cells=meta['parameters']['n_cells'],l_sim=1500,l_top_ground=448),actual


def prepare(case,p,stage,extra=200,*,campaign=None,junction=None,recess_t1=0):
    native=ROOT/'sonnet_files/taper_exploration'/campaign if campaign else NATIVE
    generated=ROOT/'geometry/taper/output'/campaign if campaign else GDS
    output=ROOT/'output_data/taper_exploration'/campaign if campaign else OUT
    dest=native/case/stage;gd=generated/case/stage;out=output/case/stage
    for folder in (dest,gd,out):folder.mkdir(parents=True,exist_ok=True)
    polys,meta=build(p,stage,extra,junction=junction,recess_t1=recess_t1)
    project=dest/f"generated_ncell_{meta['parameters']['n_cells']}.sonx"
    path=gd/'geometry.gds'
    if project.exists():
        if json.loads(path.with_suffix('.json').read_text())!=meta:raise ValueError('Changed geometry needs a new case')
    else:
        lib=gdstk.Library(unit=1e-6,precision=1e-9);lib.new_cell('TOP').add(*polys)
        lib.write_gds(str(path));path.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
        import_native(meta['parameters']['n_cells'],gds_path=path,output_dir=dest,geometry_validator=validate)
    if stage in ('direct','idc'):
        # Target-response iterations may shorten just T1's facing coupling
        # leg. Validate the rest of T1/downstream against the original metal.
        # The generated GDS/native checks still validate the changed leg.
        from geometry.taper.locked_taper2 import T1
        allowed=[]
        if p.x_junction_ground:
            if not campaign or not campaign.startswith(('capacitance_guided_','target_response_')):
                raise ValueError('Junction ground placement requires an explicit parameter-iteration campaign')
            # Only the existing grounded strip moves. The cap/meander shapes,
            # cell spacing, and all other downstream metal remain identical.
            # The GDS/native validation independently checks the entire geometry
            # against build(), including the intended strip inside this region.
            right=meta['cell_centers'][2]+T1.w/2
            allowed.append(gdstk.rectangle((right,0),(right+p.d_t1_regular,1500)))
        if recess_t1:
            if not campaign or not campaign.startswith('target_response_'):
                raise ValueError('Local ground recess requires a target-response campaign')
            left=meta['cell_centers'][2]-T1.w/2
            allowed.append(gdstk.rectangle((left,448-recess_t1),(left+T1.w,448)))
        if p.taper1 is not None and p.taper1.l_left!=T1.l_left:
            if not campaign or not campaign.startswith('target_response_'):
                raise ValueError('T1 coupling-leg changes require a target-response campaign')
            if p.taper1.t_left!=T1.t_left or p.taper1.l_left>T1.l_left:
                raise ValueError('Only shortening the original T1 left leg is supported')
            left=meta['cell_centers'][2]-T1.w/2
            allowed.append(gdstk.rectangle((left,448+8+p.taper1.l_left),
                                     (left+T1.t_left,448+8+T1.l_left)))
        require_native_locked(project,meta,allowed_regions=allowed)
    preview(polys,meta,[],dest/'geometry_review.png')
    return project,out,meta


def measure(case,p,stage='Lcal',extra=200,*,campaign=None,junction=None):
    project,out,meta=prepare(case,p,stage,extra,campaign=campaign,junction=junction)
    marker=out/'operator_prepared.json'
    if not marker.exists():
        add_capacitor(project,meta)
        tree=ET.parse(project);root=tree.getroot();comp=next(root.iter(tag('Component')))
        comp.set('Type','PORTS');comp.attrib.pop('Units',None)
        for q in list(comp):
            if q.tag==tag('Value'):comp.remove(q)
        sweeps=root.find(tag('Sweeps'))
        for q in list(sweeps):sweeps.remove(q)
        fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
        for f in (4,5,5.5,6,7):
            ET.SubElement(fs,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(f),Stop=str(f),Step='1')
        for q in root.iter(tag('File')):
            if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','FALSE')
        write_xml(tree,project);marker.write_text(json.dumps(dict(project_sha256=digest(project))))
    elif json.loads(marker.read_text())['project_sha256']!=digest(project):raise ValueError('Calibration changed')
    return calibrate(project,out,case,auto_loads=True)


def response(case,p,l,stage='direct',extra=200):
    if abs(l/.74-1)>.02:raise ValueError('Inductance must remain within 2% before full response')
    project,out,meta=prepare(case,p,stage,extra)
    marker=out/'run_status.json'
    if not marker.exists() or json.loads(marker.read_text()).get('status')!='complete':
        tree=ET.parse(project);root=tree.getroot();sweeps=root.find(tag('Sweeps'))
        for q in list(sweeps):sweeps.remove(q)
        fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
        ET.SubElement(fs,tag('Sweep'),Type='ABS_ENTRY',On='TRUE',Adaptive='TRUE',Start='4.65',Stop='5.4',Target='7500')
        for q in root.iter(tag('File')):
            if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','TRUE')
        write_xml(tree,project);configure(project,4.65,5.4,7500);source=simulate(project,out)
    else:
        record=json.loads(marker.read_text())
        if record['project_sha256']!=digest(project) or record['snp_sha256']!=digest(project.with_suffix('.snp')):
            raise ValueError('Completed response changed')
        source=project.with_suffix('.snp')
    return review(source,out/'fit',l,idc=stage=='idc')


def parameters():
    return from_dict(json.loads((ROOT/'geometry/taper/output/compact_centered/plain_w700_wire68_gap6/direct/geometry.json').read_text())['second'])


if __name__=='__main__':
    import sys
    p=parameters();case='frame700_wire68_left214'
    if '--prepare' in sys.argv:prepare(case,p,'direct')
    else:
        c=measure(case,p)
        response(case,p,c['L0_nH'])
