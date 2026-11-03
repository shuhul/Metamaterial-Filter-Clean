"""Frank's T1-to-regular coupling controls; no T2 or input IDC.

Retain the exact right coupling leg, gap ground strip and regular/downstream
geometry. Remove T1's resonator, comparing its full frame with just the driven
coupling electrode. Every EM run is serialized and limited to four threads.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import xml.etree.ElementTree as ET
from types import SimpleNamespace
import gdstk
import numpy as np
from geometry.taper.bridged_fixtures import ROOT,build as original_build,preview
from geometry.prepare_sonnet import main as prepare_native,tag,points,write_xml
from geometry.taper.sonnet import simulate
from geometry.taper.run_bridged_stage import digest

NAME='frank_coupling_controls_2026-10-23'
NATIVE=ROOT/'sonnet_files/taper_exploration'/NAME
GDS=ROOT/'geometry/taper/output'/NAME
OUT=ROOT/'output_data/taper_exploration'/NAME
REFERENCE=ROOT/'sonnet_files/taper_exploration/bridged_feed_tests/01_t1_input_bridged/generated_ncell_2.sonx'
CASES=('minimal_feed10','minimal_feed14','frame_without_inductor')


def build(case):
    if case not in CASES:raise ValueError('Unknown control')
    polys,meta,overlays=original_build('01_t1_input_bridged')
    p=meta['taper'];left=meta['cell_centers'][1]-p['w']/2;right=left+p['w']
    # Stop exactly at the bottom-ground boundary and cap underside. Preserve
    # the 6 um right electrode in full, including its top corner.
    cut_left=left+p['t_left'] if case=='frame_without_inductor' else left
    cuts=[gdstk.rectangle((cut_left,-618),(right-p['t_right'],-8-p['t_top']))]
    if case=='minimal_feed10':
        cuts.append(gdstk.rectangle((left,-8-p['t_top']),(right-p['t_right'],-18)))
    polys=gdstk.boolean(polys,cuts,'not',precision=.001,layer=1)
    meta.update(layout_kind='frank_coupling_only',control=case,T1_inductor_removed=True,
                T1_left_cap_leg_removed=case!='frame_without_inductor',
                retained_coupling_leg=dict(x_left=right-p['t_right'],x_right=right,t=p['t_right'],l=p['l_right']),
                title='T1 Removed: Coupling Electrode Only' if case.startswith('minimal') else 'T1 Frame Retained: Inductor Removed',
                feed_width_at_removed_T1=10 if case=='minimal_feed10' else p['t_top'],
                expected_metal_components=len(polys),
                simulation_status='prepared_control',source_reference=str(REFERENCE),
                source_reference_sha256=digest(REFERENCE))
    return polys,json.loads(json.dumps(meta)),overlays


def validate(path,meta):
    expected,rebuilt,_=build(meta['control'])
    if rebuilt!=meta:raise ValueError('Control metadata changed')
    lib=gdstk.read_gds(str(path));tops=lib.top_level()
    if len(tops)!=1 or tops[0].name!='TOP' or not np.isclose(lib.unit,1e-6):raise ValueError('Invalid GDS')
    actual=tops[0].get_polygons()
    if any(p.layer!=1 or p.datatype!=0 for p in actual) or gdstk.boolean(actual,expected,'xor',precision=.001):
        raise ValueError('GDS does not match control')
    p=meta['parameters']
    return SimpleNamespace(n_cells=p['n_cells'],l_sim=p['l_sim'],l_top_ground=p['l_top_ground']),actual


def check_native(project,meta):
    # Include the complete coupling electrode, gap/ground strip and all metal
    # to its right in a strict XOR. No exceptions in the retained region.
    x=meta['retained_coupling_leg']['x_left']
    window=gdstk.rectangle((x,0),(meta['w_sim'],1500))
    metal=lambda path:[gdstk.Polygon(points(q)) for q in ET.parse(path).iter(tag('PlanarPolygon'))]
    a=gdstk.boolean(metal(project),[window],'and',precision=.001)
    b=gdstk.boolean(metal(REFERENCE),[window],'and',precision=.001)
    if gdstk.boolean(a,b,'xor',precision=.001):raise ValueError('Coupling electrode or downstream geometry changed')
    return dict(retained_electrode_gap_ground_regular_end_XOR_empty=True,mesh_um=[2,2],
                external_port_ohm=50,reference_planes_um=400,reference_project_sha256=digest(REFERENCE))


def prepare(case):
    dest=NATIVE/case;gd=GDS/case;out=OUT/case
    for folder in (dest,gd,out):folder.mkdir(parents=True,exist_ok=True)
    polys,meta,overlays=build(case);project=dest/'generated_ncell_2.sonx'
    if not project.exists():
        lib=gdstk.Library(unit=1e-6,precision=1e-9);lib.new_cell('TOP').add(*polys)
        path=gd/'geometry.gds';lib.write_gds(str(path))
        path.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
        prepare_native(2,gds_path=path,output_dir=dest,geometry_validator=validate)
    elif json.loads((gd/'geometry.json').read_text())!=meta:raise ValueError('Changed control geometry')
    checks=check_native(project,meta)
    (out/'geometry_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    preview(polys,meta,overlays,dest/'geometry_review.png')
    return project,out,meta


def run(case):
    project,out,meta=prepare(case)
    marker=out/'run_status.json'
    if marker.exists() and json.loads(marker.read_text()).get('status')=='complete':
        old=json.loads(marker.read_text())
        if old['project_sha256']!=digest(project) or old['snp_sha256']!=digest(project.with_suffix('.snp')):
            raise ValueError('Completed control changed')
        return project.with_suffix('.snp')
    tree=ET.parse(project);sweeps=tree.find(tag('Sweeps'))
    for q in list(sweeps):sweeps.remove(q)
    fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
    ET.SubElement(fs,tag('Sweep'),Type='ABS_ENTRY',On='TRUE',Adaptive='TRUE',Start='4.75',Stop='5.3',Target='5500')
    for q in tree.iter(tag('File')):
        if q.get('Type')=='TOUCH':q.set('Name',project.with_suffix('.snp').name);q.set('UseInterpolatedData','TRUE')
    write_xml(tree,project)
    return simulate(project,out)


def validate_peak_flanks(case='minimal_feed10'):
    """Two held-out direct EM points check ABS interpolation near half power."""
    import shutil
    from circuit_common import read_snp
    from geometry.taper.response_cases import W
    source=NATIVE/case/'generated_ncell_2.sonx'
    dest=NATIVE/(case+'_direct_validation');dest.mkdir(exist_ok=True)
    out=OUT/(case+'_direct_validation');out.mkdir(exist_ok=True)
    marker=out/'validation.json'
    if marker.exists():return json.loads(marker.read_text())
    project=dest/source.name
    if not project.exists():
        shutil.copytree(source.parent/'sondata',dest/'sondata',dirs_exist_ok=True)
        tree=ET.parse(source);sweeps=tree.find(tag('Sweeps'))
        for q in list(sweeps):sweeps.remove(q)
        fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
        for f in (5.0123,5.0163):
            ET.SubElement(fs,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(f),Stop=str(f),Step='1')
        for q in tree.iter(tag('File')):
            if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','FALSE')
        write_xml(tree,project)
    meta=json.loads((GDS/case/'geometry.json').read_text());check_native(project,meta)
    direct=W['em'](project,2);ref=read_snp(source.with_suffix('.snp'))
    indices=np.array([np.argmin(abs(ref.f-f)) for f in direct.f])
    if np.max(abs(ref.f[indices]-direct.f))>1:raise ValueError('No exact matching ABS frequency')
    error=abs(direct.s-ref.s[indices])
    r=dict(project=str(project),project_sha256=digest(project),snp_sha256=digest(project.with_suffix('.snp')),
           source_project_sha256=digest(source),source_snp_sha256=digest(source.with_suffix('.snp')),
           frequencies_GHz=(direct.f/1e9).tolist(),S21_max_error=float(error[:,1,0].max()),
           all_S_max_error=float(error.max()),passed=bool(error.max()<.003),
           direct_magnitudes=abs(direct.s[:,1,0]).tolist(),ABS_magnitudes=abs(ref.s[indices,1,0]).tolist())
    marker.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r),flush=True)
    return r


if __name__=='__main__':
    import sys
    case=sys.argv[1] if len(sys.argv)>1 else 'minimal_feed10'
    if len(sys.argv)>2 and sys.argv[2]=='prepare':prepare(case)
    else:run(case)
