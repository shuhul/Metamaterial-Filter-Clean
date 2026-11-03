"""No-IDC width-only meander controls for minimum_inductance.py.

Only the T2 inductor changes. The outer frame, side capacitor and native T1
polygons stay fixed. Centered 4 um trace, 6 um gaps, 58 runs, 600 um height.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from dataclasses import replace, asdict
from types import SimpleNamespace
import json
import xml.etree.ElementTree as ET
import gdstk
import numpy as np
from geometry.meander import MeanderParams, make_meander
from geometry.placement import place
from geometry.prepare_sonnet import main as prepare_native, tag, write_xml
from geometry.taper.compact_stages import ROOT, build as compact_build
from geometry.taper.second_layout import from_dict
from geometry.taper.bridged_fixtures import preview, require_native_locked
from geometry.taper.sonnet import add_capacitor, simulate, configure
from geometry.taper.measure_compact_inductance import calibrate
from geometry.taper.minimum_inductance import OUT, assess
from geometry.taper.compact_fit import fit_direct
from geometry.taper.run_bridged_stage import digest

NAME='minimum_inductance_2026-10-23'
NATIVE=ROOT/'sonnet_files/taper_exploration'/NAME
GDS=ROOT/'geometry/taper/output'/NAME
BASE=ROOT/'geometry/taper/output/compact_centered/w500_top10_wire52_gap6/direct/geometry.json'


def build(width,stage):
    if width%4 or not 16<=width<=440:
        raise ValueError('Centered meander width must be 16..440 um, multiple of 4')
    p=from_dict(json.loads(BASE.read_text())['second'])
    m=MeanderParams(w=width,l=600,l_top=10,n=58,t=4,d=6)
    # Use the existing compact fixture constructor, then replace only its wire.
    p=replace(p,taper=replace(p.taper,meander=m))
    polys,meta,overlays=compact_build(p,stage)
    x=meta['cell_centers'][1]
    cut=gdstk.rectangle((x-2,-618),(x+2,-18))
    polys=gdstk.boolean(polys,[cut],'not',precision=.001,layer=1)
    coil=place(make_meander(m,name='T2_COIL'),origin=(x,-18),anchor='top')
    polys.extend(coil.get_polygons())
    polys=gdstk.boolean(polys,[],'or',precision=.001,layer=1)
    meta.update(layout_kind='minimum_L_meander',meander=asdict(m),
                title=f'Direct-Fed T2 Meander: Width {width} um, No IDC',
                expected_metal_components=len(polys),width_only_trial=True,
                original_fixture=str(BASE),original_fixture_sha256=digest(BASE))
    return polys,json.loads(json.dumps(meta)),overlays


def validate(path,meta):
    expected,rebuilt,_=build(meta['meander']['w'],meta['stage'])
    if rebuilt!=meta:
        raise ValueError('Meander metadata changed')
    lib=gdstk.read_gds(str(path));tops=lib.top_level()
    if len(tops)!=1 or tops[0].name!='TOP' or not np.isclose(lib.unit,1e-6):
        raise ValueError('Invalid GDS')
    actual=tops[0].get_polygons()
    if any(q.layer!=1 or q.datatype!=0 for q in actual) or gdstk.boolean(actual,expected,'xor',precision=.001):
        raise ValueError('Unexpected meander geometry')
    p=meta['parameters']
    return SimpleNamespace(n_cells=p['n_cells'],l_sim=p['l_sim'],l_top_ground=p['l_top_ground']),actual


def prepare(width,stage):
    name=f'meander_w{width}'
    dest=NATIVE/name/stage;gd=GDS/name/stage;out=OUT/name/stage
    for folder in (dest,gd,out):folder.mkdir(parents=True,exist_ok=True)
    polys,meta,overlays=build(width,stage)
    count=meta['parameters']['n_cells']
    project=dest/f'generated_ncell_{count}.sonx'
    if not project.exists():
        lib=gdstk.Library(unit=1e-6,precision=1e-9);lib.new_cell('TOP').add(*polys)
        gds=gd/'geometry.gds';lib.write_gds(str(gds))
        gds.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
        prepare_native(count,gds_path=gds,output_dir=dest,geometry_validator=validate)
    elif json.loads((gd/'geometry.json').read_text())!=meta:
        raise ValueError('Changed geometry requires a new case')
    if stage=='direct':require_native_locked(project,meta)
    preview(polys,meta,overlays,dest/'geometry_review.png')
    return project,out,meta


def measure(width):
    project,out,meta=prepare(width,'Lcal')
    marker=out/'operator_prepared.json'
    if not marker.exists():
        (out/'prepared_source.sonx').write_bytes(project.read_bytes())
        add_capacitor(project,meta)
        tree=ET.parse(project);comp=next(tree.iter(tag('Component')))
        comp.set('Type','PORTS');comp.attrib.pop('Units',None)
        for q in list(comp):
            if q.tag==tag('Value'):comp.remove(q)
        sweeps=tree.find(tag('Sweeps'))
        for q in list(sweeps):sweeps.remove(q)
        fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
        for f in (4,5,5.5,6,7):
            ET.SubElement(fs,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(f),Stop=str(f),Step='1')
        for q in tree.findall(f"{tag('OutputFiles')}//{tag('File')}"):
            if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','FALSE')
        write_xml(tree,project);marker.write_text(json.dumps(dict(project_sha256=digest(project))))
    elif json.loads(marker.read_text())['project_sha256']!=digest(project):
        raise ValueError('Changed L operator')
    result=out/'results.json'
    if result.exists():return json.loads(result.read_text())
    return calibrate(project,out,f'meander_w{width}',auto_loads=True,desired_GHz=(5.0,4.7,4.4,4.1))


def response(width,inductance,*,fit_legacy=True):
    project,out,meta=prepare(width,'direct')
    marker=out/'run_status.json'
    if marker.exists() and json.loads(marker.read_text()).get('status')=='complete':
        saved=json.loads(marker.read_text())
        if saved['project_sha256']!=digest(project) or saved['snp_sha256']!=digest(project.with_suffix('.snp')):
            raise ValueError('Changed completed response')
        snp=project.with_suffix('.snp')
    else:
        tree=ET.parse(project);sweeps=tree.find(tag('Sweeps'))
        for q in list(sweeps):sweeps.remove(q)
        fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
        ET.SubElement(fs,tag('Sweep'),Type='ABS_ENTRY',On='TRUE',Adaptive='TRUE',Start='4.65',Stop='5.4',Target='7500')
        for q in tree.findall(f"{tag('OutputFiles')}//{tag('File')}"):
            if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','TRUE')
        write_xml(tree,project);configure(project,4.65,5.4,7500)
        snp=simulate(project,out)
    if not fit_legacy:return snp,out
    fitpath=out/'fit/fit.json'
    if fitpath.exists():r=json.loads(fitpath.read_text())
    else:r=fit_direct(snp,out/'fit',inductance,data_label=f'Sonnet meander width {width} um (no IDC)')
    assessment=assess(r,np.load(out/'fit/response.npz'))
    (out/'assessment.json').write_text(json.dumps(assessment,indent=2)+'\n')
    return r,assessment


if __name__=='__main__':
    import sys
    width=int(sys.argv[1])
    if len(sys.argv)>2 and sys.argv[2]=='prepare':prepare(width,'direct')
    else:
        r=measure(width)
        response(width,r['L0_nH'])
