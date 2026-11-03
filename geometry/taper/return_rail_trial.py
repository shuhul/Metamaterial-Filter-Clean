"""Local T2 ground returns, without a ground pedestal or a shortened wire.

Two separate narrow rails flank the full-height straight wire. Their upper
ends clear the top bar; neither crosses the signal wire or the strong junction.
"""
from dataclasses import replace
from types import SimpleNamespace
import json
import xml.etree.ElementTree as ET
import numpy as np
import gdstk
from geometry.taper.progressive_trial import parameters,review
from geometry.taper.progressive_fixtures import ROOT,GDS,NATIVE,OUT,prepare as prepare_output,run as run_output
from geometry.taper.compact_stages import build as compact_build
from geometry.taper.second_layout import from_dict
from geometry.prepare_sonnet import main as import_native,tag,write_xml
from geometry.taper.sonnet import add_capacitor
from geometry.taper.measure_compact_inductance import calibrate
from geometry.taper.run_bridged_stage import digest


def configuration(wire=8,gap=100,rail=40,clear=12):
    p=parameters(frame=500,wire=wire,pairs=48,span=86,gap=2,left=600,right=536,junction=6)
    p=replace(p,taper=replace(p.taper,t_right=20))
    return p,dict(gap=gap,t=rail,top_clearance=clear)


def build(p,rails,stage):
    if any(v<=0 or v%2 for v in rails.values()):raise ValueError('Rail dimensions must fit the 2 um grid')
    polys,meta,_=compact_build(p,stage)
    x=meta['cell_centers'][1];inside=p.taper.meander.t/2+rails['gap'];outside=inside+rails['t']
    if outside>=min(p.taper.w/2-p.taper.t_left,p.taper.w/2-p.taper.t_right)-2:
        raise ValueError('Ground return touches a cap leg')
    top=-18-rails['top_clearance']
    bounds=[((x-outside,-618),(x-inside,top)),((x+inside,-618),(x+outside,top))]
    polys.extend(gdstk.rectangle(*b) for b in bounds)
    polys=gdstk.boolean(polys,[],'or',precision=.001,layer=1)
    meta.update(layout_kind='local_return_rail_trial',return_rails=rails,return_rail_bounds=bounds,
                expected_metal_components=len(polys),raised_ground_box=False,
                title='T2: Full Straight Wire with Local Ground Returns')
    return polys,json.loads(json.dumps(meta))


def validate(path,meta):
    expected,rebuilt=build(from_dict(meta['second']),meta['return_rails'],meta['stage'])
    if rebuilt!=meta:raise ValueError('Changed rail metadata')
    lib=gdstk.read_gds(str(path));tops=lib.top_level()
    if len(tops)!=1 or tops[0].name!='TOP' or not np.isclose(lib.unit,1e-6):raise ValueError('Invalid rail library')
    actual=tops[0].get_polygons()
    if gdstk.boolean(actual,expected,'xor',precision=.001):raise ValueError('Rail GDS mismatch')
    q=meta['parameters'];return SimpleNamespace(n_cells=q['n_cells'],l_sim=q['l_sim'],l_top_ground=q['l_top_ground']),actual


def source(name,p,rails,stage):
    folder=GDS/name/'source'/stage;folder.mkdir(parents=True,exist_ok=True)
    path=folder/'geometry.gds';polys,meta=build(p,rails,stage)
    if path.exists():
        if json.loads(path.with_suffix('.json').read_text())!=meta:raise ValueError('Changed existing rail source')
    else:
        lib=gdstk.Library(unit=1e-6,precision=1e-9);lib.new_cell('TOP').add(*polys);lib.write_gds(str(path))
        path.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
    return path,meta


def prepare(name,p,rails,stage='direct'):
    path,meta=source(name,p,rails,stage)
    return prepare_output(name+'/'+stage,str(path),3,302)


def measure(name,p,rails):
    path,meta=source(name,p,rails,'idc_Lcal')
    folder=NATIVE/name/'Lcal';out=OUT/name/'Lcal';out.mkdir(parents=True,exist_ok=True)
    project=folder/'generated_ncell_1.sonx'
    if not project.exists():import_native(1,gds_path=path,output_dir=folder,geometry_validator=validate)
    marker=out/'prepared.json'
    if not marker.exists():
        add_capacitor(project,meta);tree=ET.parse(project);comp=next(tree.iter(tag('Component')))
        comp.set('Type','PORTS');comp.attrib.pop('Units',None)
        for child in list(comp):
            if child.tag==tag('Value'):comp.remove(child)
        sweeps=tree.find(tag('Sweeps'))
        for child in list(sweeps):sweeps.remove(child)
        fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
        for f in (4,5,5.5,6,7):ET.SubElement(fs,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(f),Stop=str(f),Step='1')
        for q in tree.iter(tag('File')):
            if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','FALSE')
        write_xml(tree,project);marker.write_text(json.dumps(dict(project_sha256=digest(project)),indent=2)+'\n')
    elif json.loads(marker.read_text())['project_sha256']!=digest(project):raise ValueError('Rail L operator changed')
    if (out/'results.json').exists():return json.loads((out/'results.json').read_text())
    return calibrate(project,out,name,auto_loads=True,desired_GHz=(5.,4.7,4.4,4.1))


def run(name,p,rails):
    prepare(name,p,rails)
    l=measure(name,p,rails)['L0_nH']
    if abs(l/.74-1)>.15:raise ValueError(f'Return-rail L={l:.6f} nH: adjust wire/rail spacing before full response')
    path,_=source(name,p,rails,'direct')
    snp,out=run_output(name+'/direct',str(path),3,302)
    from geometry.taper.progressive_fitting import fit
    from geometry.taper.staged_fitting import LT1,LR
    baseline=json.loads((OUT/'01_T1_regular/twoC/fit.json').read_text())
    r=fit(snp,out/'junction_check',[l,LT1,LR],[0,baseline['C_fF'][0],190.4],
          [None,56.7,8.2865],[('G',1)],title='T2 Local Returns: No-IDC Coupling Check')
    if not r['curve_pass'] or abs(r['G_fF'][1]/56.7-1)>.1:
        raise ValueError('No-IDC coupling check misses 10% exploratory gate; do not add the IDC yet')
    path,_=source(name,p,rails,'idc');prepare(name,p,rails,'idc')
    run_output(name+'/idc',str(path),3,302)
    return review(name+'/idc',l,junction_cg=r['G_fF'][1])


if __name__=='__main__':
    import sys
    name='rails_w500_wire8_gap100_t40_idc48_span86'
    p,rails=configuration()
    if '--prepare-only' in sys.argv:prepare(name,p,rails)
    else:run(name,p,rails)
