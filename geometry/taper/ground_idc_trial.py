"""Prepare and validate an optional physical ground-IDC trial.

Importing this module runs nothing. It is used only after a verified coupled
response indicates a useful ground-capacitance change; it is not a final
design generator or an acceptance shortcut.
"""
from dataclasses import asdict
import json
from types import SimpleNamespace
import xml.etree.ElementTree as ET
import gdstk
import numpy as np
import matplotlib.pyplot as plt
from geometry.taper.compact_stages import ROOT,NATIVE,GDS,OUT,build as base_build
from geometry.taper.ground_idc import GroundIDCParams,apply
from geometry.taper.second_layout import from_dict
from geometry.taper.bridged_fixtures import preview,require_native_locked
from geometry.prepare_sonnet import main as prepare_native,tag,write_xml
from geometry.taper.sonnet import add_capacitor,configure,simulate
from geometry.taper.run_bridged_stage import digest
from geometry.taper.measure_compact_inductance import calibrate


def build(p,ground,stage):
    if p.w_upper_plate:raise ValueError('Rejected upper strips cannot be used')
    polys,meta,overlays=base_build(p,stage)
    polys,meta=apply(polys,meta,ground);meta['layout_kind']='compact_ground_idc_trial'
    return polys,meta,overlays


def validate(gds,meta):
    polys,rebuilt,_=build(from_dict(meta['second']),GroundIDCParams(**meta['ground_idc']),meta['stage'])
    if meta!=rebuilt:raise ValueError('Ground-IDC metadata changed')
    lib=gdstk.read_gds(str(gds));tops=lib.top_level()
    if len(tops)!=1 or tops[0].name!='TOP' or not np.isclose(lib.unit,1e-6,rtol=1e-9,atol=0):raise ValueError('Invalid GDS library')
    actual=tops[0].get_polygons()
    if any(q.layer!=1 or q.datatype!=0 for q in actual):raise ValueError('Wrong layer')
    if gdstk.boolean(actual,polys,'xor',precision=.001):raise ValueError('GDS differs from trial parameters')
    p=meta['parameters']
    return SimpleNamespace(n_cells=p['n_cells'],l_sim=p['l_sim'],l_top_ground=p['l_top_ground']),actual


def prepare(name,p,ground,stage='idc'):
    if stage not in ('direct','idc','Lcal','idc_Lcal'):raise ValueError('Unknown trial stage')
    dest=NATIVE/name/stage;gd=GDS/name/stage;out=OUT/name/stage
    for path in (dest,gd,out):path.mkdir(parents=True,exist_ok=True)
    polys,meta,overlays=build(p,ground,stage)
    project=dest/f"generated_ncell_{meta['parameters']['n_cells']}.sonx"
    if project.exists():
        if json.loads((gd/'geometry.json').read_text())!=meta:raise ValueError('Changed trial needs a new name')
    else:
        lib=gdstk.Library(unit=1e-6,precision=1e-9);lib.new_cell('TOP').add(*polys)
        gds=gd/'geometry.gds';lib.write_gds(str(gds));gds.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
        prepare_native(meta['parameters']['n_cells'],gds_path=gds,output_dir=dest,geometry_validator=validate)
    if stage in ('direct','idc'):require_native_locked(project,meta)
    fig,ax=plt.subplots(figsize=(16,7));preview(polys,meta,overlays,None,ax)
    ax.set_ylim(-665,max(180,ground.l+40));fig.tight_layout();fig.savefig(dest/'geometry_review.png',dpi=190);plt.close(fig)
    return project,out,meta


def response(name,p,ground,stage='idc'):
    project,out,meta=prepare(name,p,ground,stage);marker=out/'run_status.json'
    if marker.exists() and json.loads(marker.read_text()).get('status')=='complete':
        saved=json.loads(marker.read_text())
        if saved['project_sha256']!=digest(project) or saved['snp_sha256']!=digest(project.with_suffix('.snp')):raise ValueError('Completed trial changed')
        return project.with_suffix('.snp'),out
    configure(project,4.65,5.4,7500)
    return simulate(project,out),out


def measure(name,p,ground):
    project,out,meta=prepare(name,p,ground,'idc_Lcal');marker=out/'operator_prepared.json'
    if not marker.exists():
        (out/'prepared_source.sonx').write_bytes(project.read_bytes());add_capacitor(project,meta)
        tree=ET.parse(project);comp=next(tree.iter(tag('Component')))
        comp.set('Type','PORTS');comp.attrib.pop('Units',None)
        for q in list(comp):
            if q.tag==tag('Value'):comp.remove(q)
        sweeps=tree.find(tag('Sweeps'))
        for q in list(sweeps):sweeps.remove(q)
        fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
        for f in (4,5,5.5,6,7):ET.SubElement(fs,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(f),Stop=str(f),Step='1')
        for q in tree.iter(tag('File')):
            if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','FALSE')
        write_xml(tree,project);marker.write_text(json.dumps(dict(project_sha256=digest(project)),indent=2)+'\n')
    elif json.loads(marker.read_text())['project_sha256']!=digest(project):raise ValueError('Calibration operator changed')
    # Large physical shunt C can put the original 6.4 GHz calibration target
    # above the unloaded resonance. Lower positive loads use the SAME measured
    # 4–7 GHz operator, with 5.5 GHz still held out of interpolation.
    return calibrate(project,out,name,auto_loads=True,desired_GHz=(5.0,4.7,4.4,4.1))
