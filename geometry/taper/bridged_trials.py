"""Prepare and run new variants of the approved feed-driven geometries."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
from pathlib import Path
import xml.etree.ElementTree as ET
import gdstk
import numpy as np
from geometry.taper.bridged_fixtures import ROOT,build,validate,preview,CASES,require_native_locked
from geometry.prepare_sonnet import main as prepare_native
from geometry.taper.run_bridged_stage import OUT,digest
from geometry.taper.sonnet import configure,simulate
from geometry.prepare_sonnet import tag,write_xml


def prepare(name,p,stage=3,*,omit_wire=False):
    case=CASES[stage-1]
    native=ROOT/'sonnet_files/taper_exploration/bridged_trials'/name/case
    gd=ROOT/'geometry/taper/output/bridged_trials'/name/case
    out=OUT/'trials'/name/case
    project=native/'generated_ncell_3.sonx'
    polygons,meta,overlays=build(case,second=p,omit_t2_wire=omit_wire)
    if project.exists():
        stored=json.loads((gd/'geometry.json').read_text())
        if meta!=stored:raise ValueError('Changed trial parameters require a new name')
    else:
        gd.mkdir(parents=True,exist_ok=True);out.mkdir(parents=True,exist_ok=True)
        lib=gdstk.Library(unit=1e-6,precision=1e-9);lib.new_cell('TOP').add(*polygons)
        gds=gd/'geometry.gds';lib.write_gds(str(gds))
        gds.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
        prepare_native(3,gds_path=gds,output_dir=native,geometry_validator=validate)
    require_native_locked(project,meta)
    preview(polygons,meta,overlays,native/'geometry_review.png')
    return project,out


def run(name,p,stage=3,start=4.65,stop=5.4,target=7500,*,omit_wire=False):
    project,out=prepare(name,p,stage,omit_wire=omit_wire)
    out.mkdir(parents=True,exist_ok=True)
    marker=out/'run_status.json';record=json.loads(marker.read_text()) if marker.exists() else {}
    if record.get('status')=='complete':
        if digest(project)!=record['project_sha256'] or digest(project.with_suffix('.snp'))!=record['snp_sha256']:
            raise ValueError('Completed trial source/export changed')
        sweeps=ET.parse(project).findall(f"{tag('Sweeps')}//{tag('Sweep')}")
        if len(sweeps)!=1 or [float(sweeps[0].get(k)) for k in ('Start','Stop','Target')]!=[start,stop,target]:
            raise ValueError('Changed sweep requires a separate analysis record')
        return project.with_suffix('.snp'),out
    if not (out/'prepared_source.sonx').exists():(out/'prepared_source.sonx').write_bytes(project.read_bytes())
    # A fast direct-frequency screen may precede the adaptive validation.
    tree=ET.parse(project);blocks=tree.find(tag('Sweeps'))
    for q in list(blocks):blocks.remove(q)
    fs=ET.SubElement(ET.SubElement(blocks,tag('Set'),On='TRUE'),tag('Frequencies'))
    ET.SubElement(fs,tag('Sweep'),Type='ABS_ENTRY',On='TRUE',Adaptive='TRUE',Start=str(start),Stop=str(stop),Target=str(target))
    for q in tree.findall(f"{tag('OutputFiles')}//{tag('File')}"):
        if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','TRUE')
    write_xml(tree,project);configure(project,start,stop,target)
    return simulate(project,out),out


def screen(name,p,frequencies=(5.4,),*,omit_wire=False,stage=3):
    """One or two EM points to reject a bad direction before a full sweep."""
    from geometry.taper.response_cases import W
    project,out=prepare(name,p,stage,omit_wire=omit_wire)
    out.mkdir(parents=True,exist_ok=True)
    marker=out/'run_status.json'
    if marker.exists() and json.loads(marker.read_text()).get('status')=='complete':
        raise ValueError('Do not replace a validated full sweep with a frequency screen')
    tree=ET.parse(project);blocks=tree.find(tag('Sweeps'))
    for q in list(blocks):blocks.remove(q)
    fs=ET.SubElement(ET.SubElement(blocks,tag('Set'),On='TRUE'),tag('Frequencies'))
    for f in frequencies:ET.SubElement(fs,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(f),Stop=str(f),Step='1')
    for q in tree.findall(f"{tag('OutputFiles')}//{tag('File')}"):
        if q.get('Type')=='TOUCH':q.set('Name',project.with_suffix('.snp').name);q.set('UseInterpolatedData','FALSE')
    write_xml(tree,project)
    net=W['em'](project,2)
    stamp='_'.join(str(f).replace('.','p') for f in frequencies)
    saved=out/f'screen_{stamp}.snp';saved.write_bytes(project.with_suffix('.snp').read_bytes())
    report=dict(frequencies_GHz=(net.f/1e9).tolist(),S21_magnitude=abs(net.s[:,1,0]).tolist(),
        S21_phase_deg=np.angle(net.s[:,1,0],deg=True).tolist(),
        project_sha256=digest(project),snp_sha256=digest(saved),physical_acceptance=False)
    (out/f'screen_{stamp}.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True)
    return report


if __name__=='__main__':
    from geometry.taper.second_layout import from_dict
    p=from_dict(json.loads((ROOT/'geometry/taper/output/bridged_t2_w52_flat.json').read_text())['second'])
    result=run('w52_span290_right548',p)
    print(result,flush=True)
