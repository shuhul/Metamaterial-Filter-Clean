"""Guarded, cached execution for the centred no-IDC stages."""
import json
from copy import deepcopy
import xml.etree.ElementTree as ET
import gdstk
from geometry.taper.compact_stages import prepare
from geometry.taper.measure_compact_inductance import calibrate
from geometry.taper.sonnet import add_capacitor,simulate,configure
from geometry.taper.run_bridged_stage import digest
from geometry.prepare_sonnet import tag,write_xml,BASE,BIN,points,run as native_command


def add_known_junction(project,meta,value):
    """Diagnostic ideal C between facing top-bar edges; changes no metal."""
    tree=ET.parse(project);root=tree.getroot()
    existing=[q for q in root.iter(tag('Component')) if q.get('MacroID')=='compactKnownCg']
    if existing:
        if float(existing[0].findtext(tag('Value')))!=value/1000:raise ValueError('Changed ideal control C')
        return
    x0=meta['cell_centers'][1]+meta['second']['taper']['w']/2
    x1=meta['cell_centers'][2]-meta['second']['taper1']['w']/2
    y=460
    polys=[gdstk.Polygon(points(q)) for q in root.iter(tag('PlanarPolygon'))]
    if gdstk.inside([(x0-1,y),(x0+1,y),(x1-1,y),(x1+1,y)],polys)!=(True,False,False,True):
        raise ValueError('Known C pins do not face the intended top-bar edges')
    source=ET.parse(BASE/'sonnet_files/filter_v6.sonx')
    comp=deepcopy(next(q for q in source.iter(tag('Component')) if q.findtext(tag('Value'))=='C_add'))
    comp.set('MacroID','compactKnownCg');comp.set('Id',str(max(int(q.get('Id')) for q in root.iter(tag('PlanarPolygon')))+1))
    comp.find(tag('Value')).text=str(value/1000)
    comp.find(tag('LabelPos')).text=f'({(x0+x1)/2},{y-12})'
    comp.find(tag('SymbolBox')).text=f'({x0+2},{y-2},{x1-2},{y+2})'
    for i,(pin,x,orientation) in enumerate(zip(comp.findall(tag('Pin')),(x0,x1),('LEFT','RIGHT'))):
        pin.set('MacroID',f'compactCgPin{i}');pin.set('Orientation',orientation)
        pin.find(tag('Center')).text=f'({x},{y})'
    level=next(q for q in root.find(tag('Geometry')).findall(tag('Level')) if q.get('Number')=='0')
    level.find(tag('Objects')).append(comp);write_xml(tree,project)
    native_command([BIN/'soncmd.exe','-Read',project],project.parent/'known_cap_read.log')


def measure_L(name,p,*,junction=None,stage='Lcal'):
    if stage not in ('Lcal','idc_Lcal'):raise ValueError('Invalid calibration stage')
    project,out,meta=prepare(name,p,stage,junction=junction)
    prepared=out/'operator_prepared.json'
    if not prepared.exists():
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
        write_xml(tree,project);prepared.write_text(json.dumps(dict(project_sha256=digest(project))))
    elif json.loads(prepared.read_text())['project_sha256']!=digest(project):raise ValueError('L operator changed')
    return calibrate(project,out,name,auto_loads=True)


def response(name,p,*,omit_wire=False,start=4.65,stop=5.4,target=7500,ideal_cg_fF=None,junction=None,stage='direct'):
    if stage not in ('direct','idc'):raise ValueError('Invalid response stage')
    if stage=='idc' and (ideal_cg_fF is not None or omit_wire):raise ValueError('IDC stage requires the physical wire and junction')
    project,out,meta=prepare(name,p,stage,omit_wire=omit_wire,junction=junction)
    marker=out/'run_status.json';record=json.loads(marker.read_text()) if marker.exists() else {}
    if record.get('status')=='complete':
        if record['project_sha256']!=digest(project) or record['snp_sha256']!=digest(project.with_suffix('.snp')):
            raise ValueError('Completed response changed')
        sweeps=ET.parse(project).findall(f"{tag('Sweeps')}//{tag('Sweep')}")
        if len(sweeps)!=1 or [float(sweeps[0].get(k)) for k in ('Start','Stop','Target')]!=[start,stop,target]:
            raise ValueError('Changed sweep requires a new case')
        intent=out/'known_cap_intent.json'
        old_cap=json.loads(intent.read_text())['ideal_Cg_fF'] if intent.exists() else None
        if old_cap!=ideal_cg_fF:raise ValueError('Changed ideal-component control requires a new case')
        return project.with_suffix('.snp'),out
    if not (out/'prepared_source.sonx').exists():(out/'prepared_source.sonx').write_bytes(project.read_bytes())
    if ideal_cg_fF is not None:
        add_known_junction(project,meta,ideal_cg_fF)
        (out/'known_cap_intent.json').write_text(json.dumps(dict(ideal_Cg_fF=ideal_cg_fF,
            diagnostic_only=True,physical_Cg_electrode_removed=p.taper.l_right==12,IDC_present=False),indent=2)+'\n')
    tree=ET.parse(project);sweeps=tree.find(tag('Sweeps'))
    for q in list(sweeps):sweeps.remove(q)
    fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
    ET.SubElement(fs,tag('Sweep'),Type='ABS_ENTRY',On='TRUE',Adaptive='TRUE',Start=str(start),Stop=str(stop),Target=str(target))
    for q in tree.findall(f"{tag('OutputFiles')}//{tag('File')}"):
        if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','TRUE')
    write_xml(tree,project);configure(project,start,stop,target)
    return simulate(project,out),out


def screen(name,p,frequencies=(5.,),*,omit_wire=False,bridge_junction=False,junction=None,stage='direct'):
    from geometry.taper.response_cases import W
    if stage not in ('direct','idc'):raise ValueError('Invalid screen stage')
    if stage=='idc' and omit_wire:raise ValueError('IDC screen retains the physical wire')
    project,out,meta=prepare(name,p,stage,omit_wire=omit_wire,bridge_junction=bridge_junction,junction=junction)
    marker=out/'screen_completed.json'
    if marker.exists():
        record=json.loads(marker.read_text())
        saved=out/'screen.snp';saved_project=out/'screen_project.sonx'
        if not saved.exists() or not saved_project.exists():
            if record['project_sha256']!=digest(project) or record['snp_sha256']!=digest(project.with_suffix('.snp')):
                raise ValueError('Old screen was superseded by a full export; use the completed full response')
            saved.write_bytes(project.with_suffix('.snp').read_bytes());saved_project.write_bytes(project.read_bytes())
        if record['project_sha256']!=digest(saved_project) or record['snp_sha256']!=digest(saved):raise ValueError('Screen archive changed')
        if record['frequencies_GHz']!=list(frequencies):raise ValueError('Use a separate screen case')
        return W['read_network'](saved,2),out
    completed=out/'run_status.json'
    if completed.exists() and json.loads(completed.read_text()).get('status')=='complete':
        raise ValueError('Do not replace a completed full response with a screen')
    tree=ET.parse(project);sweeps=tree.find(tag('Sweeps'))
    for q in list(sweeps):sweeps.remove(q)
    fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
    for f in frequencies:ET.SubElement(fs,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(f),Stop=str(f),Step='1')
    for q in tree.findall(f"{tag('OutputFiles')}//{tag('File')}"):
        if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','FALSE')
    write_xml(tree,project);net=W['em'](project,2)
    (out/'screen.snp').write_bytes(project.with_suffix('.snp').read_bytes())
    (out/'screen_project.sonx').write_bytes(project.read_bytes())
    marker.write_text(json.dumps(dict(project_sha256=digest(project),snp_sha256=digest(project.with_suffix('.snp')),frequencies_GHz=list(frequencies)),indent=2)+'\n')
    return net,out
