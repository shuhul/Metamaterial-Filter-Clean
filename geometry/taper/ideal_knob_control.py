"""Three known capacitor loads to isolate useful T2 design directions.

No planar metal changes. These are diagnostic native projects, not generated
physical designs. Exporting the completed operator makes load optimization
independent of additional EM field solves.
"""
from copy import deepcopy
import json,shutil
import xml.etree.ElementTree as ET
from geometry.taper.compact_stages import ROOT,NATIVE,OUT,GDS
from geometry.prepare_sonnet import tag,write_xml,BIN,run as native_command
from geometry.taper.sonnet import add_capacitor,configure,simulate
from geometry.taper.run_bridged_stage import digest
from geometry.taper.ideal_load_search import export,optimize


def pair(tree,name,x0,x1,first,value_fF):
    root=tree.getroot();prototype=next(root.iter(tag('Component')));comp=deepcopy(prototype)
    ids=[int(q.get('Id')) for q in root.iter() if q.get('Id','').isdigit()]
    comp.set('Id',str(max(ids)+1));comp.set('MacroID',name);comp.set('Name',name)
    comp.find(tag('Value')).text=str(value_fF/1000)
    comp.find(tag('LabelPos')).text=f'({(x0+x1)/2},480)'
    comp.find(tag('SymbolBox')).text=f'({x0+2},458,{x1-2},462)'
    for i,(pin,x,orientation) in enumerate(zip(comp.findall(tag('Pin')),(x0,x1),('LEFT','RIGHT'))):
        pin.set('MacroID',f'{name}Pin{i}');pin.set('PortNum',str(first+i));pin.set('Orientation',orientation)
        pin.find(tag('Center')).text=f'({x},460)'
    level=next(q for q in root.find(tag('Geometry')).findall(tag('Level')) if q.get('Number')=='0')
    level.find(tag('Objects')).append(comp)


def run():
    name='plain_w500_wire52_idc30_span172_l504'
    source=NATIVE/name/'idc/generated_ncell_3.sonx';meta=json.loads((GDS/name/'idc/geometry.json').read_text())
    dest=ROOT/'sonnet_files/taper_exploration/ideal_load_diagnostics/plain_three_caps'
    out=ROOT/'output_data/taper_exploration/ideal_load_diagnostics/plain_three_caps'
    dest.mkdir(parents=True,exist_ok=True);out.mkdir(parents=True,exist_ok=True)
    project=dest/'diagnostic_three_caps.sonx';known=[900,5,10]
    if not project.exists():
        shutil.copyfile(source,project);before=ET.parse(source)
        add_capacitor(project,meta);tree=ET.parse(project)
        var=next(q for q in tree.iter(tag('Variable')) if q.get('Name')=='C_add');var.set('Value','0.9')
        p=meta['second'];right=meta['cell_centers'][1]+p['taper']['w']/2
        t1left=meta['cell_centers'][2]-p['taper1']['w']/2
        pair(tree,'extraSideC',right,t1left,5,known[1])
        capleft=meta['cell_centers'][1]-p['taper']['w']/2
        pair(tree,'extraInputC',capleft-p['idc_span']+20,capleft,7,known[2])
        polys=lambda r:[ET.tostring(q) for q in r.iter(tag('PlanarPolygon'))]
        if polys(before)!=polys(tree):raise ValueError('Planar geometry or external ports changed')
        write_xml(tree,project);configure(project,4.65,5.4,7500)
        native_command([BIN/'soncmd.exe','-Read',project],dest/'native_read.log')
        (out/'intent.json').write_text(json.dumps(dict(source=str(source),source_sha256=digest(source),
            diagnostic_only=True,physical_design=False,all_planar_polygons_and_external_ports_unchanged=True,
            known_caps_fF=known,capacitor_order=['T2 to ground','extra T2 to T1','extra input to T2']),indent=2)+'\n')
        (dest/'DIAGNOSTIC_ONLY.md').write_text('# Ideal load diagnostic only\n\nThe three capacitors are ideal Sonnet elements used to test design directions. No physical layout has been completed here. Every planar polygon and external port matches the plain source.\n')
    if any(q.get('Type')!='CAP' for q in ET.parse(project).iter(tag('Component'))):
        raise ValueError('This archived project is a ports-only crash diagnostic, not a runnable CAP sweep. Do not overwrite it.')
    marker=out/'run_status.json'
    if marker.exists() and json.loads(marker.read_text()).get('status')=='complete':
        r=json.loads(marker.read_text())
        if digest(project)!=r['project_sha256'] or digest(project.with_suffix('.snp'))!=r['snp_sha256']:
            raise ValueError('Completed diagnostic changed')
    else:simulate(project,out)
    n=export(project,out/'operator',known)
    return optimize(n,out/'operator')


if __name__=='__main__':run()
