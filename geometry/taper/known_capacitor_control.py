"""Optional known-capacitor perturbations of the assembled IDC fixture.

Keep the physical wire and every metal polygon. Add only a short-gap node
probe. Validate it open before using its numerical capacitor terminations.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from copy import deepcopy
import json
import shutil
import xml.etree.ElementTree as ET
import numpy as np
import gdstk
import skrf
from scipy.signal import find_peaks
from geometry.taper.clearance_extraction import ROOT,NATIVE,OUT,GDS
from geometry.prepare_sonnet import tag,points,write_xml,BIN,run as native_command
from geometry.taper.measure_compact_inductance import W
from geometry.taper.run_bridged_stage import digest


def prepare(case,*,ground_reference='AUTO'):
    if ground_reference not in ('AUTO','FLOAT'):raise ValueError('Unknown component ground reference')
    source=NATIVE/case/'idc/generated_ncell_3.sonx'
    prototype=NATIVE/case/'idc_Lcal/generated_ncell_1.sonx'
    name='known_capacitor_control'+('_floating' if ground_reference=='FLOAT' else '')
    dest=NATIVE/case/name;dest.mkdir(exist_ok=True)
    out=OUT/case/name;out.mkdir(exist_ok=True)
    project=dest/'load_probe.sonx'
    if not project.exists():
        tree=ET.parse(source);root=tree.getroot()
        comp=deepcopy(next(ET.parse(prototype).iter(tag('Component'))))
        if comp.get('Type')!='PORTS':raise ValueError('Expected two independent calibration terminals')
        comp.set('Id',str(max(int(q.get('Id')) for q in root.iter(tag('PlanarPolygon')))+100))
        comp.set('MacroID','knownCapProbe');comp.set('Name','KnownCapProbe')
        comp.find(tag('GndRef')).text=ground_reference
        for i,pin in enumerate(comp.findall(tag('Pin'))):pin.set('MacroID',f'knownCapPin{i}')
        level=next(q for q in root.find(tag('Geometry')).findall(tag('Level')) if q.get('Number')=='0')
        level.find(tag('Objects')).append(comp)
        for q in root.iter(tag('File')):
            if q.get('Type')=='TOUCH':q.set('Name',project.with_suffix('.snp').name)
        write_xml(tree,project);native_command([BIN/'soncmd.exe','-Read',project],dest/'read.log')
    metal=lambda p:[gdstk.Polygon(points(q)) for q in ET.parse(p).iter(tag('PlanarPolygon'))]
    if gdstk.boolean(metal(source),metal(project),'xor',precision=.001):raise ValueError('Probe changed metal')
    metadata=json.loads((GDS/case/'idc/geometry.json').read_text());x=metadata['cell_centers'][1]
    comp=next(ET.parse(project).iter(tag('Component')))
    if [q.findtext(tag('Center')) for q in comp.findall(tag('Pin'))]!=[f'({x},448)',f'({x},456)']:
        raise ValueError('Probe must occupy the original 8 um upper-ground gap')
    (out/'geometry_check.json').write_text(json.dumps(dict(source=str(source),source_sha256=digest(source),
        metal_XOR_empty=True,physical_wire_retained=True,probe_span_um=8,ground_reference=ground_reference),indent=2)+'\n')
    return source,project,out


def opened(s):
    return s[:,:2,:2]+s[:,:2,2:]@np.linalg.solve(np.eye(2)[None]-s[:,2:,2:],s[:,2:,:2])


def run(case='frame700_wire88_left214',*,ground_reference='AUTO'):
    source,project,out=prepare(case,ground_reference=ground_reference)
    record=json.loads((OUT/case/'idc/run_status.json').read_text())
    if (record.get('status')!='complete' or record['project_sha256']!=digest(source)
        or record['snp_sha256']!=digest(source.with_suffix('.snp'))):
        raise ValueError('Ordinary IDC reference is not verified complete')
    ordinary=W['read_network'](source.with_suffix('.snp'),2)
    pilot=out/'pilot_validation.json'
    if not pilot.exists():
        ix=find_peaks(abs(ordinary.s[:,1,0]),prominence=.025)[0]
        i=int(max(ix,key=lambda j:abs(ordinary.s[j,1,0])))
        f=ordinary.f[i]/1e9;tree=ET.parse(project);root=tree.getroot();sweeps=root.find(tag('Sweeps'))
        for q in list(sweeps):sweeps.remove(q)
        fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
        ET.SubElement(fs,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(f),Stop=str(f),Step='1')
        for q in root.iter(tag('File')):
            if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','FALSE')
        write_xml(tree,project);net=W['em'](project,4);delta=abs(opened(net.s)[0]-ordinary.s[i])
        check=dict(frequency_GHz=f,S21_error=float(delta[1,0]),all_S_max=float(delta.max()),
            passed=bool(delta[1,0]<.01 and delta.max()<.05),project_sha256=digest(project),snp_sha256=digest(project.with_suffix('.snp')))
        shutil.copy2(project,out/'pilot_project.sonx');shutil.copy2(project.with_suffix('.snp'),out/'pilot_data.snp')
        pilot.write_text(json.dumps(check,indent=2)+'\n')
    check=json.loads(pilot.read_text())
    if not check['passed']:raise ValueError(f'Node probe failed peak pilot; full solve not launched: {check}')
    marker=out/'full_completed.json'
    if not marker.exists():
        tree=ET.parse(project);root=tree.getroot();sweeps=root.find(tag('Sweeps'))
        for q in list(sweeps):sweeps.remove(q)
        fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
        ET.SubElement(fs,tag('Sweep'),Type='ABS_ENTRY',On='TRUE',Adaptive='TRUE',Start='4.65',Stop='5.4',Target='7500')
        for q in root.iter(tag('File')):
            if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','TRUE')
        write_xml(tree,project);net=W['em'](project,4)
        marker.write_text(json.dumps(dict(project_sha256=digest(project),snp_sha256=digest(project.with_suffix('.snp'))),indent=2)+'\n')
    else:
        if json.loads(marker.read_text())!=dict(project_sha256=digest(project),snp_sha256=digest(project.with_suffix('.snp'))):
            raise ValueError('Completed load operator changed')
        net=W['read_network'](project.with_suffix('.snp'),4)
    if not np.allclose(net.f,ordinary.f,rtol=0,atol=1):raise ValueError('Ordinary and probe grids differ')
    delta=abs(opened(net.s)-ordinary.s)
    validation=dict(S21_RMS=float(np.sqrt(np.mean(delta[:,1,0]**2))),S21_max=float(delta[:,1,0].max()),all_S_max=float(delta.max()))
    validation['passed']=validation['S21_RMS']<.005 and validation['S21_max']<.05
    (out/'open_validation.json').write_text(json.dumps(validation,indent=2)+'\n')
    if not validation['passed']:raise ValueError(f'Loaded data rejected because open validation failed: {validation}')
    records=[]
    for cap in (0,400,800,1200):
        s=opened(net.s) if not cap else W['terminate_capacitor'](net.s,net.f,cap/1000)
        path=out/f'load_{cap:04d}fF.s2p'
        skrf.Network(frequency=skrf.Frequency.from_f(net.f,unit='hz'),s=s,z0=50).write_touchstone(str(path))
        records.append(dict(known_C_fF=cap,path=str(path),sha256=digest(path)))
    (out/'loads.json').write_text(json.dumps(records,indent=2)+'\n')
    return records


if __name__=='__main__':
    import sys
    prepare('frame700_wire88_left214') if '--prepare' in sys.argv else run()
