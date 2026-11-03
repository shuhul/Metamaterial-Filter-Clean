"""Sonnet preparation and recorded runs for taper exploration."""
from copy import deepcopy
from hashlib import sha256
from pathlib import Path
from datetime import datetime,timezone
import json
import subprocess
import os
import xml.etree.ElementTree as ET
import gdstk
import numpy as np
from geometry.taper.solver_guard import run as guarded_solver
from geometry.taper.solver_resources import solver_environment
from circuit_common import read_snp
from geometry.prepare_sonnet import tag,points,write_xml,BIN,BASE


def configure(project,start=4,stop=5.8,target=900):
    project=Path(project).resolve();tree=ET.parse(project);r=tree.getroot()
    for block in r.findall(f"{tag('Sweeps')}/{tag('Set')}"):
        if block.get('On')=='TRUE':
            for sweep in block.findall(f"{tag('Frequencies')}/{tag('Sweep')}"):
                if sweep.get('On')=='TRUE':sweep.set('Start',str(start));sweep.set('Stop',str(stop));sweep.set('Target',str(target))
    for e in r.findall(f"{tag('OutputFiles')}//{tag('File')}"):
        if e.get('Type')=='TOUCH':e.set('Name',project.with_suffix('.snp').name)
    write_xml(tree,project)
    return project


def simulate(project,output,*,cache_only=False):
    project=Path(project).resolve();output=Path(output).resolve();output.mkdir(parents=True,exist_ok=True)
    record={'started_utc':datetime.now(timezone.utc).isoformat(),'project':str(project),'project_sha256':sha256(project.read_bytes()).hexdigest(),'status':'running'}
    (output/'run_status.json').write_text(json.dumps(record,indent=2)+'\n')
    try:
        with (output/'em.log').open('w') as log:
            command=[str(BIN/'em.exe'),'-v']+(['-AbsNoDiscrete'] if cache_only else [])+[str(project)]
            r=guarded_solver(command,cwd=project.parent,stdout=log,env=solver_environment())
        r.check_returncode()
    except Exception as exc:
        record.update(status='failed',ended_utc=datetime.now(timezone.utc).isoformat(),
                      reason=str(exc),returncode=getattr(exc,'returncode',None))
        (output/'run_status.json').write_text(json.dumps(record,indent=2)+'\n')
        raise
    text=(output/'em.log').read_text()
    if 'Em errors detected: 0' not in text or 'Em simulation completed' not in text:raise RuntimeError('Incomplete Sonnet run')
    if cache_only and 'Total time for 0 frequencies:' not in text:raise RuntimeError('Cache-only check did not confirm zero new EM frequencies')
    snp=project.with_suffix('.snp')
    if not snp.exists():raise RuntimeError('Missing S-parameter export')
    net=read_snp(snp)
    root=ET.parse(project).getroot()
    active=[q for b in root.findall(f"{tag('Sweeps')}/{tag('Set')}") if b.get('On')=='TRUE' for q in b.findall(f"{tag('Frequencies')}/{tag('Sweep')}") if q.get('On')=='TRUE']
    if len(active)!=1:raise ValueError('Expected one active sweep')
    sweep=active[0]
    matches=np.allclose([net.f[0],net.f[-1]],[float(sweep.get('Start'))*1e9,float(sweep.get('Stop'))*1e9],rtol=0,atol=1)
    if not matches or len(net.f)<float(sweep.get('Target')):
        reason='Sonnet export does not match the requested range/grid; cache-only mode may have retained an older export'
        record.update(status='failed',reason=reason)
        (output/'run_status.json').write_text(json.dumps(record,indent=2)+'\n')
        raise RuntimeError(reason)
    (output/'response.snp').write_bytes(snp.read_bytes());(output/'project.sonx').write_bytes(project.read_bytes())
    record.update(status='complete',ended_utc=datetime.now(timezone.utc).isoformat(),snp_sha256=sha256(snp.read_bytes()).hexdigest())
    (output/'run_status.json').write_text(json.dumps(record,indent=2)+'\n')
    return snp


def add_capacitor(project,metadata):
    project=Path(project);tree=ET.parse(project);r=tree.getroot();geometry=r.find(tag('Geometry'))
    source=ET.parse(BASE/'sonnet_files/filter_v6.sonx')
    comp=deepcopy(next(c for c in source.iter(tag('Component')) if c.findtext(tag('Value'))=='C_add'))
    x=2*np.floor(metadata.get('inductance_pin_x',metadata['cell_centers'][1])/2+.5);y=448
    gap=int(metadata.get('unit',{}).get('d_top',8))
    if gap < 2 or gap % 2:
        raise ValueError('Calibration component needs a positive even ground gap')
    half=gap//2
    actual=[gdstk.Polygon(points(p)) for p in r.iter(tag('PlanarPolygon'))]
    if gdstk.inside([(x,y-1),(x,y+1),(x,y+gap-1),(x,y+gap+1)],actual)!=(True,False,False,True):raise ValueError('Invalid capacitor pin geometry')
    comp.set('Id',str(max(int(p.get('Id')) for p in r.iter(tag('PlanarPolygon')))+1));comp.set('MacroID','taperCadd')
    comp.find(tag('LabelPos')).text=f'({x+6},{y+half})';comp.find(tag('SymbolBox')).text=f'({x-2},{y+half-2},{x+2},{y+half+2})'
    for i,pin in enumerate(comp.findall(tag('Pin'))):
        pin.set('MacroID',f'taperPin{i}');pin.find(tag('Center')).text=f'({x},{y if i==0 else y+gap})'
    level=next(l for l in geometry.findall(tag('Level')) if l.get('Number')=='0');level.find(tag('Objects')).append(comp)
    var=next(v for v in geometry.findall(tag('Variable')) if v.get('Name')=='C_add');var.set('Value','0.0001')
    write_xml(tree,project);configure(project,3.5,5.8,2300)
    return project
