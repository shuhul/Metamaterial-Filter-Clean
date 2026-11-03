"""Measure the assembled T2/T1 inductance matrix without changing their metal.

Two low-frequency points separate inductive and capacitive admittance terms.
A direct point at an actual passband peak checks that opening the temporary
measurement ports reproduces the ordinary physical simulation. Diagnostic
only: the primary fixed-L, free-capacitance S21 fit is not replaced.
"""
from copy import deepcopy
import json
import xml.etree.ElementTree as ET
import gdstk
import numpy as np
from geometry.prepare_sonnet import tag,points,write_xml,BIN,run as native_command
from geometry.taper.corrected_iteration import ROOT,OUT
from geometry.taper.measure_compact_inductance import AP,W
from geometry.taper.run_bridged_stage import digest


def prepare(name):
    source=ROOT/'sonnet_files/taper_exploration/compact_centered'/('corrected_'+name)/'idc/generated_ncell_3.sonx'
    meta=json.loads((ROOT/'geometry/taper/output/compact_centered'/('corrected_'+name)/'idc/geometry.json').read_text())
    reference=W['read_network'](source.with_suffix('.snp'),2)
    peak=int(np.argmax(abs(reference.s[:,1,0])))
    frequencies=[.2,.4,float(reference.f[peak]/1e9)]
    folder=source.parent.parent/'node_matrix_control';folder.mkdir(exist_ok=True)
    out=OUT/name/'node_matrix_control';out.mkdir(parents=True,exist_ok=True)
    project=folder/'T2_T1_node_probes.sonx'
    if not project.exists():
        tree=ET.parse(source);root=tree.getroot()
        original=[gdstk.Polygon(points(q)) for q in root.iter(tag('PlanarPolygon'))]
        prototype=next(ET.parse(source.parent.parent/'idc_Lcal/generated_ncell_1.sonx').iter(tag('Component')))
        maxid=max(int(q.get('Id')) for q in root.iter(tag('PlanarPolygon')))
        objects=next(q for q in root.find(tag('Geometry')).findall(tag('Level')) if q.get('Number')=='0').find(tag('Objects'))
        for j,x in enumerate(meta['cell_centers'][1:3]):
            if gdstk.inside([(x,447),(x,449),(x,455),(x,457)],original)!=(True,False,False,True):
                raise ValueError('Node probe is not on the intended ground and cap edges')
            comp=deepcopy(prototype);comp.set('MacroID',f'nodeMatrix{j}');comp.set('Name',f'NodeMatrix{j}')
            comp.set('Id',str(maxid+100+j));comp.find(tag('LabelPos')).text=f'({x+6},452)'
            comp.find(tag('SymbolBox')).text=f'({x-2},450,{x+2},454)'
            for i,(pin,y) in enumerate(zip(comp.findall(tag('Pin')),(448,456))):
                pin.set('MacroID',f'nodeMatrix{j}Pin{i}');pin.set('PortNum',str(3+2*j+i))
                pin.find(tag('Center')).text=f'({x},{y})'
            objects.append(comp)
        actual=[gdstk.Polygon(points(q)) for q in root.iter(tag('PlanarPolygon'))]
        if gdstk.boolean(actual,original,'xor',precision=.001):raise ValueError('Node probes changed metal')
        sweeps=root.find(tag('Sweeps'))
        for q in list(sweeps):sweeps.remove(q)
        fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
        for f in frequencies:
            ET.SubElement(fs,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(f),Stop=str(f),Step='1')
        for q in root.iter(tag('File')):
            if q.get('Type')=='TOUCH':q.set('Name',project.with_suffix('.snp').name);q.set('UseInterpolatedData','FALSE')
        write_xml(tree,project);native_command([BIN/'soncmd.exe','-Read',project],folder/'read.log')
    record=dict(source=str(source),source_sha256=digest(source),source_snp_sha256=digest(source.with_suffix('.snp')),
                project_sha256=digest(project),metal_unchanged=True,frequencies_GHz=frequencies,diagnostic_only=True)
    marker=out/'prepared.json'
    if marker.exists() and json.loads(marker.read_text())!=record:raise ValueError('Node probe sources changed')
    marker.write_text(json.dumps(record,indent=2)+'\n')
    return project,out,source


def run(name):
    project,out,source=prepare(name);marker=out/'operator_completed.json'
    if marker.exists():
        saved=json.loads(marker.read_text())
        if saved!=dict(project_sha256=digest(project),snp_sha256=digest(project.with_suffix('.snp'))):
            raise ValueError('Completed node operator changed')
        net=W['read_network'](project.with_suffix('.snp'),6)
    else:
        net=W['em'](project,6)
        marker.write_text(json.dumps(dict(project_sha256=digest(project),snp_sha256=digest(project.with_suffix('.snp'))),indent=2)+'\n')
    reduced,z0=AP['reduce_pairs'](net.s);y=AP['s_to_y'](reduced,z0)
    opened=AP['y_to_s'](AP['open_centers'](y),[50,50])
    reference=W['read_network'](source.with_suffix('.snp'),2)
    high=net.f>1e9;ix=np.array([np.argmin(abs(reference.f-f)) for f in net.f[high]])
    if not np.allclose(reference.f[ix],net.f[high],rtol=0,atol=1):raise ValueError('Peak-validation frequencies differ')
    error=abs(opened[high]-reference.s[ix]);passed=bool(np.max(error[:,1,0])<.05)
    if not passed:raise ValueError(f'Node-probe open validation failed: {error.max()}')
    low=net.f<1e9;k,c,res=AP['extract_kc'](net.f[low],y[low]);l=np.linalg.inv(k[2:,2:])
    result=dict(diagnostic_only=True,source=str(source),node_order=['Input','Output end','T2','T1'],
                open_peak_S21_max_error=float(error[:,1,0].max()),open_peak_all_S_max_error=float(error.max()),
                probe_validation_passed=passed,low_frequency_GHz=(net.f[low]/1e9).tolist(),
                K_inverse_nH=k.tolist(),C_matrix_fF=c.tolist(),internal_L_matrix_nH=l.tolist(),
                coupling_coefficient=float(l[0,1]/np.sqrt(l[0,0]*l[1,1])),
                C_T2_T1_matrix_fF=float(-c[2,3]),ground_row_sums_fF=c.sum(axis=1).tolist(),
                limitation='Low-frequency terminal-matrix diagnostic, not a replacement for the primary S21 fit or isolated L. Unprobed regular branch contributes to effective shunt capacitance.')
    (out/'node_matrix.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result
