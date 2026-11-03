"""Test whether isolated-L results transfer across different ground returns.

Keep full-chain ground, box, T2 metal and feed positions, but remove downstream
resonators. A floating dummy electrode retains a local capacitive neighbour.
The paired probe is across the existing 8 um top-ground gap. This diagnostic
does not certify an assembled circuit or change any production reference.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
from pathlib import Path
from types import SimpleNamespace
import xml.etree.ElementTree as ET
import gdstk
import numpy as np
from geometry.taper.compact_stages import build, ROOT
from geometry.taper.second_layout import from_dict
from geometry.taper.bridged_fixtures import preview
from geometry.prepare_sonnet import main as prepare_native, tag, write_xml
from geometry.taper.sonnet import add_capacitor
from geometry.taper.measure_compact_inductance import AP, W
from geometry.taper.run_bridged_stage import digest


def prepare(source_case, *, wire_width=None):
    source=ROOT/'geometry/taper/output/compact_centered'/source_case/'direct/geometry.json'
    p=from_dict(json.loads(source.read_text())['second'])
    if wire_width is not None:
        from dataclasses import replace
        p=replace(p,taper=replace(p.taper,meander=replace(p.taper.meander,t=wire_width,w=4*(3*wire_width//4+2))))
    full,meta,_=build(p)
    width=meta['w_sim'];left1=meta['cell_centers'][2]-p.taper1.w/2
    t2left=meta['cell_centers'][1]-p.taper.w/2
    # Explicit ground union from second_layout, identical to the full fixture.
    ground=[gdstk.rectangle((0,0),(width,448)),gdstk.rectangle((0,-1052),(width,-618)),
            gdstk.rectangle((0,-1052),(412,-24)),gdstk.rectangle((width-412,-1052),(width,-24))]
    for center in (left1+p.taper1.w+p.d_t1_regular/2,meta['cell_centers'][3]+151+28):
        ground.append(gdstk.rectangle((center-2,-618),(center+2,0)))
    mask=ground+[gdstk.rectangle((0,-618),(left1,-8)),
                 gdstk.rectangle((width-412,-18),(width,-8))]
    polys=gdstk.boolean(full,mask,'and',precision=.001,layer=1)
    # Weak input coupling so the capacitor-load pole can be measured.
    cut=gdstk.rectangle((t2left-8,-18),(t2left,-8))
    polys=gdstk.boolean(polys,[cut],'not',precision=.001,layer=1)
    polys.append(gdstk.rectangle((left1,-544),(left1+20,-8),layer=1))
    suffix='' if wire_width is None else f'_wire{wire_width}'
    name=source_case+suffix
    out=ROOT/'output_data/taper_exploration/ground_return_control'/name
    dest=ROOT/'sonnet_files/taper_exploration/ground_return_control'/name
    out.mkdir(parents=True,exist_ok=True);dest.mkdir(parents=True,exist_ok=True)
    # Prove the T2 region and all ground-only regions are preserved.
    t2window=gdstk.rectangle((t2left,-618),(t2left+p.taper.w,-8))
    if gdstk.boolean(gdstk.boolean(full,[t2window],'and'),gdstk.boolean(polys,[t2window],'and'),'xor',precision=.001):
        raise ValueError('T2 metal changed')
    if gdstk.boolean(ground,gdstk.boolean(polys,ground,'and'),'xor',precision=.001):
        raise ValueError('Ground metal changed')
    meta.update(title='Ground-Return Control: Full Ground, Isolated T2',
                expected_metal_components=len(gdstk.boolean(polys,[],'or',precision=.001)),
                ground_return_control=True,full_ground_preserved=True,T2_metal_preserved=True)
    gds=dest/'geometry.gds'
    if not gds.exists():
        lib=gdstk.Library(unit=1e-6,precision=1e-9);lib.new_cell('TOP').add(*polys);lib.write_gds(str(gds))
        gds.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
    elif json.loads(gds.with_suffix('.json').read_text())!=meta:
        raise ValueError('Changed control requires new name')
    def validate(path,record):
        actual=gdstk.read_gds(str(path)).top_level()[0].get_polygons()
        if record!=meta or gdstk.boolean(actual,polys,'xor',precision=.001):raise ValueError('Control GDS mismatch')
        return SimpleNamespace(n_cells=3,l_sim=1500,l_top_ground=448),actual
    project=dest/'generated_ncell_3.sonx'
    if not project.exists():
        prepare_native(3,gds_path=gds,output_dir=dest,geometry_validator=validate)
        add_capacitor(project,meta)
        tree=ET.parse(project);comp=next(tree.iter(tag('Component')));comp.set('Type','PORTS');comp.attrib.pop('Units',None)
        for q in list(comp):
            if q.tag==tag('Value'):comp.remove(q)
        sweeps=tree.find(tag('Sweeps'))
        for q in list(sweeps):sweeps.remove(q)
        fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
        for f in (.1,.3,1,4,5.5,7):
            ET.SubElement(fs,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(f),Stop=str(f),Step='1')
        for q in tree.findall(f"{tag('OutputFiles')}//{tag('File')}"):
            if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','FALSE')
        write_xml(tree,project)
    preview(polys,meta,[],dest/'geometry_review.png')
    return project,out


def measure(source_case, *, wire_width=None):
    project,out=prepare(source_case,wire_width=wire_width)
    marker=out/'operator_completed.json'
    if marker.exists():
        record=json.loads(marker.read_text())
        if record['project_sha256']!=digest(project) or record['snp_sha256']!=digest(project.with_suffix('.snp')):
            raise ValueError('Control source changed')
        net=W['read_network'](project.with_suffix('.snp'),4)
    else:
        net=W['em'](project,4)
        marker.write_text(json.dumps(dict(project_sha256=digest(project),snp_sha256=digest(project.with_suffix('.snp'))),indent=2))
    s,z=AP['reduce_pairs'](net.s);y=AP['s_to_y'](s,z)
    effective=y[:,2,2]-(y[:,2:3,:2]@np.linalg.solve(y[:,:2,:2]+np.eye(2)[None]/50,y[:,:2,2:3]))[:,0,0]
    # Apparent L at each frequency; capacitive contamination is visible in its
    # variation, and is fitted separately using low and microwave bands.
    rows=[]
    for indices in (net.f<=1e9,net.f>=4e9):
        w=2*np.pi*net.f[indices];basis=np.c_[-1/w,w]
        a=np.linalg.lstsq(basis*np.array([1e9,1e-12]),effective[indices].imag,rcond=None)[0]
        inverse_l=a[0]*1e9;c=a[1]*1e-12
        rows.append(dict(range_GHz=(net.f[indices][[0,-1]]/1e9).tolist(),L_nH=float(1e9/inverse_l),
                         effective_C_fF=float(c*1e15),max_Y_residual_S=float(np.max(abs(basis@np.array([inverse_l,c])-effective[indices].imag)))))
    result=dict(source=str(project),frequency_GHz=(net.f/1e9).tolist(),fits=rows,
                diagnostic_only=True,full_ground_preserved=True,assembly_certified=False)
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2),flush=True)
    return result


if __name__=='__main__':
    measure('w500_top10_wire52_gap6')
