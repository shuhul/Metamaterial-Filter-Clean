"""Small low-frequency two-port IDC check; not a whole-cell C0 extraction.

Preserves finger count, overlap, gap, trunks and local substrate/metal/grid.
Different terminal surroundings make this a coupling sanity check, not a
replacement for the full T2/T1/regular fit.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from dataclasses import asdict,replace
from types import SimpleNamespace
from pathlib import Path
import json
import xml.etree.ElementTree as ET
import gdstk
import numpy as np
from geometry.taper.locked_taper2 import ROOT
from geometry.taper.second_layout import from_dict
from geometry.prepare_sonnet import main as prepare,tag,write_xml
from geometry.taper.sonnet import configure,simulate
from circuit_common import read_snp


def build(p):
    left=426;right=left+p.idc_span;width=right+p.taper.t_left+426
    first=-8-p.taper.t_top-p.d_idc_top;pitch=p.idc_t+p.idc_d
    bottom=first-(2*p.idc_pairs-1)*pitch-p.idc_t
    parts=[gdstk.rectangle((0,0),(width,448)),gdstk.rectangle((0,-1052),(width,-618)),
           gdstk.rectangle((0,-1052),(412,-24)),gdstk.rectangle((width-412,-1052),(width,-24)),
           gdstk.rectangle((left,bottom),(left+20,-8)),
           gdstk.rectangle((right,-8-p.taper.l_left),(right+p.taper.t_left,-8)),
           gdstk.rectangle((0,-18),(left+20,-8)),gdstk.rectangle((right,-18),(width,-8))]
    for row in range(2*p.idc_pairs):
        y=first-row*pitch
        x0,x1=(left+20+p.idc_d,right+p.taper.t_left) if row%2 else (left,right-p.idc_d)
        parts.append(gdstk.rectangle((x0,y-p.idc_t),(x1,y)))
    top=gdstk.Cell('TOP');top.add(*gdstk.boolean(parts,[],'or',precision=.001,layer=1))
    meta={'layout_kind':'idc_component_diagnostic','second':asdict(p),'w_sim':width,
          'expected_metal_components':4,'parameters':{'n_cells':1,'l_sim':1500,'l_top_ground':448,
          'w_side':412,'t_feed':10,'unit':{'d_top':8,'cap':{'t_top':14}}}}
    return top,meta


def validate(gds,meta):
    top,expected=build(from_dict(meta['second']))
    if expected!=meta:raise ValueError('IDC diagnostic metadata mismatch')
    lib=gdstk.read_gds(str(gds));tops=lib.top_level()
    if len(tops)!=1 or tops[0].name!='TOP':raise ValueError('Invalid IDC diagnostic cell')
    polygons=tops[0].get_polygons()
    if gdstk.boolean(polygons,top.get_polygons(),'xor',precision=.001):raise ValueError('IDC diagnostic metal mismatch')
    if len(gdstk.boolean(polygons,[],'or',precision=.001))!=4:raise ValueError('IDC diagnostic connectivity failed')
    return SimpleNamespace(n_cells=1,l_sim=1500,l_top_ground=448),polygons


def run(source_case,**changes):
    p=from_dict(json.loads((ROOT/'geometry/taper/output'/f'{source_case}.json').read_text())['second'])
    p=replace(p,**changes)
    name=f'idc_check_pairs{p.idc_pairs}_span{p.idc_span}_t{p.idc_t}_d{p.idc_d}'
    out=ROOT/'output_data/taper_exploration'/name;out.mkdir(parents=True,exist_ok=True)
    if (out/'idc_capacitance.json').exists():
        old_gds=ROOT/'geometry/taper/output'/f'{name}.gds'
        expected,_=build(p);actual=gdstk.read_gds(str(old_gds)).top_level()[0].get_polygons()
        if gdstk.boolean(expected.get_polygons(),actual,'xor',precision=.001):raise ValueError('Changed IDC diagnostic geometry needs a new case name')
        return json.loads((out/'idc_capacitance.json').read_text())
    gds=ROOT/'geometry/taper/output'/f'{name}.gds';top,meta=build(p)
    library=gdstk.Library(unit=1e-6,precision=1e-9);library.add(top);library.write_gds(str(gds))
    gds.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
    dest=ROOT/'sonnet_files/taper_exploration'/name
    prepare(1,gds_path=gds,output_dir=dest,geometry_validator=validate)
    project=dest/'generated_ncell_1.sonx';configure(project,.1,.3,2)
    tree=ET.parse(project)
    for block in tree.findall(f"{tag('Sweeps')}/{tag('Set')}"):
        if block.get('On')=='TRUE':
            for sweep in block.findall(f"{tag('Frequencies')}/{tag('Sweep')}"):
                if sweep.get('On')=='TRUE':sweep.set('Type','SWEEP');sweep.set('Adaptive','FALSE');sweep.set('Step','0.1')
    write_xml(tree,project)
    snp=simulate(project,out);net=read_snp(snp)
    identity=np.tile(np.eye(2,dtype=complex),(len(net.f),1,1))
    y=(identity-net.s)@np.linalg.inv(identity+net.s)/50
    coupling=-y[:,1,0].imag/(2*np.pi*net.f)*1e15
    result={'purpose':'IDC coupling sanity check only, terminal surroundings differ from full fixture',
            'frequency_GHz':(net.f/1e9).tolist(),'mutual_Cg_fF':coupling.tolist(),
            'mean_mutual_Cg_fF':float(np.mean(coupling)),'pairs':p.idc_pairs,
            'overlap_um':p.idc_span-20-2*p.idc_d,'finger_width_um':p.idc_t,'gap_um':p.idc_d,
            'project':str(project),'snp':str(snp)}
    (out/'idc_capacitance.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));return result


if __name__=='__main__':run('straight4_w400_xm76_idc30_span180_ground160_rise338_right536')
