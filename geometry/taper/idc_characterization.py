"""Independent two-terminal IDC check; no resonator inductors in this fixture."""
from pathlib import Path
from types import SimpleNamespace
from dataclasses import asdict
import json
import gdstk
import numpy as np
from geometry.taper.second_layout import from_dict
from geometry.prepare_sonnet import main as prepare
from geometry.taper.sonnet import configure,simulate
from circuit_common import read_snp

ROOT=Path(__file__).resolve().parents[2]


def build(p):
    t=p.taper;left=426+p.idc_span;center=left+t.w/2;right=left+t.w
    width=right+426;bottom=-618;parts=[]
    def rect(x0,y0,x1,y1):parts.append(gdstk.rectangle((x0,y0),(x1,y1)))
    rect(left,-8-t.t_top,right,-8)
    rect(left,-8-t.l_left,left+t.t_left,-8)
    rect(right-t.t_right,-8-t.l_right,right,-8)
    first=-8-t.t_top-p.d_idc_top;pitch=p.idc_t+p.idc_d
    last=first-(2*p.idc_pairs-1)*pitch-p.idc_t
    for row in range(2*p.idc_pairs):
        y=first-row*pitch
        if row%2:rect(446+p.idc_d,y-p.idc_t,left+t.t_left,y)
        else:rect(426,y-p.idc_t,left-p.idc_d,y)
    rect(426,last,446,-8);rect(0,-18,446,-8);rect(center,-18,width,-8)
    rect(0,0,width,448);rect(0,-1052,width,bottom)
    rect(0,-1052,412,-24);rect(width-412,-1052,width,-24)
    top=gdstk.Cell('TOP');top.add(*gdstk.boolean(parts,[],'or',precision=.001,layer=1))
    meta={'layout_kind':'idc_characterization','second':asdict(p),'w_sim':width,
          'expected_metal_components':4,'parameters':{'n_cells':1,'l_sim':1500,'l_top_ground':448,
          'w_side':412,'t_feed':10,'unit':{'d_top':8,'cap':{'t_top':t.t_top}}}}
    return top,meta


def validate(gds,meta):
    expected,rebuilt=build(from_dict(meta['second']))
    if rebuilt!=meta:raise ValueError('Stale IDC check metadata')
    lib=gdstk.read_gds(str(gds));polys=lib.top_level()[0].get_polygons()
    if gdstk.boolean(expected.get_polygons(),polys,'xor',precision=.001):raise ValueError('IDC GDS differs')
    return SimpleNamespace(n_cells=1,l_sim=1500,l_top_ground=448),polys


def run(name,p):
    name=name+'_idc_check';top,meta=build(p)
    gds=ROOT/'geometry/taper/output'/f'{name}.gds'
    lib=gdstk.Library(unit=1e-6,precision=1e-9);lib.add(top);lib.write_gds(str(gds))
    gds.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
    dest=ROOT/'sonnet_files/taper_exploration'/name
    output=ROOT/'output_data/taper_exploration'/name
    prepare(1,gds_path=gds,output_dir=dest,geometry_validator=validate)
    project=configure(dest/'generated_ncell_1.sonx',1,9,1600)
    snp=simulate(project,output);net=read_snp(snp)
    records=[]
    for ghz in (1,3,5,7,9):
        i=int(np.argmin(abs(net.f-ghz*1e9)));y=net.y[i];w=2*np.pi*net.f[i]
        records.append({'frequency_GHz':float(net.f[i]/1e9),'Cg_fF':float(-y[0,1].imag/w*1e15),
                        'C_input_ground_fF':float((y[0,0]+y[0,1]).imag/w*1e15),
                        'C_node_ground_fF':float((y[1,1]+y[0,1]).imag/w*1e15)})
    result={'purpose':'Independent IDC admittance check; changed downstream ground environment, not a final mixed-cell extraction',
            'parameters':asdict(p),'project':str(project),'snp':str(snp),'capacitances':records}
    (output/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result
