"""Frank-style extraction fixtures ending at a directly driven resonator claw.

Only the output fixture changes. Retained T1 metal is copied exactly from the
reviewed GDS; no shield, ground pedestal, or extra electrode is introduced.
"""
import json
from types import SimpleNamespace
import gdstk
import numpy as np
from geometry.taper.compact_stages import ROOT
from geometry.taper.bridged_fixtures import preview
from geometry.prepare_sonnet import main as import_native
from geometry.taper.run_bridged_stage import digest
from geometry.taper.sonnet import configure,simulate

NAME='progressive_extraction_2026-10-23'
NATIVE=ROOT/'sonnet_files/taper_exploration'/NAME
GDS=ROOT/'geometry/taper/output'/NAME
OUT=ROOT/'output_data/taper_exploration'/NAME


def build(source,output_index,output_width,omit_wire=False):
    source=ROOT/source
    meta=json.loads(source.with_suffix('.json').read_text())
    original=gdstk.read_gds(str(source)).top_level()[0].get_polygons()
    if omit_wire:
        p=meta['second'];x=meta['cell_centers'][1];t=p['taper']['meander']['t']
        original=gdstk.boolean(original,[gdstk.rectangle((x-t/2,-618),(x+t/2,-18))],
                              'not',precision=.001,layer=1)
    center=meta['cell_centers'][output_index]
    edge=center+output_width/2;width=edge+426
    clip=gdstk.rectangle((0,-1052),(edge,448))
    polys=gdstk.boolean(original,[clip],'and',precision=.001,layer=1)
    polys.extend([gdstk.rectangle((edge,0),(width,448)),
                  gdstk.rectangle((edge,-1052),(width,-618)),
                  gdstk.rectangle((edge+14,-1052),(width,-24)),
                  gdstk.rectangle((center,-18),(width,-8))])
    polys=gdstk.boolean(polys,[],'or',precision=.001,layer=1)
    # All metal through the output claw's left edge must be exactly retained.
    roi=gdstk.rectangle((0,-1052),(center-output_width/2,448))
    if gdstk.boolean(gdstk.boolean(original,[roi],'and',precision=.001),
                     gdstk.boolean(polys,[roi],'and',precision=.001),'xor',precision=.001):
        raise ValueError('Retained upstream geometry changed')
    m=dict(layout_kind='progressive_claw_output',source_gds=str(source),source_sha256=digest(source),
           source_metadata_sha256=digest(source.with_suffix('.json')),output_index=output_index,
           output_width=output_width,w_sim=width,cell_centers=meta['cell_centers'][:output_index+1],
           expected_metal_components=len(polys),title='Progressive Extraction: Direct Output Claw',
           parameters=dict(n_cells=output_index,l_sim=1500,l_top_ground=448,w_side=412,
                           t_feed=10,unit=dict(d_top=8,cap=dict(t_top=14))),
           output_feed_direct=True,extra_output_end_cell=False,
           source_metadata=meta)
    if omit_wire:m['wire_removed']=True
    return polys,m


def validate(path,meta):
    expected,rebuilt=build(meta['source_gds'],meta['output_index'],meta['output_width'],meta.get('wire_removed',False))
    if rebuilt!=meta:raise ValueError('Progressive source changed')
    lib=gdstk.read_gds(str(path));tops=lib.top_level()
    if len(tops)!=1 or tops[0].name!='TOP' or not np.isclose(lib.unit,1e-6):raise ValueError('Invalid library')
    actual=tops[0].get_polygons()
    if gdstk.boolean(expected,actual,'xor',precision=.001):raise ValueError('GDS mismatch')
    return SimpleNamespace(n_cells=meta['parameters']['n_cells'],l_sim=1500,l_top_ground=448),actual


def prepare(name,source,index,width,omit_wire=False):
    polys,meta=build(source,index,width,omit_wire)
    dest=NATIVE/name;gd=GDS/name;out=OUT/name
    for p in (dest,gd,out):p.mkdir(parents=True,exist_ok=True)
    gds=gd/'geometry.gds';project=dest/f'generated_ncell_{index}.sonx'
    if not project.exists():
        lib=gdstk.Library(unit=1e-6,precision=1e-9);lib.new_cell('TOP').add(*polys);lib.write_gds(str(gds))
        gds.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
        import_native(index,gds_path=gds,output_dir=dest,geometry_validator=validate)
        configure(project,4.65,5.4,7500)
    elif json.loads(gds.with_suffix('.json').read_text())!=meta:raise ValueError('Changed existing fixture')
    preview(polys,meta,[],dest/'geometry_review.png')
    return project,out


def run(name,source,index,width,omit_wire=False):
    project,out=prepare(name,source,index,width,omit_wire)
    marker=out/'run_status.json'
    if marker.exists():
        status=json.loads(marker.read_text())
        if status['status']=='complete':
            if status['project_sha256']!=digest(project) or status['snp_sha256']!=digest(project.with_suffix('.snp')):
                raise ValueError('Changed completed project')
            return project.with_suffix('.snp'),out
    return simulate(project,out),out


CASES={
 '01_T1_regular':('geometry/taper/output/finger536_t20_right310_input6_j56.gds',2,302),
 '02_T2_direct_T1':('geometry/taper/output/compact_centered/w500_top10_wire52_gap6/direct/geometry.gds',2,328),
 '02a_T2_absent_T1':('geometry/taper/output/compact_centered/w500_top10_wire52_gap6/direct/geometry.gds',2,328,True),
 '03_T2_IDC_T1':('geometry/taper/output/compact_centered/plain_w500_wire52_idc26_r584_t24_gap6/idc/geometry.gds',2,328),
 '04_T2_IDC_T1_regular':('geometry/taper/output/compact_centered/plain_w500_wire52_idc26_r584_t24_gap6/idc/geometry.gds',3,302)}
CASES['05_T2_w700_direct_T1']=('geometry/taper/output/compact_centered/plain_w700_wire68_gap6/direct/geometry.gds',2,328)


if __name__=='__main__':
    import sys
    name=sys.argv[1]
    if '--prepare-only' in sys.argv:prepare(name,*CASES[name])
    else:run(name,*CASES[name])
