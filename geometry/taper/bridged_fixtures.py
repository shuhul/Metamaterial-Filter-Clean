r"""Prepare feed-driven taper fixtures with upstream coupling bypassed.

Run: .venv\Scripts\python.exe -B -m geometry.taper.bridged_fixtures
Uses the standard GDS importer and native read check. Never runs EM.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from dataclasses import replace
from hashlib import sha256
from pathlib import Path
from types import SimpleNamespace
import json
import gdstk
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from geometry.prepare_sonnet import main as prepare_native
from geometry.taper.layout import build as build_t1
from geometry.taper.second_layout import build as build_t2,from_dict
from geometry.taper.locked_taper2 import T1,VERIFIED,require_native_locked

ROOT=Path(__file__).resolve().parents[2]
NAME='bridged_feed_tests'
NATIVE=ROOT/'sonnet_files/taper_exploration'/NAME
GDS_ROOT=ROOT/'geometry/taper/output'/NAME
WIRE_WIDTH=32  # Inspection candidate only; not a measured 0.74 nH result.
CASES=('01_t1_input_bridged','02_t1_input_restored',
       '03_t2_idc_removed_direct_feed','04_t2_idc_restored')


def t2_parameters():
    source=ROOT/'geometry/taper/output/response_match_span290_pairs31_rise32_right548_tab116.json'
    p=from_dict(json.loads(source.read_text())['second'])
    # Use a uniform centred wire to the original, flat bottom ground.
    m=replace(p.taper.meander,w=4*(3*WIRE_WIDTH//4+1),t=WIRE_WIDTH)
    return replace(p,taper=replace(p.taper,meander=m),taper1=T1,
                   x_meander=0,l_ground_rise=0,t_ground_tail=None,
                   straight_inductor=True)


def build(name,*,second=None,omit_t2_wire=False):
    if name not in CASES:raise ValueError('Unknown bridged-feed inspection case')
    t2=name.startswith(('03_','04_'))
    if not t2:
        if second is not None:raise ValueError('T2 parameters supplied for a T1 fixture')
        top,base=build_t1(T1,d_input=6,d_junction=56,input_style='finger',l_input=536,t_input=20)
    else:
        p=second if second is not None else t2_parameters()
        top,base=build_t2(p)
    if omit_t2_wire and name!=CASES[2]:raise ValueError('Wire removal is a direct-T2 diagnostic only')
    polygons=top.get_polygons()
    overlays=[]
    if name=='01_t1_input_bridged':
        # Remove the old upstream top bar and downward coupling finger.
        # Replace their entire footprint with the original 10 um feed width.
        input_left=base['cell_centers'][0]-T1.w/2
        cap_left=base['cell_centers'][1]-T1.w/2
        window=gdstk.rectangle((input_left,-8-base['l_input']),(cap_left,-8))
        polygons=[piece for q in polygons for piece in
                  gdstk.boolean([q],[window],'not',precision=.001,layer=1)]
        bridge=gdstk.rectangle((input_left-2,-18),(cap_left+2,-8),layer=1)
        polygons.append(bridge);overlays=[bridge]
    elif name=='03_t2_idc_removed_direct_feed':
        cap_left=base['cell_centers'][1]-p.taper.w/2
        input_left=cap_left-p.idc_span
        window=gdstk.rectangle((input_left,-8-p.taper.l_left),(cap_left,-8))
        # Remove both combs and the input trunk. Cap metal on the right of
        # cap_left is unchanged; no unused IDC teeth or trunks remain.
        polygons=[piece for q in polygons for piece in
                  gdstk.boolean([q],[window],'not',precision=.001,layer=1)]
        bridge=gdstk.rectangle((input_left-2,-18),(cap_left+2,-8),layer=1)
        polygons.append(bridge);overlays=[bridge]
    if omit_t2_wire:
        x=base['cell_centers'][1]+p.x_meander
        opening=gdstk.rectangle((x-p.taper.meander.t/2,-618),(x+p.taper.meander.t/2,-8-p.taper.t_top))
        polygons=[piece for q in polygons for piece in
                  gdstk.boolean([q],[opening],'not',precision=.001,layer=1)]
    titles={CASES[0]:'T1 Direct Feed — Upstream Electrode Removed',
            CASES[1]:'T1 Input Restored — Full Reference Fixture',
            CASES[2]:'T2 Direct Feed — IDC Removed, Downstream Chain Retained',
            CASES[3]:'T2 IDC Restored — Matching Integration Control'}
    meta=dict(base)
    meta.update(layout_kind='bridged_feed_inspection',case=name,
        source_layout_kind=base['layout_kind'],title=titles[name],
        expected_metal_components=len(gdstk.boolean(polygons,[],'or',precision=.001)),
        simulation_status='not_run_awaiting_geometry_review',
        upstream_coupling_bypassed=name in (CASES[0],CASES[2]),
        raised_ground_box=False,regular_cap_tabs=False,
        t2_wire_width_um=p.taper.meander.t if t2 else None,
        ports='Two normal 50-ohm box-wall feed ports; 400 um reference planes',
        overlay_bounds=[q.bounding_box() for q in overlays])
    if omit_t2_wire:
        meta['T2_wire_removed']=True
        meta['title']='T2 Direct Feed — Inductor Removed for Coupling Control'
    # Convert tuples to the same representation used in the JSON sidecar.
    return polygons,json.loads(json.dumps(meta)),overlays


def validate(gds,meta):
    second=from_dict(meta['second']) if meta['source_layout_kind']=='second_taper' else None
    expected,rebuilt,_=build(meta['case'],second=second,omit_t2_wire=meta.get('T2_wire_removed',False))
    if rebuilt!=meta:raise ValueError('Inspection geometry metadata changed')
    library=gdstk.read_gds(str(gds));tops=library.top_level()
    if len(tops)!=1 or tops[0].name!='TOP' or not np.isclose(library.unit,1e-6):
        raise ValueError('Invalid GDS library')
    actual=tops[0].get_polygons()
    if any(q.layer!=1 or q.datatype!=0 for q in actual):raise ValueError('Invalid metal layer')
    if gdstk.boolean(actual,expected,'xor',precision=.001):raise ValueError('GDS does not match intended fixture')
    p=meta['parameters']
    return SimpleNamespace(n_cells=p['n_cells'],l_sim=p['l_sim'],l_top_ground=p['l_top_ground']),actual


def preview(polygons,meta,overlays,path,ax=None):
    standalone=ax is None
    if standalone:fig,ax=plt.subplots(figsize=(15,6))
    for q in polygons:ax.add_patch(Polygon(q.points,facecolor='#388563',edgecolor='none'))
    for q in overlays:ax.add_patch(Polygon(q.points,facecolor='#1976d2',edgecolor='none'))
    if overlays:
        (a,b),(c,d)=overlays[0].bounding_box()
        ax.annotate('Direct metal connection',((a+c)/2,-13),((a+c)/2,45),ha='center',fontsize=10,
                    arrowprops=dict(arrowstyle='->',color='#1976d2'),color='#1256a0')
    for number,x,direction in ((1,0,1),(2,meta['w_sim'],-1)):
        ax.annotate('',xy=(x+direction*400,-13),xytext=(x,-13),
                    arrowprops=dict(arrowstyle='->',color='black'))
        ax.text(x+direction*80,30,f'Port {number}',ha='center',fontsize=10)
    ax.set(xlim=(-25,meta['w_sim']+25),ylim=(-665,100),xlabel='x (um)',ylabel='y (um)',title=meta['title'])
    ax.set_aspect('equal')
    if standalone:
        fig.tight_layout();fig.savefig(path,dpi=190);plt.close(fig)


def main():
    if NATIVE.exists() and any(NATIVE.rglob('*.snp')):
        raise ValueError('This folder contains simulated results; use a new NAME before regenerating')
    selection=ROOT/'geometry/taper/selected_design.json'
    selection_hash=sha256(selection.read_bytes()).hexdigest()
    prepared=[]
    for name in CASES:
        polygons,meta,overlays=build(name)
        gd=GDS_ROOT/name;gd.mkdir(parents=True,exist_ok=True)
        path=gd/'geometry.gds'
        lib=gdstk.Library(unit=1e-6,precision=1e-9);lib.new_cell('TOP').add(*polygons);lib.write_gds(str(path))
        path.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
        dest=NATIVE/name;dest.mkdir(parents=True,exist_ok=True)
        prepare_native(meta['parameters']['n_cells'],gds_path=path,output_dir=dest,geometry_validator=validate)
        project=dest/f"generated_ncell_{meta['parameters']['n_cells']}.sonx"
        # Both T2 files retain exact native T1/right-hand selected geometry.
        if name.startswith(('03_','04_')):require_native_locked(project,meta)
        preview(polygons,meta,overlays,dest/'geometry_review.png')
        prepared.append(dict(case=name,project=str(project),geometry=str(path),solver_run=False,
                             project_sha256=sha256(project.read_bytes()).hexdigest()))
    if sha256(selection.read_bytes()).hexdigest()!=selection_hash:raise ValueError('Selected T1 record changed')
    (NATIVE/'manifest.json').write_text(json.dumps(dict(status='inspection_only_no_simulations',cases=prepared),indent=2)+'\n')
    fig,axes=plt.subplots(4,1,figsize=(16,18))
    for name,ax in zip(CASES,axes):
        polygons,meta,overlays=build(name);preview(polygons,meta,overlays,None,ax)
    fig.suptitle('Feed-Driven Bridged Fixtures — Geometry Review Only',fontsize=20)
    fig.tight_layout(rect=(0,0,1,.97));fig.savefig(NATIVE/'all_fixtures.png',dpi=170);plt.close(fig)
    print(f'Prepared four feed-driven inspection projects. No EM simulations.\n{NATIVE}',flush=True)


if __name__=='__main__':main()
