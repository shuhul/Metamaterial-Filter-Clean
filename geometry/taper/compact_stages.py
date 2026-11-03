"""Compact centred T2 tests: verify the directly driven wire before adding IDC.

Direct and isolated-L fixtures have NO IDC teeth, trunk, or unused node bus.
The selected T1 and all downstream metal are preserved in the direct fixture.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
from dataclasses import dataclass,asdict
from pathlib import Path
from types import SimpleNamespace
import gdstk
import numpy as np
from geometry.taper.second_layout import build as original_build,from_dict
from geometry.taper.bridged_fixtures import ROOT,preview,require_native_locked
from geometry.prepare_sonnet import main as prepare_native

NATIVE=ROOT/'sonnet_files/taper_exploration/compact_centered'
GDS=ROOT/'geometry/taper/output/compact_centered'
OUT=ROOT/'output_data/taper_exploration/compact_centered'


@dataclass(frozen=True)
class JunctionIDC:
    """Small 56.7 fF T2-to-T1 coupler, distinct from the 459 fF input IDC."""
    n: int=10
    t: int=2
    d: int=2
    l_finger: int=126
    d_top: int=10

    def __post_init__(self):
        if not isinstance(self.n,int) or self.n<2:raise ValueError('Junction pairs must be an integer >=2')
        if any(v<=0 or v%2 for v in (self.t,self.d,self.l_finger,self.d_top)):
            raise ValueError('Junction dimensions must be positive on the 2 um grid')


def build(p,stage='direct',*,omit_wire=False,bridge_junction=False,junction=None):
    if stage not in ('direct','Lcal','idc','idc_Lcal'):raise ValueError('Unknown compact stage')
    if p.taper.t_top!=10 or p.x_meander or not p.straight_inductor or p.l_ground_rise or p.t_ground_tail is not None:
        raise ValueError('Use a centred uniform straight wire, 10 um bar, and flat bottom ground')
    if omit_wire and stage not in ('direct','idc'):raise ValueError('Wire-removed control requires a response stage')
    top,meta=original_build(p,isolation=stage in ('Lcal','idc_Lcal'));polys=top.get_polygons()
    cap_left=meta['cell_centers'][1]-p.taper.w/2;input_left=cap_left-p.idc_span
    overlays=[]
    if stage in ('direct','Lcal'):
        # Clear both combs, input trunk, and the whole unused left node bus.
        opening=gdstk.rectangle((input_left,-8-p.taper.l_left),
                                (cap_left+p.taper.t_left,-18))
        polys=[piece for q in polys for piece in gdstk.boolean([q],[opening],'not',precision=.001,layer=1)]
        end=cap_left-8 if stage=='Lcal' else cap_left+2
        feed=gdstk.rectangle((input_left+18,-18),(end,-8),layer=1)
        polys.append(feed);overlays=[feed] if stage=='direct' else []
    if omit_wire:
        x=meta['cell_centers'][1];t=p.taper.meander.t
        opening=gdstk.rectangle((x-t/2,-618),(x+t/2,-18))
        polys=[piece for q in polys for piece in gdstk.boolean([q],[opening],'not',precision=.001,layer=1)]
    if bridge_junction:
        if stage!='direct':raise ValueError('Physical bridge is a direct-stage control only')
        right=cap_left+p.taper.w;left1=meta['cell_centers'][2]-p.taper1.w/2
        polys.append(gdstk.rectangle((right-2,-18),(left1+2,-8),layer=1))
    if junction is not None:
        if bridge_junction:raise ValueError('Do not combine a short with a junction capacitor')
        right=cap_left+p.taper.w;left1=meta['cell_centers'][2]-p.taper1.w/2
        gap=left1-right
        if not gap/2<junction.l_finger<=gap-2:raise ValueError('Junction fingers must overlap and clear the opposite bus')
        first=-18-junction.d_top
        last=first-(2*junction.n-1)*(junction.t+junction.d)-junction.t
        if last < -8-p.taper.l_right:raise ValueError('T2 junction bus is too short')
        for row in range(2*junction.n):
            y=first-row*(junction.t+junction.d)
            a,b=(right-2,right+junction.l_finger) if row%2==0 else (left1-junction.l_finger,left1+2)
            polys.append(gdstk.rectangle((a,y-junction.t),(b,y),layer=1))
    meta.update(layout_kind='compact_centered_t2',stage=stage,IDC_present=stage=='idc',
        unused_left_bus_removed=stage in ('direct','Lcal'),wire_removed=omit_wire,
        title='Compact T2: '+{'direct':'Direct Feed, No IDC','Lcal':'Isolated Wire Calibration','idc':'IDC Restored','idc_Lcal':'Wire Calibration with Node IDC Comb'}[stage],
        overlay_bounds=[q.bounding_box() for q in overlays],
        expected_metal_components=len(gdstk.boolean(polys,[],'or',precision=.001)))
    if bridge_junction:meta['junction_physical_bridge']=True
    if stage=='idc_Lcal':meta['node_IDC_comb_present']=True
    if junction is not None:
        meta.update(junction_idc=asdict(junction),input_IDC_present=stage=='idc',junction_IDC_present=True,IDC_present=True,
                    title='Compact T2: Short Junction Capacitor, '+('Input IDC Restored' if stage=='idc' else 'No Input IDC'))
    return polys,json.loads(json.dumps(meta)),overlays


def validate(gds,meta):
    junction=JunctionIDC(**meta['junction_idc']) if 'junction_idc' in meta else None
    polys,rebuilt,_=build(from_dict(meta['second']),meta['stage'],omit_wire=meta['wire_removed'],bridge_junction=meta.get('junction_physical_bridge',False),junction=junction)
    if rebuilt!=meta:raise ValueError('Stale compact-stage metadata')
    lib=gdstk.read_gds(str(gds));tops=lib.top_level()
    if len(tops)!=1 or tops[0].name!='TOP' or not np.isclose(lib.unit,1e-6,rtol=1e-9,atol=0):raise ValueError('Invalid GDS library')
    actual=tops[0].get_polygons()
    if any(q.layer!=1 or q.datatype!=0 for q in actual):raise ValueError('Invalid GDS layer')
    if gdstk.boolean(actual,polys,'xor',precision=.001):raise ValueError('GDS differs from parameters')
    p=meta['parameters']
    return SimpleNamespace(n_cells=p['n_cells'],l_sim=p['l_sim'],l_top_ground=p['l_top_ground']),actual


def prepare(name,p,stage='direct',*,omit_wire=False,bridge_junction=False,junction=None):
    if p.w_upper_plate:
        raise ValueError('Top-gap electrode variants were rejected by the user. Use the plain side-coupled design; archived geometry remains readable via build().')
    dest=NATIVE/name/stage;gd=GDS/name/stage;out=OUT/name/stage
    for folder in (dest,gd,out):folder.mkdir(parents=True,exist_ok=True)
    polys,meta,overlays=build(p,stage,omit_wire=omit_wire,bridge_junction=bridge_junction,junction=junction)
    project=dest/f"generated_ncell_{meta['parameters']['n_cells']}.sonx"
    if project.exists():
        if json.loads((gd/'geometry.json').read_text())!=meta:raise ValueError('Use a new name for changed geometry')
    else:
        lib=gdstk.Library(unit=1e-6,precision=1e-9);lib.new_cell('TOP').add(*polys)
        gds=gd/'geometry.gds';lib.write_gds(str(gds))
        gds.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
        prepare_native(meta['parameters']['n_cells'],gds_path=gds,output_dir=dest,geometry_validator=validate)
    if stage not in ('Lcal','idc_Lcal'):require_native_locked(project,meta)
    preview(polys,meta,overlays,dest/'geometry_review.png')
    return project,out,meta
