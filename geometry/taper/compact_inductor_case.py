"""T2 variant with a compact inductor below the input IDC.

All parameters are edited here; no command-line arguments are needed.
The capacitor frame, uniform IDC, and verified downstream cells are retained.
The long interior wire is replaced, not left as a second shunt path.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from dataclasses import asdict,replace
from pathlib import Path
from types import SimpleNamespace
import json
import gdstk
import numpy as np
from geometry.meander import MeanderParams,make_meander
from geometry.taper.second_layout import build,from_dict
from geometry.taper.locked_taper2 import require_locked,T1
from geometry.taper.response_cases import ROOT,case,A

CASE_NAME='response_compact_below_idc_w32_l90_parts'
RUN_EM=True
INDUCTOR=MeanderParams(w=32,l=90,l_top=8,n=3,t=4,d=4)
X_FROM_CAP_LEFT=-78
BRIDGE_T=4

def build_compact(p,m,x_from_left=-78,bridge_t=4,*,isolation=False):
    require_locked(p,allow_t2_junction=True)
    if not p.straight_inductor or p.l_ground_rise:raise ValueError('Start from a straight wire without a pedestal')
    if any(v%2 for v in (m.w,m.l,m.l_top,m.t,m.d,x_from_left,bridge_t)) or m.w%4:
        raise ValueError('Compact inductor must fit the 2 um grid')
    base,meta=build(p,isolation=isolation)
    cap_left=meta['cell_centers'][1]-p.taper.w/2
    original_x=meta['cell_centers'][1]+p.x_meander
    old_wire=gdstk.rectangle((original_x-p.taper.meander.t/2,-618),
                            (original_x+p.taper.meander.t/2,-8-p.taper.t_top))
    t2_cell=base.references[0].cell
    old_parts=[q.copy().translate(*base.references[0].origin) for q in t2_cell.polygons]
    if len(old_parts)!=1 or gdstk.boolean(old_parts,[old_wire],'xor',precision=.001):
        raise ValueError('Cannot identify the original independent T2 wire polygon')
    t2_cell.remove(*list(t2_cell.polygons))
    # Preserve component polygon boundaries. Unioning the entire grounded
    # network fractures every meander on import and hides the probe wire.
    metal=base.get_polygons()
    x=cap_left+x_from_left;top=-618+m.l
    # Keep the whole component below the lowest IDC finger, with 8 um clearance.
    lowest=-8-p.taper.t_top-p.d_idc_top-(2*p.idc_pairs-1)*(p.idc_t+p.idc_d)-p.idc_t
    if top+bridge_t>lowest-8:raise ValueError('Compact inductor bridge touches the IDC')
    if not cap_left-p.idc_span+24 < x-m.w/2 < x+m.w/2 < cap_left-4:
        raise ValueError('Compact inductor must fit below the IDC and clear both trunks')
    coil=make_meander(m,name='T2_COMPACT_INDUCTOR')
    new_polys=[q.copy().translate(x-m.w/2,-618) for q in coil.get_polygons()]
    bridge=gdstk.rectangle((x-m.t/2,top),(cap_left+p.taper.t_left,top+bridge_t),layer=1)
    metal=metal+new_polys+[bridge]
    if len(gdstk.boolean(metal,[],'or',precision=.001))!=meta['expected_metal_components']:
        raise ValueError('Unexpected compact-case connectivity')
    # Every downstream polygon must remain exactly the baseline geometry.
    if not isolation:
        left=meta['cell_centers'][2]-T1.w/2
        region=gdstk.rectangle((left,-1052),(meta['w_sim'],448))
        if gdstk.boolean(gdstk.boolean(metal,[region],'and'),
                         gdstk.boolean(base.get_polygons(),[region],'and'),'xor',precision=.001):
            raise ValueError('Compact inductor altered downstream metal')
    meta['layout_kind']='second_taper_compact_inductor'
    meta['external_inductor']={'meander':asdict(m),'x_from_cap_left':x_from_left,'t_bridge':bridge_t}
    meta['polygon_partition']='preserve_components_v2'
    meta['series_probe_T2_x']=x
    # Ideal calibration capacitors still attach to the cap bar at its center.
    meta['inductance_pin_x']=meta['cell_centers'][1]
    return metal,meta

def write_compact(name,p,m=INDUCTOR,x_from_left=X_FROM_CAP_LEFT,bridge_t=BRIDGE_T,*,isolation=False):
    metal,meta=build_compact(p,m,x_from_left,bridge_t,isolation=isolation)
    path=ROOT/'geometry/taper/output'/f'{name}.gds'
    if path.exists():
        if json.loads(path.with_suffix('.json').read_text())!=meta:raise ValueError('Changed geometry needs a new name')
    else:
        lib=gdstk.Library(unit=1e-6,precision=1e-9);top=lib.new_cell('TOP');top.add(*metal)
        lib.write_gds(str(path));path.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
    validate_compact(path,meta)
    out=ROOT/'output_data/taper_exploration'/name;out.mkdir(parents=True,exist_ok=True)
    plt=A['plt'];from matplotlib.patches import Polygon
    fig,ax=plt.subplots(figsize=(14,6))
    for q in metal:ax.add_patch(Polygon(q.points,facecolor='#388563',edgecolor='none'))
    ax.set(xlim=(350,meta['w_sim']-350),ylim=(-650,55),xlabel='x (um)',ylabel='y (um)',
        title='Isolated T2 Compact Inductor' if isolation else 'T2 Compact Inductor Below Uniform IDC; T1 Unchanged')
    ax.set_aspect('equal');fig.tight_layout();fig.savefig(out/'geometry_review.png',dpi=200);plt.close(fig)
    return path,out

def validate_compact(gds,meta):
    cfg=meta['external_inductor'];p=from_dict(meta['second'])
    expected,rebuilt=build_compact(p,MeanderParams(**cfg['meander']),cfg['x_from_cap_left'],cfg['t_bridge'],isolation=meta['isolation'])
    if meta!=rebuilt:raise ValueError('Stale compact-case metadata')
    lib=gdstk.read_gds(str(gds));tops=lib.top_level()
    if len(tops)!=1 or tops[0].name!='TOP' or not np.isclose(lib.unit,1e-6):raise ValueError('Invalid compact GDS library')
    actual=tops[0].get_polygons()
    if any(q.layer!=1 or q.datatype for q in actual):raise ValueError('Unexpected GDS layer/datatype')
    if gdstk.boolean(actual,expected,'xor',precision=.001):raise ValueError('Compact GDS differs from parameters')
    return SimpleNamespace(n_cells=meta['parameters']['n_cells'],l_sim=1500,l_top_ground=448),actual

def main():
    meta=json.loads((ROOT/'geometry/taper/output/response_w400_wireleft120_uniform30_gap6_pad24.json').read_text())
    p=replace(from_dict(meta['second']),l_ground_rise=0)
    gds,out=write_compact(CASE_NAME,p)
    if RUN_EM:case(CASE_NAME,p,prepared_gds=gds,geometry_validator=validate_compact)
    else:print(gds,out/'geometry_review.png')

if __name__=='__main__':main()
