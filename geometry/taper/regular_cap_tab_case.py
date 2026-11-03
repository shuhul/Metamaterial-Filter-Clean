"""Separate response diagnostic: add a small signal-cap tab to regular cells.

T1 is kept exactly unchanged. Original unit-cell files are not overwritten.
The upper ground boundary stays continuous. Tabs are 2 um tall and chamfered.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from dataclasses import asdict,dataclass
from types import SimpleNamespace
import json
import gdstk
import numpy as np
from geometry.taper.second_layout import build,from_dict
from geometry.taper.locked_taper2 import require_locked,require_native_locked,T1
from geometry.taper.response_cases import ROOT,case,A

@dataclass(frozen=True)
class RegularCapTab:
    w:float=100
    t:float=2

def build_tabbed(p,tab=RegularCapTab()):
    require_locked(p,allow_t2_junction=True)
    if tab.w<0 or tab.w>280 or tab.w%4 or tab.t not in (2,4):raise ValueError('Tabs require width 0..280 in 4 um steps, height 2 or 4 um')
    base,meta=build(p);parts=base.get_polygons();regions=[]
    for center in meta['cell_centers'][3:]:
        x=2*np.floor(center/2+.5)
        if tab.w:
            q=gdstk.Polygon([(x-tab.w/2-tab.t,-8),(x-tab.w/2,-8+tab.t),
                              (x+tab.w/2,-8+tab.t),(x+tab.w/2+tab.t,-8)],layer=1)
            parts.append(q)
            regions.append(gdstk.rectangle((x-tab.w/2-tab.t,456-tab.t),(x+tab.w/2+tab.t,456)))
    if len(gdstk.boolean(parts,[],'or',precision=.001))!=2:raise ValueError('Unexpected tabbed-fixture connectivity')
    meta['layout_kind']='second_taper_regular_cap_tab';meta['regular_cap_tab']=asdict(tab)
    return parts,meta,regions

def validate_tabbed(gds,meta):
    expected,rebuilt,_=build_tabbed(from_dict(meta['second']),RegularCapTab(**meta['regular_cap_tab']))
    if rebuilt!=meta:raise ValueError('Stale tabbed-fixture metadata')
    lib=gdstk.read_gds(str(gds));top=lib.top_level()
    if len(top)!=1 or top[0].name!='TOP' or not np.isclose(lib.unit,1e-6):raise ValueError('Invalid library')
    actual=top[0].get_polygons()
    if any(q.layer!=1 or q.datatype for q in actual):raise ValueError('Invalid layer/datatype')
    if gdstk.boolean(actual,expected,'xor',precision=.001):raise ValueError('GDS differs from tabbed specification')
    return SimpleNamespace(n_cells=3,l_sim=1500,l_top_ground=448),actual

def validate_native_tabbed(project,meta):
    _,_,regions=build_tabbed(from_dict(meta['second']),RegularCapTab(**meta['regular_cap_tab']))
    # Only the explicitly generated regular-cap tabs may differ. This mask
    # cannot reach T1, its junction, the regular meander, or ground.
    t1right=meta['cell_centers'][2]+T1.w/2
    if any(q.bounding_box()[0][0]<=t1right for q in regions):raise ValueError('Tab mask overlaps T1')
    return require_native_locked(project,meta,allowed_regions=regions)

def write_tabbed(name,p,tab=RegularCapTab()):
    parts,meta,_=build_tabbed(p,tab);path=ROOT/'geometry/taper/output'/f'{name}.gds'
    if path.exists():
        if json.loads(path.with_suffix('.json').read_text())!=meta:raise ValueError('Changed geometry needs a new name')
    else:
        lib=gdstk.Library(unit=1e-6,precision=1e-9);top=lib.new_cell('TOP');top.add(*parts);lib.write_gds(str(path))
        path.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
    validate_tabbed(path,meta)
    out=ROOT/'output_data/taper_exploration'/name;out.mkdir(parents=True,exist_ok=True)
    plt=A['plt'];from matplotlib.patches import Polygon
    fig,ax=plt.subplots(figsize=(14,6))
    for q in parts:ax.add_patch(Polygon(q.points,facecolor='#388563',edgecolor='none'))
    for center in meta['cell_centers'][3:]:
        ax.annotate(f'{tab.w:g} x {tab.t:g} um cap tab',xy=(center,-8+tab.t),xytext=(center,30),ha='center',fontsize=9,
                    arrowprops={'arrowstyle':'->'})
    ax.set(xlim=(350,meta['w_sim']-350),ylim=(-650,55),xlabel='x (um)',ylabel='y (um)',title='Separate Regular-Cap Correction Test; T1 Unchanged')
    ax.set_aspect('equal');fig.tight_layout();fig.savefig(out/'geometry_review.png',dpi=200);plt.close(fig)
    return path,out

def run(name,p,tab=RegularCapTab(),*,frequencies=(4.5,4.9,5.,5.5),analyze_screen=True):
    gds,_=write_tabbed(name,p,tab)
    return case(name,p,prepared_gds=gds,geometry_validator=validate_tabbed,native_validator=validate_native_tabbed,
                frequencies=frequencies,analyze_screen=analyze_screen)

if __name__=='__main__':
    p=from_dict(json.loads((ROOT/'geometry/taper/output/response_w400_xm160_span206_rise0_right578.json').read_text())['second'])
    run('response_best_regular_tab100',p,frequencies=(5.,),analyze_screen=False)
