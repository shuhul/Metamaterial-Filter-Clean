"""T2 with vertical IDC fingers connected at the common upper LC node."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from dataclasses import dataclass,asdict,replace
from types import SimpleNamespace
import json
import gdstk
import numpy as np
from geometry.taper.second_layout import build,from_dict
from geometry.taper.locked_taper2 import require_locked,T1
from geometry.taper.response_cases import ROOT,A,case

@dataclass(frozen=True)
class VerticalIDC:
    n:int=16             # TOTAL alternating fingers; n/2 fingers per electrode.
    l:float=480          # Vertical separation between the two bus top edges.
    t:float=4
    d:float=4
    t_bus:float=4        # Same width as fingers: no wide plate over the IDC.

CASE_NAME='response_vertical16_l480_w328_rise358_v2'
IDC=VerticalIDC()
RUN_EM=True

def build_vertical(p,idc=IDC,*,isolation=False):
    require_locked(p,allow_t2_junction=True)
    if not p.straight_inductor:raise ValueError('This variant uses a separate straight shunt wire')
    if idc.n<4 or idc.n%2 or any(v<=0 or v%2 for v in (idc.l,idc.t,idc.d,idc.t_bus)):
        raise ValueError('IDC needs an even finger count and positive 2 um grid dimensions')
    base,meta=build(p,isolation=isolation);left=meta['cell_centers'][1]-p.taper.w/2
    input_left=left-p.idc_span;start=input_left+20+idc.d
    finish=start+(idc.n-1)*(idc.t+idc.d)+idc.t
    if finish>left-4:raise ValueError('Vertical fingers must clear the capacitor-frame edge')
    bus_top=-8;lower_top=bus_top-idc.l;bottom=lower_top-idc.t_bus
    if bottom<=-608 or idc.l<40:raise ValueError('IDC must clear the lower ground')
    wire=meta['cell_centers'][1]+p.x_meander
    if p.l_ground_rise and wire-p.w_ground_riser/2<left+p.taper.t_left:
        raise ValueError('The ground return must stay outside the region being replaced')
    # Remove the old horizontal comb and its tall signal bus/left cap leg.
    # Keep component polygons separate for accurate, efficient GAP probing.
    erase=gdstk.rectangle((input_left,-8-p.taper.l_left),(left+p.taper.t_left,-8))
    parts=[]
    for poly in base.get_polygons():parts.extend(gdstk.boolean([poly],[erase],'not',precision=.001,layer=1))
    def rect(a,b,c,d):parts.append(gdstk.rectangle((a,b),(c,d),layer=1))
    # The retained original feed already reaches input_left; the trunk below
    # joins its right edge. Do not duplicate the port edge at x=0.
    rect(input_left,bottom,input_left+20,-8)
    rect(input_left,bottom,finish,lower_top)
    # The node bus joins the capacitor bar directly at the top LC node.
    rect(start,bus_top-idc.t_bus,left+p.taper.t_left,bus_top)
    for row in range(idc.n):
        x=start+row*(idc.t+idc.d)
        if row%2==0:rect(x,lower_top+idc.d,x+idc.t,bus_top)
        elif not isolation:rect(x,bottom,x+idc.t,bus_top-idc.t_bus-idc.d)
    if len(gdstk.boolean(parts,[],'or',precision=.001))!=meta['expected_metal_components']:
        raise ValueError('Unexpected vertical-IDC connectivity')
    if not isolation:
        t1left=meta['cell_centers'][2]-T1.w/2
        region=gdstk.rectangle((t1left,-1052),(meta['w_sim'],448))
        if gdstk.boolean(gdstk.boolean(parts,[region],'and'),gdstk.boolean(base.get_polygons(),[region],'and'),'xor',precision=.001):
            raise ValueError('Vertical IDC changed downstream geometry')
    meta['layout_kind']='second_taper_vertical_idc';meta['vertical_idc']=asdict(idc)
    meta['polygon_partition']='vertical_ports_v2'
    meta['series_probe_T2_x']=wire;meta['inductance_pin_x']=wire
    return parts,meta

def validate_vertical(gds,meta):
    parts,expected=build_vertical(from_dict(meta['second']),VerticalIDC(**meta['vertical_idc']),isolation=meta['isolation'])
    if meta!=expected:raise ValueError('Stale vertical-IDC metadata')
    lib=gdstk.read_gds(str(gds));tops=lib.top_level()
    if len(tops)!=1 or tops[0].name!='TOP' or not np.isclose(lib.unit,1e-6):raise ValueError('Invalid library')
    actual=tops[0].get_polygons()
    if any(q.layer!=1 or q.datatype for q in actual):raise ValueError('Unexpected layer/datatype')
    if gdstk.boolean(parts,actual,'xor',precision=.001):raise ValueError('GDS differs from vertical-IDC specification')
    return SimpleNamespace(n_cells=meta['parameters']['n_cells'],l_sim=1500,l_top_ground=448),actual

def write_vertical(name,p,idc=IDC,*,isolation=False):
    parts,meta=build_vertical(p,idc,isolation=isolation)
    gds=ROOT/'geometry/taper/output'/f'{name}.gds'
    if gds.exists():
        if json.loads(gds.with_suffix('.json').read_text())!=meta:raise ValueError('Changed geometry requires a new case name')
    else:
        lib=gdstk.Library(unit=1e-6,precision=1e-9);top=lib.new_cell('TOP');top.add(*parts)
        lib.write_gds(str(gds));gds.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
    validate_vertical(gds,meta)
    out=ROOT/'output_data/taper_exploration'/name;out.mkdir(parents=True,exist_ok=True)
    plt=A['plt'];from matplotlib.patches import Polygon,Rectangle
    fig,ax=plt.subplots(figsize=(14,6))
    for q in parts:ax.add_patch(Polygon(q.points,facecolor='#388563',edgecolor='none'))
    if p.l_ground_rise:
        x=meta['series_probe_T2_x'];ax.add_patch(Rectangle((x-p.w_ground_riser/2,-618),p.w_ground_riser,p.l_ground_rise,fc='#b7bdc3'))
        ax.text(x,-618+p.l_ground_rise/2,'GND',ha='center',va='center',rotation=90,fontsize=9)
    ax.set(xlim=(350,meta['w_sim']-350),ylim=(-650,55),xlabel='x (um)',ylabel='y (um)',
        title='T2 Vertical IDC: Common Upper LC Node; T1 Unchanged' if not isolation else 'Isolated T2 Vertical IDC')
    ax.set_aspect('equal');fig.tight_layout();fig.savefig(out/'geometry_review.png',dpi=200);plt.close(fig)
    return gds,out

def parameters():
    p=from_dict(json.loads((ROOT/'geometry/taper/output/response_w400_wireleft120_uniform30_gap6_pad24.json').read_text())['second'])
    return replace(p,x_meander=-84,w_ground_riser=80,l_ground_rise=358,
                   taper=replace(p.taper,w=328))

def main():
    p=parameters();gds,out=write_vertical(CASE_NAME,p)
    if RUN_EM:case(CASE_NAME,p,prepared_gds=gds,geometry_validator=validate_vertical)
    else:print(gds,out/'geometry_review.png')

if __name__=='__main__':main()
