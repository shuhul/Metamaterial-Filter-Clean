"""First taper and mixed taper/regular layout; all dimensions are microns."""
from dataclasses import dataclass, asdict, field
from pathlib import Path
import json
from types import SimpleNamespace
import gdstk
import numpy as np
from geometry.meander import MeanderParams, make_meander
from geometry.straight_wire import StraightWireParams
from geometry.unit_cell import UnitCellParams, make_unit_cell
from geometry.placement import place

@dataclass(frozen=True)
class TaperParams:
    w: float=328
    t_top: float=14
    t_left: float=24
    t_right: float=6
    l_left: float=584
    l_right: float=600
    meander: MeanderParams=field(default_factory=lambda:MeanderParams(w=220))

    def __post_init__(self):
        values=[self.w,self.t_top,self.t_left,self.t_right,self.l_left,self.l_right]
        if not all(np.isfinite(v) and v>0 and v%2==0 for v in values):
            raise ValueError('Taper dimensions must be positive multiples of 2 um')
        if any(getattr(self.meander,k)%2 for k in ('w','l_top','l','t','d')):
            raise ValueError('Meander dimensions must fit the 2 um grid')
        if self.meander.w%4:raise ValueError('Use meander.w divisible by 4 for centered 4 um leads')
        if self.w%4:raise ValueError('Use taper.w divisible by 4 to keep center leads on the grid')
        if self.w/2-self.t_left<=self.meander.w/2 or self.w/2-self.t_right<=self.meander.w/2:
            raise ValueError('Cap legs must clear meander')
        if min(self.l_left,self.l_right)<=self.t_top or max(self.l_left,self.l_right)>=self.t_top+self.meander.l:
            raise ValueError('Cap legs must clear bottom ground')


def make_taper(p:TaperParams,name='TAPER',*,meander_offset_x=0,t_ground_tail=None,straight_inductor=False):
    """Origin centered 8 um above cap top; asymmetric legs, no ground pads."""
    if isinstance(p.meander,StraightWireParams) and not straight_inductor:
        raise ValueError('Straight-wire parameters cannot render a meander')
    parts=[gdstk.rectangle((-p.w/2,-8-p.t_top),(p.w/2,-8)),
           gdstk.rectangle((-p.w/2,-8-p.l_left),(-p.w/2+p.t_left,-8)),
           gdstk.rectangle((p.w/2-p.t_right,-8-p.l_right),(p.w/2,-8))]
    cap=gdstk.Cell(name+'_CAP');cap.add(*gdstk.boolean(parts,[],'or',precision=.001,layer=1))
    cell=gdstk.Cell(name);cell.add(gdstk.Reference(cap))
    if straight_inductor:
        # One thin vertical conductor. A local grounded pedestal determines
        # the exposed length; overlapping metal becomes part of the ground.
        cell.add(gdstk.rectangle((meander_offset_x-p.meander.t/2,-8-p.t_top-p.meander.l),
                                (meander_offset_x+p.meander.t/2,-8-p.t_top),layer=1))
    else:
        m=make_meander(p.meander,name=name+'_MEANDER')
        cell.add(place(m,origin=(meander_offset_x,-8-p.t_top),anchor='top'))
    if t_ground_tail is not None:
        # A wider grounded return below the last run reduces the long lead's L.
        bottom=-8-p.t_top-p.meander.l
        body_bottom=-8-p.t_top-p.meander.l_top-p.meander.l_body
        cell.add(gdstk.rectangle((meander_offset_x-t_ground_tail/2,bottom),
                                (meander_offset_x+t_ground_tail/2,body_bottom),layer=1))
    return cell


def build(p:TaperParams, *, isolation=False,d_input=6,d_junction=56,input_style='end_cell',l_input=None,t_input=24,
          w_top_slot=0,l_top_slot=100):
    """Mixed: taper end + taper resonator + regular resonator + regular end.

    No grounded strip in the strong input gap; strips remain at weak gaps.
    Isolation: three equal taper units, strong gaps, with middle resonator.
    All cap tops and the entire upper ground boundary stay horizontal.
    """
    if not all(np.isfinite(v) and v>=2 and v%2==0 for v in (d_input,d_junction)):
        raise ValueError('Coupling gaps must be positive multiples of 2 um')
    if d_junction<8:raise ValueError('Junction gap must leave clearance around its 4 um ground strip')
    if p.t_top<10:raise ValueError('End cap top must contain the 10 um feed')
    if input_style not in ('end_cell','finger'):raise ValueError('Unknown input coupler style')
    l_input=p.l_left if l_input is None else l_input
    if not all(np.isfinite(v) and v>0 and v%2==0 for v in (l_input,t_input)):
        raise ValueError('Input coupler dimensions must fit the 2 um grid')
    if l_input>=p.t_top+p.meander.l or t_input>=p.w:
        raise ValueError('Input finger must fit in its cap width and clear bottom ground')
    regular=UnitCellParams(meander=MeanderParams(w=186))
    taper=make_taper(p);reg=make_unit_cell(regular,name='REGULAR')
    widths=[p.w]*3 if isolation else [p.w,p.w,302,302]
    gaps=[d_input,d_input] if isolation else [d_input,d_junction,56]
    side=412;sep=14;top_height=448;height=1500;bottom=-8-p.t_top-p.meander.l
    if bottom!=-618:raise ValueError('Initial taper layouts must end at the regular bottom, y=-618 um')
    lefts=[side+sep]
    for w,g in zip(widths,gaps):lefts.append(lefts[-1]+w+g)
    centers=[x+w/2 for x,w in zip(lefts,widths)]
    width=lefts[-1]+widths[-1]+sep+side
    top=gdstk.Cell('TOP')
    for i,c in enumerate(centers):
        if i==0 and not isolation and input_style=='finger':
            # Match the native taper-unit fixture: feed bar and coupling finger, no inductor.
            top.add(gdstk.rectangle((c-p.w/2,-8-p.t_top),(c+p.w/2,-8),layer=1),
                    gdstk.rectangle((c+p.w/2-t_input,-8-l_input),(c+p.w/2,-8),layer=1))
        else:top.add(gdstk.Reference(taper if isolation or i<2 else reg,origin=(c,0)))
    parts=[gdstk.rectangle((0,0),(width,top_height)),
           gdstk.rectangle((0,top_height-height),(width,bottom)),
           gdstk.rectangle((0,top_height-height),(side,-24)),
           gdstk.rectangle((width-side,top_height-height),(width,-24))]
    if not isolation:
        for i in (1,2):
            x=lefts[i]+widths[i]+gaps[i]/2
            parts.append(gdstk.rectangle((x-2,bottom),(x+2,0)))
    metal=gdstk.boolean(parts,[],'or',precision=.001,layer=1)
    if w_top_slot:
        if w_top_slot%4 or not 4<=w_top_slot<=p.w-4 or l_top_slot%2 or not 2<=l_top_slot<=top_height-4:
            raise ValueError('Invalid upper ground slot')
        holes=[gdstk.rectangle((c-w_top_slot/2,2),(c+w_top_slot/2,2+l_top_slot)) for c in centers]
        metal=gdstk.boolean(metal,holes,'not',precision=.001,layer=1)
    ground=gdstk.Cell('GROUND');ground.add(*metal)
    top.add(gdstk.Reference(ground))
    top.add(gdstk.rectangle((0,-18),(centers[0],-8),layer=1),gdstk.rectangle((centers[-1],-18),(width,-8),layer=1))
    metadata={'expected_metal_components':2 if isolation or input_style=='finger' else 1,'layout_kind':'first_taper','w_sim':width,'cell_centers':centers,
              'input_style':input_style,'l_input':l_input,'t_input':t_input,'taper':asdict(p),'isolation':isolation,'d_input':d_input,'d_junction':d_junction,
              'parameters':{'n_cells':1 if isolation else 2,'l_sim':height,'l_top_ground':top_height,
                            'w_side':side,'t_feed':10,'unit':{'d_top':8,'cap':{'t_top':p.t_top}}},
              'regular_reference':{'meander_w':186,'d_cell':56,'L0_nH':4.922723945068658,'C0_fF':188.8781163861302}}
    metadata.update(w_top_slot=w_top_slot,l_top_slot=l_top_slot)
    return top,metadata


def write(p,basename,**kwargs):
    top,metadata=build(p,**kwargs);path=Path(basename).with_suffix('.gds');path.parent.mkdir(parents=True,exist_ok=True)
    lib=gdstk.Library(unit=1e-6,precision=1e-9);lib.add(top,*top.dependencies(True));lib.write_gds(str(path))
    path.with_suffix('.json').write_text(json.dumps(metadata,indent=2)+'\n')
    return path


def validate(gds,metadata):
    args=dict(metadata['taper']);args['meander']=MeanderParams(**args['meander'])
    top,expected=build(TaperParams(**args),isolation=metadata['isolation'],d_input=metadata['d_input'],d_junction=metadata['d_junction'],input_style=metadata.get('input_style','end_cell'),l_input=metadata.get('l_input'),t_input=metadata.get('t_input',24),w_top_slot=metadata.get('w_top_slot',0),l_top_slot=metadata.get('l_top_slot',100))
    for key in ('expected_metal_components','input_style','l_input','t_input','w_top_slot','l_top_slot'):
        if key not in metadata:expected.pop(key)
    if expected!=metadata:raise ValueError('Taper metadata does not match rebuilt dimensions')
    lib=gdstk.read_gds(str(gds));tops=lib.top_level()
    if len(tops)!=1 or tops[0].name!='TOP' or not np.isclose(lib.unit,1e-6,rtol=1e-9,atol=0):raise ValueError('Invalid taper GDS library')
    polys=tops[0].get_polygons()
    if any(p.layer!=1 or p.datatype!=0 for p in polys):raise ValueError('Unexpected layer/datatype')
    if gdstk.boolean(polys,top.get_polygons(),'xor',precision=.001):raise ValueError('Taper GDS does not match parameters')
    p=metadata['parameters'];return SimpleNamespace(n_cells=p['n_cells'],l_sim=p['l_sim'],l_top_ground=p['l_top_ground']),polys
