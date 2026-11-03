"""Taper 2 + taper 1 + regular fixture with an explicit two-terminal IDC."""
from dataclasses import dataclass, asdict, field
from pathlib import Path
from types import SimpleNamespace
import json
import gdstk
import numpy as np
from geometry.meander import MeanderParams
from geometry.straight_wire import StraightWireParams
from geometry.taper.layout import TaperParams, make_taper
from geometry.unit_cell import UnitCellParams, make_unit_cell


@dataclass(frozen=True)
class SecondParams:
    taper: TaperParams = field(default_factory=lambda: TaperParams(
        w=328,t_top=14,t_left=24,t_right=6,l_left=584,l_right=310,
        meander=MeanderParams(w=32,l=596,l_top=10,n=58,t=4,d=6)))
    idc_pairs: int = 20
    idc_t: int = 4
    idc_d: int = 4
    idc_span: int = 328  # Trunk-to-taper distance; overlap = span - 20 - 2*d.
    d_idc_top: int = 10  # Clearance below the underside of the cap top bar.
    equal_idc_sides: bool = False  # Extend the source bus to the node bus bottom.
    x_meander: int = 0
    t_ground_tail: int | None = None
    t_coupling_finger: int = 0  # Optional t2-node finger inside t1's left gap.
    l_coupling_finger: int = 566
    taper1: TaperParams | None = None  # Optional embedded variant; original remains the default.
    n_coupling_teeth: int = 0
    w_coupling_teeth: int = 100
    w_upper_plate: int = 0  # Driven t2 strip in the 8 um gap over t1's cap bar.
    l_coupling_bridge: int = 606
    l_cal_output: int = 536
    w_t1_ground_slot: int = 0
    l_t1_ground_slot: int = 100
    d_t2_t1: int = 6
    l_t1_left: int = 584
    t_t1_left: int = 24
    d_t1_regular: int = 56
    t_junction_ground: int = 4  # Only the T1-to-regular grounded strip.
    x_junction_ground: int = 0  # Offset of that existing strip; positive toward regular cell.
    l_ground_rise: int = 0  # Optional local ground pedestal under T2 only.
    w_ground_riser: int = 40
    straight_inductor: bool = False  # T2 only; T1 always keeps its verified meander.
    w_top_left_extension: int = 0  # Extend T2 top over the IDC without widening the box.
    t_inner_ground: int = 0  # Optional narrow lateral shield; never under the wire.
    x_inner_ground: int = 0  # Offset from the T2 frame centre.
    d_inner_ground_top: int = 4  # Shield clearance below the top bar.
    w_inner_ground_hole: int = 0  # Optional opening, leaving a grounded shield loop.
    w_ground_pad: int = 0  # Optional identical bottom-ground pads, one inside each cap leg.
    l_ground_pad: int = 586
    d_ground_pad: int = 12  # Clearance from inner cap-leg edge to pad outer edge.

    def __post_init__(self):
        if 8+self.taper.t_top+self.taper.meander.l!=618:
            raise ValueError('T2 meander must meet the shared bottom ground: 8 + t_top + meander.l = 618')
        if not isinstance(self.idc_pairs,int) or self.idc_pairs<1:
            raise ValueError('IDC pairs must be a positive integer')
        for key in ('idc_t','idc_d','idc_span','d_idc_top','d_t2_t1','l_t1_left','t_t1_left','d_t1_regular'):
            value=getattr(self,key)
            if value<=0 or value%2:raise ValueError(f'{key} must be positive and on the 2 um grid')
        if self.idc_span<40:raise ValueError('IDC span too short')
        if self.idc_span<=20+2*self.idc_d:raise ValueError('IDC fingers need positive overlap')
        if self.w_top_left_extension<0 or self.w_top_left_extension%2:
            raise ValueError('Top extension must be nonnegative and on the 2 um grid')
        if self.w_top_left_extension and self.w_top_left_extension>self.idc_span-20-self.idc_d:
            raise ValueError('Top extension must clear the input trunk')
        if self.d_t1_regular<8:raise ValueError('Weak gap must clear its grounded strip')
        if self.t_junction_ground not in (2,4):raise ValueError('Junction ground thickness must be 2 or 4 um')
        if (not np.isfinite(self.x_junction_ground) or self.x_junction_ground%2 or
            abs(self.x_junction_ground)>(self.d_t1_regular-self.t_junction_ground)/2-2):
            raise ValueError('Junction ground offset must be on grid and clear both cap legs by at least 2 um')
        if self.l_ground_rise<0 or self.l_ground_rise%2 or self.l_ground_rise>=self.taper.meander.l_bot:
            raise ValueError('Ground pedestal must be on grid and stay below the meander body')
        if self.w_ground_riser%4 or (self.l_ground_rise and self.w_ground_riser<self.taper.meander.t):
            raise ValueError('Ground pedestal width must cover the trace and be a multiple of 4 um')
        if self.l_ground_rise and (self.x_meander-self.w_ground_riser/2<=-self.taper.w/2+self.taper.t_left or self.x_meander+self.w_ground_riser/2>=self.taper.w/2-self.taper.t_right):
            raise ValueError('Ground pedestal must clear both cap legs')
        if self.x_meander%2:raise ValueError('Meander offset must be on the 2 um grid')
        m=self.taper.meander
        if isinstance(m,StraightWireParams) and not self.straight_inductor:
            raise ValueError('StraightWireParams requires straight_inductor=True')
        if self.t_ground_tail is not None and (self.t_ground_tail%4 or not m.t<=self.t_ground_tail<=m.w):
            raise ValueError('Ground tail must fit inside meander width, be at least trace thickness and a multiple of 4 um')
        # A straight inductor has no horizontal meander runs. Its actual
        # footprint is the trace (or explicit tail), not the unused coil w.
        if self.straight_inductor and m.t%4:
            raise ValueError('Centered straight trace width must be a multiple of 4 um')
        inductor_width=max(m.t,self.t_ground_tail or 0) if self.straight_inductor else m.w
        if self.w_ground_pad:
            if any(v <= 0 or v % 2 for v in (self.w_ground_pad, self.l_ground_pad, self.d_ground_pad)):
                raise ValueError('Ground pads must use positive dimensions on the 2 um grid')
            if self.l_ground_pad >= m.l:
                raise ValueError('Ground pads must clear the underside of the top bar')
            inner = self.taper.w/2-max(self.taper.t_left,self.taper.t_right)-self.d_ground_pad-self.w_ground_pad
            if inner <= abs(self.x_meander)+inductor_width/2:
                raise ValueError('Ground pads must clear the full inductor')
        if self.t_inner_ground:
            if self.t_inner_ground<4 or self.t_inner_ground%4 or self.x_inner_ground%2 or self.d_inner_ground_top<2 or self.d_inner_ground_top%2:
                raise ValueError('Inner ground shield must clear the cap and use the 2 um grid')
            a=self.x_inner_ground-self.t_inner_ground/2
            b=self.x_inner_ground+self.t_inner_ground/2
            if a<=-self.taper.w/2+self.taper.t_left or b>=self.taper.w/2-self.taper.t_right:
                raise ValueError('Inner ground shield touches a cap leg')
            if not (b<=self.x_meander-inductor_width/2-4 or a>=self.x_meander+inductor_width/2+4):
                raise ValueError('Inner ground shield must stay beside the inductor')
        if self.w_inner_ground_hole and (not self.t_inner_ground or self.w_inner_ground_hole%4 or not 4<=self.w_inner_ground_hole<=self.t_inner_ground-8):
            raise ValueError('Ground shield opening must leave at least 4 um on both sides')
        if self.x_meander-inductor_width/2<=-self.taper.w/2+self.taper.t_left or self.x_meander+inductor_width/2>=self.taper.w/2-self.taper.t_right:
            raise ValueError('Offset meander touches a cap leg')
        if self.t_coupling_finger:
            t1=self.taper1
            left_t=t1.t_left if t1 else self.t_t1_left
            left_l=t1.l_left if t1 else self.l_t1_left
            coil_left=(t1.w-t1.meander.w)/2 if t1 else 44
            if self.l_coupling_bridge%2 or not left_l+18<=self.l_coupling_bridge<=606 or coil_left-left_t-self.t_coupling_finger<8:
                raise ValueError('Coupling finger needs at least 4 um clearance to t1 leg and coil')
            if self.t_coupling_finger%2 or self.t_coupling_finger<2:
                raise ValueError('Coupling finger width must be positive and even')
            if self.l_coupling_finger%2 or not 4<=self.l_coupling_finger<=self.l_coupling_bridge-40:
                raise ValueError('Coupling finger length must be even and clear the cap top')
            if 8+self.taper.l_right<self.l_coupling_bridge-4:raise ValueError('Right leg must reach the coupling-finger bridge')
        t1w=self.taper1.w if self.taper1 else 328
        if self.w_t1_ground_slot and (self.w_t1_ground_slot%4 or not 4<=self.w_t1_ground_slot<=t1w-4 or self.l_t1_ground_slot%2 or not 2<=self.l_t1_ground_slot<=444):
            raise ValueError('Invalid T1 ground slot')
        if self.w_upper_plate and (self.w_upper_plate%2 or not 2<=self.w_upper_plate<=t1w-2):
            raise ValueError('Upper plate must fit over t1 on the 2 um grid')
        if self.l_cal_output%2 or not 16<=self.l_cal_output<610:raise ValueError('Invalid calibration output finger')
        if self.n_coupling_teeth:
            t1=self.taper1 or TaperParams(w=328,t_left=self.t_t1_left,l_left=self.l_t1_left,
                                       meander=MeanderParams(w=240,l=596,l_top=10,n=58,t=4,d=6))
            m=t1.meander
            if not self.t_coupling_finger:raise ValueError('Coupling teeth require the driven finger')
            if not isinstance(self.n_coupling_teeth,int) or not 1<=self.n_coupling_teeth<=m.n//2-1:
                raise ValueError('Coupling teeth must use open even-row gaps')
            if m.d<6 or (m.d-2)%4:raise ValueError('Need room for a 2 um tooth with even clearances')
            post_x=2*np.floor((t1.t_left+(t1.w-m.w)/2-self.t_coupling_finger)/4+.5)
            if self.w_coupling_teeth<=0 or self.w_coupling_teeth%2 or post_x+self.t_coupling_finger+self.w_coupling_teeth>(t1.w+m.w)/2-m.t-2:
                raise ValueError('Coupling tooth must clear the right U-turn')
            first_top=-8-t1.t_top-m.l_top-2*(m.t+m.d)-m.t-(m.d-2)/2
            last_bottom=first_top-2-2*(self.n_coupling_teeth-1)*(m.t+m.d)
            if first_top>-self.l_coupling_bridge+self.l_coupling_finger or last_bottom<-self.l_coupling_bridge:
                raise ValueError('Coupling teeth must join the driven finger')
        if 8+self.taper.t_top+self.d_idc_top+(2*self.idc_pairs-1)*(self.idc_t+self.idc_d)+self.idc_t>8+self.taper.l_left:
            raise ValueError('IDC fingers must remain within the taper left leg')


def build(p:SecondParams, *, isolation=False):
    """Isolation retains the taper comb but uses weak feed coupling for L calibration.

    No inductor is placed in either calibration feed. Mixed runs include the
    complete t1 and regular resonators and the usual feed-connected output end.
    """
    side=412;sep=14;height=1500;top_height=448;bottom=-618
    input_left=side+sep
    t2_left=input_left+p.idc_span
    t2_x=t2_left+p.taper.w/2
    t1=p.taper1 or TaperParams(w=328,t_top=14,t_left=p.t_t1_left,t_right=6,
                             l_left=p.l_t1_left,l_right=310,
                             meander=MeanderParams(w=240,l=596,l_top=10,n=58,t=4,d=6))
    if 8+t1.t_top+t1.meander.l!=618:raise ValueError('T1 must end at the shared bottom ground')
    t1_left=t2_left+p.taper.w+p.d_t2_t1
    t1_x=t1_left+t1.w/2
    reg_left=t1_left+t1.w+p.d_t1_regular
    reg_x=reg_left+151
    end_left=reg_left+302+56;end_x=end_left+151
    if isolation:
        output_right=t1_left+t1.w
    else:output_right=end_left+302
    width=output_right+sep+side
    top=gdstk.Cell('TOP')
    top.add(gdstk.Reference(make_taper(p.taper,name='TAPER_2',meander_offset_x=p.x_meander,t_ground_tail=p.t_ground_tail,straight_inductor=p.straight_inductor),origin=(t2_x,0)))
    # The taper comb is part of the t2 conductor, even in its L calibration.
    parts=[]
    if p.w_ground_pad:
        # Identical pads alongside the straight inductor, as in the regular cell.
        # Each touches the original bottom ground and clears its adjacent cap leg.
        left=t2_left+p.taper.t_left+p.d_ground_pad
        right=t2_left+p.taper.w-p.taper.t_right-p.d_ground_pad
        parts.extend([gdstk.rectangle((left,bottom),(left+p.w_ground_pad,bottom+p.l_ground_pad)),
                      gdstk.rectangle((right-p.w_ground_pad,bottom),(right,bottom+p.l_ground_pad))])
    if p.t_inner_ground:
        edge=-8-p.taper.t_top-p.d_inner_ground_top
        shield=gdstk.rectangle((t2_x+p.x_inner_ground-p.t_inner_ground/2,bottom),
                               (t2_x+p.x_inner_ground+p.t_inner_ground/2,edge))
        if p.w_inner_ground_hole:
            hole=gdstk.rectangle((t2_x+p.x_inner_ground-p.w_inner_ground_hole/2,bottom+4),
                                 (t2_x+p.x_inner_ground+p.w_inner_ground_hole/2,edge-4))
            parts.extend(gdstk.boolean([shield],[hole],'not',precision=.001))
        else:parts.append(shield)
    if p.w_top_left_extension:
        parts.append(gdstk.rectangle((t2_left-p.w_top_left_extension,-8-p.taper.t_top),
                                     (t2_left,-8)))
    pitch=p.idc_t+p.idc_d
    first_y=-8-p.taper.t_top-p.d_idc_top
    for row in range(2*p.idc_pairs):
        y=first_y-row*pitch
        if row%2:  # Node fingers extend from taper left leg toward the input trunk.
            parts.append(gdstk.rectangle((input_left+20+p.idc_d,y-p.idc_t),(t2_left+p.taper.t_left,y)))
        elif not isolation:
            parts.append(gdstk.rectangle((input_left,y-p.idc_t),(t2_left-p.idc_d,y)))
    # Input trunk and 10 um feed are isolated from the t2 node and ground.
    trunk_bottom=first_y-(2*p.idc_pairs-1)*pitch-p.idc_t
    if p.equal_idc_sides:
        trunk_bottom=-8-p.taper.l_left
    parts.extend([gdstk.rectangle((input_left,trunk_bottom),(input_left+20,-8)),
                  gdstk.rectangle((0,-18),(input_left+20,-8))])
    if p.w_upper_plate:
        edge=t2_left+p.taper.w
        parts.extend([gdstk.rectangle((edge-4,-8),(edge,-2)),
                      gdstk.rectangle((edge-4,-6),(t1_left+p.w_upper_plate,-2))])
    if p.t_coupling_finger:
        # Go below t1's left leg (ends y=-592), above ground (starts -618).
        # The 20 um inner gap is between the 24 um cap leg and the 240 um coil.
        coil_left=(t1.w-t1.meander.w)/2
        x=t1_left+2*np.floor((t1.t_left+coil_left-p.t_coupling_finger)/4+.5)
        bridge=-p.l_coupling_bridge
        parts.extend([gdstk.rectangle((t2_left+p.taper.w-p.taper.t_right,bridge),(x+p.t_coupling_finger,bridge+4)),
                      gdstk.rectangle((x,bridge),(x+p.t_coupling_finger,bridge+p.l_coupling_finger))])
        for i in range(p.n_coupling_teeth):
            m=t1.meander;row=2*(i+1)
            row_top=-8-t1.t_top-m.l_top-row*(m.t+m.d)
            tooth_top=row_top-m.t-(m.d-2)/2
            parts.append(gdstk.rectangle((x,tooth_top-2),(x+p.t_coupling_finger+p.w_coupling_teeth,tooth_top)))
    if isolation:
        # A finger rather than a second resonator terminates the opposite port.
        parts.extend([gdstk.rectangle((t1_left,-22),(output_right,-8)),
                      gdstk.rectangle((t1_left,-8-p.l_cal_output),(t1_left+20,-8)),
                      gdstk.rectangle((t1_x,-18),(width,-8))])
        centers=[input_left+10,t2_x,t1_x]
    else:
        top.add(gdstk.Reference(make_taper(t1,name='TAPER_1'),origin=(t1_x,0)))
        regular=make_unit_cell(UnitCellParams(meander=MeanderParams(w=186)),name='REGULAR')
        top.add(gdstk.Reference(regular,origin=(reg_x,0)),gdstk.Reference(regular,origin=(end_x,0)))
        parts.append(gdstk.rectangle((end_x,-18),(width,-8)))
        centers=[input_left+10,t2_x,t1_x,reg_x,end_x]
    # Continuous upper ground; strips only at the weak regular-cell gaps.
    parts.extend([gdstk.rectangle((0,0),(width,top_height)),
                  gdstk.rectangle((0,top_height-height),(width,bottom)),
                  gdstk.rectangle((0,top_height-height),(side,-24)),
                  gdstk.rectangle((width-side,top_height-height),(width,-24))])
    if p.l_ground_rise:
        x=t2_x+p.x_meander
        parts.append(gdstk.rectangle((x-p.w_ground_riser/2,bottom),(x+p.w_ground_riser/2,bottom+p.l_ground_rise)))
    if not isolation:
        x=t1_left+t1.w+p.d_t1_regular/2+p.x_junction_ground
        if p.t_junction_ground==4:
            parts.append(gdstk.rectangle((x-2,bottom),(x+2,0)))
        else:
            # Keep both edges on the 2 um grid. For an even center, remove
            # the left half of the original strip; its regular-side edge stays.
            edge=2*np.floor(x/2+.5)
            parts.append(gdstk.rectangle((edge,bottom),(edge+2,0)))
        x=reg_left+302+28
        parts.append(gdstk.rectangle((x-2,bottom),(x+2,0)))
    metal=gdstk.boolean(parts,[],'or',precision=.001,layer=1)
    if p.w_t1_ground_slot:
        hole=gdstk.rectangle((t1_x-p.w_t1_ground_slot/2,2),(t1_x+p.w_t1_ground_slot/2,2+p.l_t1_ground_slot))
        metal=gdstk.boolean(metal,[hole],'not',precision=.001,layer=1)
    top.add(*metal)
    # `unit` below is the legacy port-template scaffold (regular end cap).
    # Actual T2/T1 capacitor dimensions are stored in `second` above it.
    meta={'layout_kind':'second_taper','second':asdict(p),'isolation':isolation,
          'w_sim':width,'cell_centers':centers,'expected_metal_components':4 if isolation else 2,
          'parameters':{'n_cells':1 if isolation else 3,'l_sim':height,'l_top_ground':top_height,
                        'w_side':side,'t_feed':10,'unit':{'d_top':8,'cap':{'t_top':14}}}}
    # Keep archived zero-shield metadata and geometry identities unchanged.
    if not p.t_inner_ground:
        for key in ('t_inner_ground','x_inner_ground','d_inner_ground_top'):meta['second'].pop(key)
    if not p.w_inner_ground_hole:meta['second'].pop('w_inner_ground_hole')
    if not p.w_ground_pad:
        for key in ('w_ground_pad','l_ground_pad','d_ground_pad'):meta['second'].pop(key)
    if not p.equal_idc_sides:meta['second'].pop('equal_idc_sides')
    if not p.x_junction_ground:meta['second'].pop('x_junction_ground')
    if p.x_meander:meta['inductance_pin_x']=t2_x+p.x_meander
    return top,meta


def from_dict(value):
    value=dict(value);taper=dict(value['taper'])
    wire_type=StraightWireParams if value.get('straight_inductor') else MeanderParams
    taper['meander']=wire_type(**taper['meander'])
    value['taper']=TaperParams(**taper)
    if value.get('taper1') is not None:
        t1=dict(value['taper1']);t1['meander']=MeanderParams(**t1['meander']);value['taper1']=TaperParams(**t1)
    return SecondParams(**value)


def write(p,basename,**kwargs):
    top,meta=build(p,**kwargs);path=Path(basename).with_suffix('.gds')
    path.parent.mkdir(parents=True,exist_ok=True)
    lib=gdstk.Library(unit=1e-6,precision=1e-9);lib.add(top,*top.dependencies(True));lib.write_gds(str(path))
    path.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
    return path


def require_same_isolation_geometry(p,reference_gds):
    """Reuse L only if the complete calibration metal is exactly unchanged.

    An embedded T1 right-leg length is not drawn in the isolation fixture.
    Its metadata can therefore differ without requiring another identical EM
    calibration. Compare the actual metal rather than assuming equivalence.
    """
    reference_gds=Path(reference_gds)
    metadata=json.loads(reference_gds.with_suffix('.json').read_text())
    if not metadata['isolation']:raise ValueError('Expected an isolation reference')
    _,polys=validate(reference_gds,metadata)
    expected,_=build(p,isolation=True)
    if gdstk.boolean(expected.get_polygons(),polys,'xor',precision=.001):
        raise ValueError('L reference has different isolation metal')


def validate(gds,meta):
    expected,rebuilt=build(from_dict(meta['second']),isolation=meta['isolation'])
    for optional in ('d_idc_top','x_meander','t_ground_tail','t_coupling_finger','l_coupling_finger','taper1',
                     'n_coupling_teeth','w_coupling_teeth','w_upper_plate','l_coupling_bridge','l_cal_output',
                     'w_t1_ground_slot','l_t1_ground_slot','t_junction_ground','l_ground_rise','w_ground_riser','straight_inductor',
                     'w_top_left_extension'):
        if optional not in meta['second']:rebuilt['second'].pop(optional)
    if rebuilt!=meta:raise ValueError('Stale second-taper metadata')
    lib=gdstk.read_gds(str(gds));tops=lib.top_level()
    if len(tops)!=1 or tops[0].name!='TOP' or not np.isclose(lib.unit,1e-6,rtol=1e-9,atol=0):
        raise ValueError('Invalid GDS library')
    polys=tops[0].get_polygons()
    if any(poly.layer!=1 or poly.datatype!=0 for poly in polys):raise ValueError('Wrong metal layer')
    if gdstk.boolean(polys,expected.get_polygons(),'xor',precision=.001):raise ValueError('GDS differs from parameters')
    count=len(gdstk.boolean(polys,[],'or',precision=.001))
    if count!=meta['expected_metal_components']:raise ValueError(f'Unexpected connectivity: {count}')
    pars=meta['parameters']
    return SimpleNamespace(n_cells=pars['n_cells'],l_sim=pars['l_sim'],l_top_ground=pars['l_top_ground']),polys
