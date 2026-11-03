"""Optional symmetric T2-to-ground IDC, confined to the T2 footprint.

This is distinct from the input IDC and from the rejected strip above T1.
No geometry is generated merely by importing this module. A physical trial
still requires independent inductance and ordinary two-port validation.
"""
from dataclasses import dataclass,asdict
import gdstk


@dataclass(frozen=True)
class GroundIDCParams:
    n: int = 14      # Finger pairs per half; two symmetric halves.
    t: int = 2       # Horizontal finger thickness, um.
    d: int = 2       # Horizontal gap, um.
    l: int = 120     # Vertical overlap, um.
    d_edge: int = 8  # Inset from each outer T2 bar edge, um.

    def __post_init__(self):
        if not isinstance(self.n,int) or self.n<2:raise ValueError('Need at least two pairs per half')
        if any(v<2 or v%2 for v in (self.t,self.d,self.l,self.d_edge)):
            raise ValueError('Ground IDC dimensions must be positive on the 2 um grid')
        if self.l<12:raise ValueError('Overlap is too short')

    @property
    def w(self):return 2*self.n*self.t+(2*self.n-1)*self.d


def apply(polygons,meta,p):
    """Cut two local ground slots and add alternating vertical teeth.

    Coordinate origin follows the GDS: upper-ground lower edge y=0, T2 bar
    y=-18..-8. The central ground corridor remains intact above the wire.
    """
    second=meta['second'];frame=second['taper'];center=meta['cell_centers'][1]
    if frame['t_top']!=10 or second.get('w_upper_plate',0):raise ValueError('Use the plain 10 um T2 bar')
    clear=frame['w']-2*(p.d_edge+p.w)
    if clear<frame['meander']['t']+8:raise ValueError('Keep ground intact above the centered wire')
    if p.d_edge<p.d+2 or p.l-4>meta['parameters']['l_top_ground']-8:
        raise ValueError('Ground slots must remain inside T2 and below the top box wall')
    left=center-frame['w']/2
    starts=[left+p.d_edge,left+frame['w']-p.d_edge-p.w]
    holes=[gdstk.rectangle((x-p.d,0),(x+p.w+p.d,p.l-4)) for x in starts]
    metal=gdstk.boolean(polygons,holes,'not',precision=.001,layer=1)
    for start in starts:
        for row in range(2*p.n):
            x=start+row*(p.t+p.d)
            bottom,top=(-10,p.l-6) if row%2==0 else (-6,p.l-2)
            metal.append(gdstk.rectangle((x,bottom),(x+p.t,top),layer=1))
    metal=gdstk.boolean(metal,[],'or',precision=.001,layer=1)
    result=dict(meta);result['ground_idc']=asdict(p)
    result['title']=meta['title']+' with T2 Ground IDC'
    result['expected_metal_components']=len(metal)
    return metal,result
