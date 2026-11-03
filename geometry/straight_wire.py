"""Straight-wire parameter adapter for existing serialized taper geometry.

The renderer uses only l and t. Legacy meander fields stay in JSON so archived
cases rebuild unchanged, but unused turn-count/width constraints do not limit
the width of a straight conductor.
"""
from dataclasses import dataclass
from math import isfinite
from geometry.meander import MeanderParams


@dataclass(frozen=True)
class StraightWireParams(MeanderParams):
    def __post_init__(self):
        for name in ('w','l_top','l','t','d'):
            if not isfinite(getattr(self,name)) or getattr(self,name)<=0:
                raise ValueError(f'straight wire {name} must be finite and positive')
        if not 0<self.l_top<self.l:raise ValueError('Legacy lead metadata must fit the wire length')
        if isinstance(self.n,bool) or not isinstance(self.n,int) or self.n<3:
            raise ValueError('Legacy run-count metadata must be an integer >=3')

    @property
    def l_body(self):
        return 0  # No horizontal runs are rendered.

    @property
    def l_bot(self):
        return self.l-self.l_top
