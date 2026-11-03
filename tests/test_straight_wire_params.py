from dataclasses import replace,asdict
import unittest
import gdstk
from geometry.meander import MeanderParams
from geometry.straight_wire import StraightWireParams
from geometry.taper.plain_gap_study import baseline
from geometry.taper.clearance_extraction import build


class StraightWireParamsTest(unittest.TestCase):
    def test_adapter_does_not_change_saved_geometry(self):
        p=baseline();old=replace(p,taper=replace(p.taper,meander=MeanderParams(**asdict(p.taper.meander))))
        a,ma=build(p,'idc',200);b,mb=build(old,'idc',200)
        self.assertEqual(ma,mb);self.assertFalse(gdstk.boolean(a,b,'xor',precision=.001))

    def test_wide_straight_wire_has_no_unused_turn_constraint(self):
        p=baseline();m=StraightWireParams(w=512,l=600,t=512,n=3,l_top=10,d=6)
        q=replace(p,taper=replace(p.taper,w=1800,meander=m))
        metal,_=build(q,'idc',200)
        self.assertTrue(metal)
        with self.assertRaises(ValueError):MeanderParams(w=1600,l=600,t=512,n=3)
        with self.assertRaises(ValueError):replace(q,straight_inductor=False)


if __name__=='__main__':unittest.main()
