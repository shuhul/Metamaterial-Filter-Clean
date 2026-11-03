"""Check the exact metal changes in feed-driven bypass fixtures, without EM."""
import unittest
import gdstk
import numpy as np
from geometry.taper import bridged_fixtures as fixtures
from geometry.taper.locked_taper2 import VERIFIED


class BridgedFixturesTests(unittest.TestCase):
    def test_lateral_shield_is_grounded_and_does_not_widen_wire(self):
        from dataclasses import replace
        p=fixtures.t2_parameters()
        p=replace(p,x_meander=-120,t_inner_ground=4)
        full,meta,_=fixtures.build(fixtures.CASES[2],second=p)
        x=meta['cell_centers'][1]
        self.assertEqual(gdstk.inside([(x,-26),(x,-400),(x,-617),(x,-24),
                                     (x-120,-600),(x-100,-600)],full),
                         (True,True,True,False,True,False))
        unshielded,_,_=fixtures.build(fixtures.CASES[2],second=replace(p,t_inner_ground=0))
        difference=gdstk.boolean(full,unshielded,'xor',precision=.001)
        shield=gdstk.rectangle((x-2,-618),(x+2,-26))
        self.assertFalse(gdstk.boolean(difference,[shield],'xor',precision=.001))

    def test_wire_control_removes_only_wire_and_keeps_direct_feed(self):
        from dataclasses import replace
        p=fixtures.t2_parameters()
        full,_,_=fixtures.build(fixtures.CASES[2],second=p)
        control,meta,_=fixtures.build(fixtures.CASES[2],second=p,omit_t2_wire=True)
        x=meta['cell_centers'][1]
        removed=gdstk.rectangle((x-p.taper.meander.t/2,-618),(x+p.taper.meander.t/2,-22))
        self.assertFalse(gdstk.boolean(gdstk.boolean(full,control,'xor',precision=.001),[removed],'xor',precision=.001))
        self.assertEqual(gdstk.inside([(436,-100),(570,-34),(570,-13)],control),(False,False,True))
        shifted=replace(p,x_meander=-100)
        other,_,_=fixtures.build(fixtures.CASES[2],second=shifted,omit_t2_wire=True)
        self.assertFalse(gdstk.boolean(control,other,'xor',precision=.001))

    def test_t1_removes_upstream_electrode_and_extends_original_feed(self):
        bridged,meta,_=fixtures.build(fixtures.CASES[0])
        original,_,_=fixtures.build(fixtures.CASES[1])
        difference=gdstk.boolean(bridged,original,'xor',precision=.001)
        allowed=gdstk.rectangle((426,-544),(760,-8))
        self.assertFalse(gdstk.boolean(difference,[allowed],'not',precision=.001))
        # Old finger is gone; old 14 um cap top becomes a 10 um feed.
        self.assertEqual(gdstk.inside([(744,-300),(600,-20),(600,-13),(757,-13)],bridged),
                         (False,False,True,True))
        self.assertEqual(gdstk.inside([(744,-300),(600,-20),(600,-13),(757,-13)],original),
                         (True,True,True,False))
        reference=gdstk.read_gds(VERIFIED['artifacts']['fixture_gds']).top_level()[0].get_polygons()
        self.assertFalse(gdstk.boolean(original,reference,'xor',precision=.001))
        self.assertEqual(meta['parameters']['n_cells'],2)

    def test_t2_removes_both_combs_and_connects_top_feed(self):
        direct,meta,_=fixtures.build(fixtures.CASES[2])
        original,_,_=fixtures.build(fixtures.CASES[3])
        self.assertEqual(gdstk.inside([(436,-100),(570,-34),(570,-13)],original),(True,True,False))
        self.assertEqual(gdstk.inside([(436,-100),(570,-34),(570,-13)],direct),(False,False,True))
        difference=gdstk.boolean(direct,original,'xor',precision=.001)
        allowed=gdstk.rectangle((426,-592),(716,-8))
        self.assertFalse(gdstk.boolean(difference,[allowed],'not',precision=.001))
        self.assertEqual(meta['parameters']['n_cells'],3)

    def test_t2_has_uniform_centred_wire_and_no_pedestal(self):
        p=fixtures.t2_parameters()
        self.assertEqual(p.l_ground_rise,0)
        self.assertIsNone(p.t_ground_tail)
        self.assertEqual(p.x_meander,0)
        for name in fixtures.CASES[2:]:
            polygons,meta,_=fixtures.build(name)
            x=meta['cell_centers'][1]
            for y in (-100,-400,-600):
                self.assertEqual(gdstk.inside([(x-15,y),(x+15,y),(x-17,y),(x+17,y)],polygons),
                                 (True,True,False,False))
            self.assertEqual(meta['second']['taper1'],fixtures.VERIFIED['geometry'])


if __name__=='__main__':unittest.main()
