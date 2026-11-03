import unittest
import gdstk
from geometry.taper.return_rail_trial import build,configuration
from geometry.taper.compact_stages import build as original


class ReturnRailGeometryTests(unittest.TestCase):
    def test_rails_clear_wire_bar_and_do_not_create_a_pedestal(self):
        p,rails=configuration();polys,meta=build(p,rails,'direct');x=meta['cell_centers'][1]
        self.assertEqual(gdstk.inside([(x,-300),(x+124,-300),(x-124,-300)],polys),(True,True,True))
        self.assertEqual(gdstk.inside([(x+50,-300),(x-50,-300),(x+124,-24),(x+10,-600)],polys),(False,False,False,False))
        self.assertEqual(gdstk.inside([(x+124,-617),(x+124,-619)],polys),(True,True))

    def test_T1_and_strong_junction_are_unchanged(self):
        p,rails=configuration();polys,meta=build(p,rails,'idc');base,_,_=original(p,'idc')
        edge=meta['cell_centers'][1]+p.taper.w/2
        roi=gdstk.rectangle((edge-20,-1052),(meta['w_sim'],448))
        a=gdstk.boolean(polys,[roi],'and',precision=.001)
        b=gdstk.boolean(base,[roi],'and',precision=.001)
        self.assertFalse(gdstk.boolean(a,b,'xor',precision=.001))


if __name__=='__main__':unittest.main()
