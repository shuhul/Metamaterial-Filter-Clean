"""Verify that T1 is removed without changing the coupling or regular cells."""
import unittest
import gdstk
from geometry.taper.frank_coupling_controls import build,original_build


class FrankCouplingControlsTests(unittest.TestCase):
    def test_all_retained_coupling_and_downstream_metal_is_identical(self):
        original,_,_=original_build('01_t1_input_bridged')
        for case in ('minimal_feed10','minimal_feed14','frame_without_inductor'):
            polygons,meta,_=build(case)
            leg=meta['retained_coupling_leg']
            clip=gdstk.rectangle((leg['x_left'],-1052),(meta['w_sim'],448))
            expected=gdstk.boolean(original,[clip],'and',precision=.001)
            actual=gdstk.boolean(polygons,[clip],'and',precision=.001)
            self.assertFalse(gdstk.boolean(expected,actual,'xor',precision=.001))
            x=meta['cell_centers'][1]
            self.assertFalse(gdstk.inside([(x,-300)],polygons)[0])
            self.assertTrue(all(gdstk.inside([(x,-13),(x,-619),(leg['x_left']+2,-300)],polygons)))
            left=meta['cell_centers'][1]-meta['taper']['w']/2
            present=gdstk.inside([(left+2,-300)],polygons)[0]
            self.assertEqual(present,case=='frame_without_inductor')

    def test_only_feed_thickness_changes_between_minimal_controls(self):
        a,ma,_=build('minimal_feed10');b,mb,_=build('minimal_feed14')
        delta=gdstk.boolean(a,b,'xor',precision=.001)
        left=ma['cell_centers'][1]-ma['taper']['w']/2
        rectangle=gdstk.rectangle((left,-22),(ma['retained_coupling_leg']['x_left'],-18))
        self.assertFalse(gdstk.boolean(delta,[rectangle],'xor',precision=.001))


if __name__=='__main__':unittest.main()
