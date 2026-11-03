"""User-required physical constraints for the stricter T1 campaign."""
import unittest
import gdstk
from geometry.taper.three_tapers_v3 import build,BASE


class GroundPlacementTests(unittest.TestCase):
    def test_wider_strip_stays_centered_on_grid_without_input_metal(self):
        for isolation,width,gap in ((False,8,84),(True,8,84),(False,6,74),(True,6,74)):
            top,meta=build(BASE,isolation=isolation,d_input=50,d_junction=gap,t_junction_ground=width)
            edge=meta['cell_centers'][1]+BASE.cap.w/2
            mid=edge+((96+width%4) if isolation else gap)/2
            half=width/2
            self.assertEqual(gdstk.inside([(mid-half+1,-300),(mid+half-1,-300),
                                           (mid-half-1,-300),(mid+half+1,-300)],top.get_polygons()),
                             (True,True,False,False))
            self.assertEqual((mid-half)%2,0)
            self.assertEqual((mid+half)%2,0)
            opening=gdstk.rectangle((426+BASE.cap.w,-618),(426+BASE.cap.w+50,0))
            self.assertFalse(gdstk.boolean(top.get_polygons(),[opening],'and',precision=.001))

    def test_input_gap_is_empty_and_regular_strip_is_centered(self):
        top,meta=build(BASE,d_input=96,d_junction=64)
        metal=top.get_polygons()
        input_edge=426+BASE.cap.w
        opening=gdstk.rectangle((input_edge,-618),(input_edge+96,0))
        self.assertFalse(gdstk.boolean(metal,[opening],'and',precision=.001))
        t1_edge=meta['cell_centers'][1]+BASE.cap.w/2
        regular_edge=meta['cell_centers'][2]-151
        mid=(t1_edge+regular_edge)/2
        self.assertEqual(mid-t1_edge,regular_edge-mid)
        self.assertEqual(gdstk.inside([(mid,-300),(mid-3,-300),(mid+3,-300)],metal),
                         (True,False,False))
        strip=gdstk.rectangle((mid-2,-618),(mid+2,0))
        self.assertFalse(gdstk.boolean([strip],metal,'not',precision=.001))
        upper=gdstk.rectangle((0,0),(meta['w_sim'],448))
        self.assertFalse(gdstk.boolean([upper],metal,'not',precision=.001))

    def test_calibration_input_sweep_keeps_core_and_output_geometry_fixed(self):
        reference,meta=build(BASE,isolation=True,d_input=96)
        left=meta['cell_centers'][1]-BASE.cap.w/2
        region=gdstk.rectangle((left,-1052),(meta['w_sim'],448))
        expected=gdstk.boolean(reference.get_polygons(),[region],'and',precision=.001)
        for gap in (40,42):
            changed,other=build(BASE,isolation=True,d_input=gap)
            delta=meta['cell_centers'][1]-other['cell_centers'][1]
            aligned=[q.copy().translate(delta,0) for q in changed.get_polygons()]
            actual=gdstk.boolean(aligned,[region],'and',precision=.001)
            self.assertFalse(gdstk.boolean(expected,actual,'xor',precision=.001))


if __name__=='__main__':unittest.main()
