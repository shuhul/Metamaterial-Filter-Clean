"""Check that the user's removed input resonator really is absent."""
import unittest
import gdstk
from geometry.taper.three_tapers_v2 import BASE,build


class CorrectedTaperTests(unittest.TestCase):
    def test_left_region_contains_only_facing_electrode_and_feed(self):
        top,meta=build(BASE,d_input=40,d_junction=52)
        polys=top.get_polygons()
        # Interior of the removed left meander/pads must be completely empty.
        empty=gdstk.rectangle((450,-610),(744,-24))
        self.assertFalse(gdstk.boolean(polys,[empty],'and',precision=.001))
        self.assertFalse(meta['input_end_resonator'])
        self.assertEqual(len(gdstk.boolean(polys,[],'or',precision=.001)),2)
        flat=gdstk.rectangle((0,0),(meta['w_sim'],448))
        self.assertFalse(gdstk.boolean([flat],polys,'not',precision=.001))

    def test_isolation_has_exactly_one_meander(self):
        top,meta=build(BASE,isolation=True)
        meanders=[c for c in top.dependencies(True) if c.name.endswith('_MEANDER')]
        self.assertEqual(len(meanders),1)
        self.assertEqual(len(gdstk.boolean(top.get_polygons(),[],'or',precision=.001)),3)

    def test_junction_ground_stays_at_regular_cell_clearance(self):
        top,meta=build(BASE,d_input=28,d_junction=80,d_regular_ground=26)
        regular_left=meta['cell_centers'][2]-151
        self.assertEqual(gdstk.inside([(regular_left-28,-300),(regular_left-25,-300)],top.get_polygons()),
                         (True,False))

    def test_even_input_gap_preserves_t1_side_ground_clearance(self):
        top,meta=build(BASE,d_input=26,d_junction=80,d_input_ground=12,d_regular_ground=26)
        t1_left=meta['cell_centers'][1]-BASE.cap.w/2
        self.assertEqual(gdstk.inside([(t1_left-14,-300),(t1_left-11,-300)],top.get_polygons()),
                         (True,False))
        self.assertEqual(len(gdstk.boolean(top.get_polygons(),[],'or',precision=.001)),2)


if __name__=='__main__':unittest.main()
