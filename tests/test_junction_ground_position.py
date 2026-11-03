"""Verify a junction-position parameter cannot alter the approved cell shapes."""
from dataclasses import replace
import json
import unittest
import gdstk
from geometry.taper.plain_gap_study import baseline,BASE_GDS
from geometry.taper.clearance_extraction import build,validate


class JunctionGroundPositionTest(unittest.TestCase):
    def test_archived_default_geometry_still_validates(self):
        metadata=json.loads(BASE_GDS.with_suffix('.json').read_text())
        validate(BASE_GDS,metadata)

    def test_only_existing_ground_strip_moves(self):
        p=baseline();old,meta=build(p,'idc',200)
        new,_=build(replace(p,x_junction_ground=8),'idc',200)
        right=meta['cell_centers'][2]+p.taper1.w/2
        roi=gdstk.rectangle((right,-618),(right+p.d_t1_regular,0))
        difference=gdstk.boolean(old,new,'xor',precision=.001)
        self.assertFalse(gdstk.boolean(difference,[roi],'not',precision=.001))
        self.assertAlmostEqual(sum(q.area() for q in difference),2*4*618,places=5)
        self.assertNotIn('x_junction_ground',meta['second'])
        self.assertFalse(gdstk.boolean([gdstk.rectangle((0,0),(meta['w_sim'],448))],new,'not',precision=.001))

    def test_no_off_grid_or_touching_positions(self):
        for offset in (1,26,-26):
            with self.assertRaises(ValueError):replace(baseline(),x_junction_ground=offset)


if __name__=='__main__':unittest.main()
