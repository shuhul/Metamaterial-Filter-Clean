"""The fixture correction is confined to two small regular-cell cap tabs."""
import unittest
import numpy as np
import gdstk
from geometry.taper.run_locked_taper2 import GEOMETRY
from geometry.taper.second_layout import build
from geometry.taper.regular_cap_tab_case import build_tabbed,RegularCapTab
from geometry.taper.locked_taper2 import T1

class RegularTabTests(unittest.TestCase):
    def test_only_regular_cap_tabs_are_added_on_grid(self):
        base,old=build(GEOMETRY);parts,meta,regions=build_tabbed(GEOMETRY,RegularCapTab(w=100,t=4))
        diff=gdstk.boolean(parts,base.get_polygons(),'xor',precision=.001)
        self.assertAlmostEqual(sum(p.area() for p in diff),2*(100*4+4**2),places=6)
        self.assertTrue(all(np.allclose(p.points/2,np.round(p.points/2)) for p in diff))
        self.assertEqual(len(gdstk.boolean(parts,[],'or',precision=.001)),2)
        left=meta['cell_centers'][2]-T1.w/2;right=meta['cell_centers'][2]+T1.w/2
        t1=gdstk.rectangle((left,-1052),(right,448))
        self.assertFalse(gdstk.boolean(gdstk.boolean(parts,[t1],'and'),gdstk.boolean(base.get_polygons(),[t1],'and'),'xor'))
        self.assertTrue(all(q.bounding_box()[0][0]>right for q in regions))

    def test_zero_tab_preserves_the_original_design(self):
        base,_=build(GEOMETRY);parts,_,regions=build_tabbed(GEOMETRY,RegularCapTab(w=0))
        self.assertFalse(gdstk.boolean(parts,base.get_polygons(),'xor',precision=.001))
        self.assertEqual(regions,[])
        with self.assertRaises(ValueError):build_tabbed(GEOMETRY,RegularCapTab(w=100,t=8))

if __name__=='__main__':unittest.main()
