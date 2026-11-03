"""Compact inductor connectivity, layers, and independent probe partition."""
from dataclasses import replace
import unittest
import gdstk
from geometry.meander import MeanderParams
from geometry.taper.run_locked_taper2 import GEOMETRY
from geometry.taper.compact_inductor_case import build_compact

class CompactTaperTests(unittest.TestCase):
    def parameters(self):
        return replace(GEOMETRY,w_top_left_extension=0,l_ground_rise=0,w_ground_riser=24,
                       x_meander=-120,taper=replace(GEOMETRY.taper,w=400,l_right=536))

    def test_replaces_original_wire_with_single_connected_compact_branch(self):
        polys,meta=build_compact(self.parameters(),MeanderParams(w=32,l=90,l_top=8,n=3,t=4,d=4))
        self.assertTrue(all(p.layer==1 and p.datatype==0 for p in polys))
        self.assertEqual(len(gdstk.boolean(polys,[],'or',precision=.001)),2)
        self.assertFalse(gdstk.inside([(678,-300)],polys)[0])
        self.assertTrue(all(gdstk.inside([(520,-600),(580,-526)],polys)))
        self.assertEqual(meta['series_probe_T2_x'],520)
        # The wire must be its own polygon; a full ground union breaks the
        # nonperturbing series-probe construction on multi-resonator devices.
        q=[p for p in polys if all(gdstk.inside([(520,-615.5),(520,-616.5)],[p]))]
        self.assertEqual(len(q),1)
        self.assertEqual(len(gdstk.boolean(q,[gdstk.rectangle((0,-616),(2500,500))],'and')),1)
        self.assertEqual(len(gdstk.boolean(q,[gdstk.rectangle((0,-1200),(2500,-616))],'and')),1)

    def test_compact_component_cannot_touch_the_IDC(self):
        with self.assertRaisesRegex(ValueError,'touches the IDC'):
            build_compact(self.parameters(),MeanderParams(w=32,l=108,l_top=8,n=3,t=4,d=4))

    def test_isolation_keeps_separate_calibration_feeds(self):
        polys,meta=build_compact(self.parameters(),MeanderParams(w=32,l=90,l_top=8,n=3,t=4,d=4),isolation=True)
        self.assertEqual(len(gdstk.boolean(polys,[],'or',precision=.001)),4)
        self.assertTrue(meta['isolation'])
        self.assertEqual(meta['parameters']['n_cells'],1)

if __name__=='__main__':unittest.main()
