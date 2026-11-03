"""Physical invariants for inspection fixtures; these tests never run EM."""
import unittest
from pathlib import Path
from unittest.mock import patch
import gdstk
import numpy as np
from geometry.taper import staged_tests as stage
from geometry.taper.layout import make_taper


class StagedTaperTests(unittest.TestCase):
    def test_straight_wire_family_changes_only_uniform_width(self):
        cases = {q['name']: q for q in stage.cases()}
        for i, width in enumerate(stage.WIRE_WIDTHS, start=3):
            item = cases[f'{i:02d}_t2_wire_w{width}']
            self.assertEqual(len(item['metal']), 5)
            wire = item['metal'][-1]
            np.testing.assert_allclose(wire.bounding_box(),
                [(stage.T2_X-width/2,-618),(stage.T2_X+width/2,-8)])
            self.assertAlmostEqual(wire.area(), width*610)
            self.assertFalse(gdstk.boolean(item['metal'][:4],stage.grounds(),'xor',precision=.001))

    def test_t1_complete_unit_is_exactly_retained(self):
        expected = make_taper(stage.T1).get_polygons()
        expected = [q.copy().translate(stage.T1_X,0) for q in expected]
        box = stage.rect(stage.T1_X-stage.T1.w/2,-618,stage.T1_X+stage.T1.w/2,0)
        for item in stage.cases():
            if item['name'] in ('02_t1_alone','10_t2_t1_with_idc'):
                actual = gdstk.boolean(item['metal'],[box],'and',precision=.001)
                self.assertFalse(gdstk.boolean(actual,expected,'xor',precision=.001))

    def test_all_source_vertices_on_grid_and_expected_connectivity(self):
        counts = {'01_t1_inductor_only':1, '02_t1_alone':1,
                  '03_t2_wire_w16':1, '04_t2_wire_w32':1, '05_t2_wire_w64':1,
                  '06_t2_cap_only':2, '07_t2_alone_w32':1,
                  '08_input_idc_only':3, '09_t2_t1_coupling_only':3,
                  '10_t2_t1_with_idc':2}
        for item in stage.cases():
            with self.subTest(case=item['name']):
                for q in item['metal']:
                    np.testing.assert_allclose(q.points/2,np.round(q.points/2))
                self.assertEqual(len(gdstk.boolean(item['metal'],[],'or',precision=.001)),counts[item['name']])

    def test_preparation_rejects_em_and_non_read_soncmd(self):
        with patch.object(stage,'run') as runner:
            for command in ([stage.BIN/'em.exe','-v'],[stage.BIN/'soncmd.exe','-Analyze']):
                with self.assertRaisesRegex(ValueError,'no simulations'):
                    stage.native_read_only(command,Path('unused.log'))
            runner.assert_not_called()


if __name__=='__main__': unittest.main()
