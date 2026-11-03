import unittest
import numpy as np
import gdstk
from geometry.taper.progressive_fitting import transmission
from geometry.taper.staged_fitting import circuit
from geometry.taper.progressive_fixtures import build,CASES


class ProgressiveExtractionTests(unittest.TestCase):
    def test_scalar_transmission_matches_matrix_cascade(self):
        f=np.linspace(4.65e9,5.4e9,111)
        for ls,cs,gs in (([.74,6.014],[100,101],[459,56.7]),
                         ([6.014,4.923],[101,190.4],[56.7,8.2865]),
                         ([.74,6.014,4.923],[0,101,190.4],[None,56.7,8.2865])):
            np.testing.assert_allclose(transmission(f,ls,cs,gs,7),
                                       circuit(f,ls,cs,gs,left_ps=7)[:,1,0],rtol=1e-10,atol=1e-11)

    def test_progressive_fixture_removes_extra_end_cell(self):
        for name in CASES:
            polys,meta=build(*CASES[name])
            self.assertFalse(meta['extra_output_end_cell'])
            self.assertTrue(meta['output_feed_direct'])
            x=meta['w_sim']
            self.assertEqual(gdstk.inside([(x-1,-13),(x-1,-30),(x-1,10)],polys),(True,True,True))
            # The output feed is isolated from both nearby ground planes.
            self.assertEqual(gdstk.inside([(x-1,-4),(x-1,-22)],polys),(False,False))


if __name__=='__main__':unittest.main()
