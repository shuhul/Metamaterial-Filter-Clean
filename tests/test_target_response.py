"""Check that acceptance uses specified electrical values, never fitted Cs."""
import unittest
import numpy as np
from geometry.taper.target_iteration import target_response
from geometry.taper.staged_fitting import circuit
from geometry.taper.target_iteration import parameters
from geometry.taper.clearance_extraction import build
from geometry.taper.target_heuristic import fast_network
import gdstk


class TargetResponseTests(unittest.TestCase):
    def test_matches_independent_matrix_cascade(self):
        f=np.linspace(4.65e9,5.4e9,1001)
        expected=circuit(f,[.74,6.06,4.9194,4.9194],
                         [100,100,190.4,190.4],[459,56.7,8.2865,8.2865],
                         left_ps=3,right_ps=4)[:,1,0]
        np.testing.assert_allclose(target_response(f,7),expected,rtol=1e-10,atol=1e-12)

    def test_delay_cannot_change_magnitude(self):
        f=np.linspace(4.65e9,5.4e9,1001)
        np.testing.assert_allclose(abs(target_response(f,80)),abs(target_response(f,0)),atol=1e-14)

    def test_fast_network_matches_matrix_model(self):
        f=np.linspace(4.65e9,5.4e9,1001)
        args=([.74,6.06,4.9194,4.9194],[174,121,190.4,190.4],
              [723,34,8.2865,8.2865])
        np.testing.assert_allclose(fast_network(f,*args,left_ps=6,right_ps=.7),
            circuit(f,*args,left_ps=6,right_ps=.7),rtol=1e-10,atol=1e-12)

    def test_equal_bus_ends_and_top_finger_clearance(self):
        p=parameters();polys,meta=build(p,'idc',200)
        left=meta['cell_centers'][1]-p.taper.w/2
        source=left-p.idc_span
        # Both buses extend to y=-592 and clear the flat bottom ground.
        self.assertEqual(gdstk.inside([(source+1,-591),(left+1,-591),
                                      (source+1,-593),(left+1,-593)],polys),
                         (True,True,False,False))
        self.assertEqual(p.taper.l_right,p.taper1.l_left)
        # At the first finger, horizontal source metal begins 2 um below
        # the node top bar and remains separated by a 4 um tip clearance.
        self.assertEqual(gdstk.inside([(source+40,-21),(source+40,-19),
                                      (left-2,-21)],polys),(True,False,False))


if __name__=='__main__':unittest.main()
