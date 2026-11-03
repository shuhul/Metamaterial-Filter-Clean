"""Vertical comb isolation, common upper node, and unique feed-port edge."""
import unittest
import gdstk
from geometry.taper.vertical_idc_case import parameters,build_vertical

class VerticalIDCTests(unittest.TestCase):
    def test_connected_comb_has_unique_feed_and_no_lower_node_bus(self):
        parts,meta=build_vertical(parameters())
        self.assertEqual(len(gdstk.boolean(parts,[],'or',precision=.001)),2)
        self.assertEqual(sum(gdstk.inside([(.5,-13)],[q])[0] for q in parts),1)
        self.assertTrue(gdstk.inside([(610,-10)],parts)[0])
        self.assertFalse(any(gdstk.inside([(610,-18),(610,-300),(448,-10),(452,-486),(460,-14)],parts)))
        self.assertTrue(all(q.layer==1 and q.datatype==0 for q in parts))

    def test_isolation_retains_separate_feed_conductors(self):
        parts,meta=build_vertical(parameters(),isolation=True)
        self.assertEqual(len(gdstk.boolean(parts,[],'or',precision=.001)),4)
        self.assertEqual(meta['parameters']['n_cells'],1)

if __name__=='__main__':unittest.main()
