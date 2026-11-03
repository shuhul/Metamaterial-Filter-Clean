"""Geometry invariants that protect the accepted T1 during T2 iteration."""
import unittest
from dataclasses import replace,asdict
import gdstk
from geometry.taper.taper2_v3 import BASE,T1,InputCoupler,build,require_locked


class Taper2GeometryTests(unittest.TestCase):
    def test_junction_probe_has_only_t1_and_the_requested_feed_rail(self):
        top,meta=build(T1,isolation=True,input_probe=True,coupler=InputCoupler(d=168,t_input=64))
        self.assertEqual(meta['unit'],asdict(T1))
        self.assertEqual(len(top.references),1)
        self.assertEqual(meta['expected_metal_components'],3)
        edge=meta['cell_lefts'][0]+T1.cap.w
        self.assertEqual(gdstk.inside([(edge-63,-300),(edge-65,-300),(edge+1,-300)],top.get_polygons()),
                         (True,False,False))

    def test_plain_and_idc_preserve_t1_and_downstream(self):
        for coupling in (InputCoupler(),InputCoupler(d=28,n=72),InputCoupler(d=28,n=72,d_t1=96)):
            with self.subTest(coupling=coupling):
                top,meta=build(BASE,coupler=coupling)
                polygons=top.get_polygons();require_locked(polygons,meta)
                self.assertEqual(len(gdstk.boolean(polygons,[],'or',precision=.001)),2)
                upper=gdstk.rectangle((0,0),(meta['w_sim'],448))
                self.assertFalse(gdstk.boolean([upper],polygons,'not',precision=.001))
                edge=meta['cell_lefts'][1]+BASE.cap.w
                opening=gdstk.rectangle((edge,-618),(meta['cell_lefts'][2],0))
                self.assertFalse(gdstk.boolean(polygons,[opening],'and',precision=.001))

    def test_idc_is_uniform_and_fits_equal_rails(self):
        coupler=InputCoupler(d=28,n=72)
        top,meta=build(BASE,coupler=coupler)
        edge=meta['cell_lefts'][0]+BASE.cap.w
        teeth=[p for p in top.polygons if p.bounding_box()[0][0]>=edge and p.bounding_box()[1][0]<=edge+28]
        self.assertEqual(len(teeth),72)
        self.assertTrue(all(abs(p.area()-26*4)<1e-8 for p in teeth))
        upper=-BASE.d_top-BASE.cap.t_top-2
        self.assertEqual(min(p.bounding_box()[0][1] for p in teeth),upper-71*8-4)
        self.assertEqual(max(p.bounding_box()[1][1] for p in teeth),upper)

    def test_rejects_shorted_or_off_grid_idc(self):
        for coupling in (InputCoupler(d=3),InputCoupler(d=4,n=72),InputCoupler(d=28,n=74)):
            with self.subTest(coupling=coupling),self.assertRaises(ValueError):build(BASE,coupler=coupling)

    def test_isolated_and_response_share_same_unit(self):
        _,cal=build(BASE,isolation=True);_,full=build(BASE)
        self.assertEqual(cal['unit'],full['unit'])

    def test_declared_junction_adjustment_preserves_cell_shapes(self):
        for gap,ground in ((64,4),(62,2),(68,8)):
            top,meta=build(BASE,coupler=InputCoupler(d=34,n=72,d_regular=gap,t_regular_ground=ground))
            require_locked(top.get_polygons(),meta)
            edge=meta['cell_lefts'][2]+meta['taper1']['cap']['w']
            mid=edge+gap/2
            self.assertEqual((gap-ground)/2,30)
            self.assertEqual(gdstk.inside([(mid,-300),(mid-ground/2-1,-300),
                                          (mid+ground/2+1,-300)],top.get_polygons()),(True,False,False))


if __name__=='__main__':unittest.main()
