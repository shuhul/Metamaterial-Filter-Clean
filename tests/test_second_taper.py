"""Geometry/connectivity checks for the new IDC and full three-stage fixture."""
import unittest
from dataclasses import replace
from tempfile import TemporaryDirectory
from pathlib import Path
import json
import gdstk
import numpy as np
from geometry.taper.second_layout import SecondParams,build,write,validate,require_same_isolation_geometry
from geometry.taper.second_experiment import model_response


class SecondTaperTests(unittest.TestCase):
    def test_straight_trace_clearance_uses_actual_width_and_even_edges(self):
        from geometry.meander import MeanderParams
        taper=replace(SecondParams().taper,w=400,
                      meander=MeanderParams(w=52,t=16,n=3,l=596,d=6))
        p=SecondParams(taper=taper,straight_inductor=True,x_meander=-160)
        top,meta=build(p)
        self.assertEqual(len(gdstk.boolean(top.get_polygons(),[],'or',precision=.001)),2)
        wire=top.references[0].cell.polygons[0]
        self.assertTrue(np.allclose(wire.points/2,np.round(wire.points/2)))
        self.assertEqual(wire.bounding_box()[1][0]-wire.bounding_box()[0][0],16)
        with self.assertRaisesRegex(ValueError,'touches a cap leg'):
            replace(p,x_meander=-170)
        with self.assertRaisesRegex(ValueError,'multiple of 4'):
            replace(p,taper=replace(taper,meander=replace(taper.meander,t=6)))

    def test_calibration_reuse_requires_identical_isolation_metal(self):
        from geometry.taper.layout import TaperParams
        p=SecondParams(taper1=TaperParams(w=328,l_right=310))
        with TemporaryDirectory() as folder:
            gds=write(p,Path(folder)/'reference',isolation=True)
            require_same_isolation_geometry(replace(p,taper1=replace(p.taper1,l_right=240)),gds)
            with self.assertRaisesRegex(ValueError,'different isolation metal'):
                require_same_isolation_geometry(replace(p,idc_span=p.idc_span+4),gds)

    def test_two_um_junction_strip_keeps_grid_and_output_ground(self):
        p=SecondParams(t_junction_ground=2,d_t1_regular=48)
        top,meta=build(p);polys=top.get_polygons()
        t1_left=426+p.idc_span+p.taper.w+p.d_t2_t1
        edge=t1_left+328+24
        self.assertEqual(gdstk.inside([(edge-1,-400),(edge+1,-400),(edge+3,-400)],polys),
                         (False,True,False))
        self.assertEqual(len(gdstk.boolean(polys,[],'or',precision=.001)),2)
        self.assertEqual(edge%2,0)
        self.assertEqual((edge+2)%2,0)

    def test_upper_plate_and_ground_rail_keep_two_um_gaps(self):
        p=SecondParams(w_upper_plate=276,w_t1_ground_slot=264)
        top,_=build(p);polys=top.get_polygons()
        x=426+p.idc_span+p.taper.w+p.d_t2_t1+100
        self.assertEqual(gdstk.inside([(x,-1),(x,-3),(x,-7),(x,-9),(x,1),(x,10),(x,103)],polys),
                         (False,True,False,True,True,False,True))
        self.assertEqual(len(gdstk.boolean(polys,[],'or',precision=.001)),2)

    def test_coupling_teeth_use_open_gaps_without_touching_t1_runs(self):
        from geometry.taper.layout import TaperParams
        from geometry.meander import MeanderParams
        t1=TaperParams(w=280,t_left=24,t_right=6,l_left=200,l_right=310,
                      meander=MeanderParams(w=152,n=72,t=2,d=6))
        p=SecondParams(idc_pairs=15,t_coupling_finger=32,n_coupling_teeth=4,taper1=t1,
                       taper=replace(SecondParams().taper,l_left=300,l_right=598))
        top,_=build(p);polys=top.get_polygons()
        x=426+p.idc_span+p.taper.w+p.d_t2_t1
        self.assertEqual(gdstk.inside([(x+100,-53),(x+100,-51),(x+100,-55),(x+100,-49)],polys),
                         (True,False,False,True))
        self.assertEqual(len(gdstk.boolean(polys,[],'or',precision=.001)),2)

    def test_driven_finger_clears_cap_coil_and_bottom_ground(self):
        p=SecondParams(idc_pairs=15,t_coupling_finger=8,
                       taper=replace(SecondParams().taper,l_left=300,l_right=598))
        top,meta=build(p);polys=top.get_polygons()
        # t1 left = first input left + IDC span + t2 cap width + gap.
        x=426+p.idc_span+p.taper.w+p.d_t2_t1
        self.assertEqual(gdstk.inside([(x+34,-300),(x+26,-300),(x+42,-300)],polys),(True,False,False))
        self.assertEqual(gdstk.inside([(x+12,-604),(x+12,-598),(x+12,-610)],polys),(True,False,False))
        self.assertEqual(len(gdstk.boolean(polys,[],'or',precision=.001)),2)

    def test_idc_clearances_and_ground_connectivity(self):
        p=SecondParams();top,meta=build(p);polys=top.get_polygons()
        self.assertEqual(len(gdstk.boolean(polys,[],'or',precision=.001)),2)
        # First input tooth, gap, then first taper tooth.
        x=600
        self.assertEqual(gdstk.inside([(x,-34),(x,-38),(x,-42)],polys),(True,False,True))
        self.assertFalse(gdstk.inside([(752,-34)],polys)[0])
        self.assertTrue(all(gdstk.inside([(x,1) for x in range(2,int(meta['w_sim'])-2,50)],polys)))
    def test_calibration_keeps_taper_teeth_and_removes_feed_teeth(self):
        top,meta=build(SecondParams(),isolation=True)
        self.assertEqual(len(gdstk.boolean(top.get_polygons(),[],'or',precision=.001)),4)
        self.assertEqual(gdstk.inside([(600,-34),(600,-42)],top.get_polygons()),(False,True))
    def test_metadata_detects_stale_comb(self):
        with TemporaryDirectory() as d:
            path=write(SecondParams(),Path(d)/'fixture')
            meta=json.loads(path.with_suffix('.json').read_text());validate(path,meta)
            meta['second']['idc_pairs']-=1
            with self.assertRaises(ValueError):validate(path,meta)
    def test_three_stage_reciprocity_and_delay_magnitude(self):
        f=np.linspace(4e9,18e9,401)
        x=np.array([100,101.49,188.45,459,56.7,8.2865,8.2865,0])
        s=model_response(f,x,.74)
        np.testing.assert_allclose(s[:,0,1],s[:,1,0])
        np.testing.assert_allclose(np.swapaxes(s.conj(),1,2)@s,np.broadcast_to(np.eye(2),s.shape),atol=1e-10)
        x[-1]=5
        np.testing.assert_allclose(abs(s[:,1,0]),abs(model_response(f,x,.74)[:,1,0]),atol=1e-14)

    def test_joint_complex_fit_recovers_strong_input_and_both_delays(self):
        import contextlib,io
        import skrf as rf
        from geometry.taper.second_experiment import fit_network
        f=np.linspace(3.8e9,9e9,1001)
        expected=[100,101.49,188.45,459,56.7,8.2865,8.2865,4.0,.6]
        s=model_response(f,expected,.74,l1=6.08)
        with TemporaryDirectory() as d,contextlib.redirect_stdout(io.StringIO()):
            path=Path(d)/'synthetic'
            rf.Network(frequency=rf.Frequency.from_f(f,unit='Hz'),s=s,z0=50).write_touchstone(str(path))
            result=fit_network(path.with_suffix('.s2p'),Path(d),.74,l1=6.08)
        np.testing.assert_allclose(result['x'],expected,rtol=1e-5,atol=1e-5)


if __name__=='__main__':unittest.main()
