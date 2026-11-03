import unittest
from dataclasses import replace
import gdstk
from geometry.taper.bridged_fixtures import t2_parameters
from geometry.taper.compact_stages import build, prepare, JunctionIDC


class CompactStageTests(unittest.TestCase):
    def test_IDC_ideal_control_removes_only_the_full_straight_wire(self):
        p=replace(self.parameters(),w_ground_pad=24)
        actual,meta,_=build(p,'idc')
        absent,am,_=build(p,'idc',omit_wire=True)
        x=meta['cell_centers'][1];t=p.taper.meander.t
        expected=gdstk.boolean(actual,[gdstk.rectangle((x-t/2,-618),(x+t/2,-18))],
                               'not',precision=.001)
        self.assertFalse(gdstk.boolean(expected,absent,'xor',precision=.001))
        self.assertTrue(am['IDC_present'])
        self.assertTrue(am['wire_removed'])

    def test_two_ground_pads_preserve_T1_and_full_straight_wire(self):
        p=self.parameters()
        bare,meta,_=build(p,'idc')
        pads,pm,_=build(replace(p,w_ground_pad=24,l_ground_pad=586,d_ground_pad=12),'idc')
        x=meta['cell_centers'][1]
        left=x-p.taper.w/2+p.taper.t_left+12
        right=x+p.taper.w/2-p.taper.t_right-12
        expected=[gdstk.rectangle((left,-618),(left+24,-32)),
                  gdstk.rectangle((right-24,-618),(right,-32))]
        self.assertFalse(gdstk.boolean(gdstk.boolean(pads,bare,'xor',precision=.001),
                                      expected,'xor',precision=.001))
        # The straight wire keeps its original full height, with no pedestal.
        self.assertEqual(gdstk.inside([(x,-600),(x,-20),(x+28,-600)],pads),(True,True,False))
        t1left=meta['cell_centers'][2]-p.taper1.w/2
        region=[gdstk.rectangle((t1left,-1052),(meta['w_sim'],448))]
        self.assertFalse(gdstk.boolean(gdstk.boolean(bare,region,'and'),
                                      gdstk.boolean(pads,region,'and'),'xor',precision=.001))

    def test_rejected_top_electrode_cannot_be_prepared_or_simulated(self):
        with self.assertRaisesRegex(ValueError,'rejected by the user'):
            prepare('must_not_be_created',replace(self.parameters(),w_upper_plate=300))

    def test_fast_direct_predictor_preserves_complex_circuit_response(self):
        import numpy as np
        from geometry.taper.compact_fit import direct_predictor,GJ
        from geometry.taper.staged_fitting import circuit,LT1,CT1,LR,CR,GOUT
        f=np.linspace(4.7e9,5.4e9,71)
        a=direct_predictor(f,.734)([64,56.9,4.2])
        b=circuit(f,[.734,LT1,LR,LR],[64,CT1,CR,CR],[None,56.9,GJ,GOUT],left_ps=4.2)[:,1,0]
        np.testing.assert_allclose(a,b,rtol=1e-10,atol=1e-12)

    def parameters(self):
        p=t2_parameters()
        return replace(p,taper=replace(p.taper,w=500,t_top=10,
            meander=replace(p.taper.meander,l=600,w=160,t=52)),idc_span=160)

    def test_direct_has_only_feed_wire_and_right_coupler(self):
        p=self.parameters();polys,meta,_=build(p)
        x=meta['cell_centers'][1];left=x-250
        self.assertEqual(gdstk.inside([(left+12,-300),(left-30,-100),
            (left-30,-13),(left+100,-13),(left+100,-19),
            (x,-600),(x-28,-600),(left+490,-300)],polys),
            (False,False,True,True,False,True,False,True))
        self.assertFalse(meta['IDC_present'])

    def test_Lcal_has_weak_input_gap_and_identical_T2_metal(self):
        p=self.parameters();direct,dm,_=build(p);cal,cm,_=build(p,'Lcal')
        x=dm['cell_centers'][1];left=x-250
        window=gdstk.rectangle((left,-618),(left+500,-8))
        a=gdstk.boolean(direct,[window],'and',precision=.001)
        b=gdstk.boolean(cal,[window],'and',precision=.001)
        self.assertFalse(gdstk.boolean(a,b,'xor',precision=.001))
        self.assertEqual(gdstk.inside([(left-4,-13),(left-10,-13)],cal),(False,True))
        self.assertFalse(cm['IDC_present'])

    def test_shield_loop_adds_only_ground_metal_beside_full_height_wire(self):
        p=self.parameters();bare,meta,_=build(p)
        shielded,_,_=build(replace(p,t_inner_ground=80,x_inner_ground=100,w_inner_ground_hole=72))
        x=meta['cell_centers'][1]
        outer=gdstk.rectangle((x+60,-618),(x+140,-22))
        hole=gdstk.rectangle((x+64,-614),(x+136,-26))
        expected=gdstk.boolean([outer],[hole],'not',precision=.001)
        delta=gdstk.boolean(shielded,bare,'xor',precision=.001)
        self.assertFalse(gdstk.boolean(delta,expected,'xor',precision=.001))
        self.assertEqual(gdstk.inside([(x+62,-300),(x+100,-300),(x+100,-24),
                                     (x,-600),(x+28,-600)],shielded),
                         (True,False,True,True,False))

    def test_junction_IDC_stays_in_gap_and_preserves_T1(self):
        p=self.parameters()
        p=replace(p,d_t2_t1=200,taper=replace(p.taper,l_right=108))
        bare,meta,_=build(p)
        teeth,tm,_=build(p,junction=JunctionIDC())
        left1=meta['cell_centers'][2]-p.taper1.w/2
        t1_window=gdstk.rectangle((left1,-618),(meta['w_sim'],0))
        a=gdstk.boolean(bare,[t1_window],'and',precision=.001)
        b=gdstk.boolean(teeth,[t1_window],'and',precision=.001)
        self.assertFalse(gdstk.boolean(a,b,'xor',precision=.001))
        self.assertFalse(tm['input_IDC_present'])
        self.assertTrue(tm['junction_IDC_present'])
        # Two alternating combs remain electrically separate in this gap.
        self.assertEqual(len(gdstk.boolean(teeth,[],'or',precision=.001)),
                         len(gdstk.boolean(bare,[],'or',precision=.001)))

    def test_top_electrode_adds_metal_only_in_empty_gap(self):
        p=self.parameters();p=replace(p,taper=replace(p.taper,l_right=12))
        bare,meta,_=build(p)
        added,_,_=build(replace(p,w_upper_plate=300))
        self.assertFalse(gdstk.boolean(bare,added,'not',precision=.001))
        edge=meta['cell_centers'][1]+p.taper.w/2
        left1=meta['cell_centers'][2]-p.taper1.w/2
        expected=[gdstk.rectangle((edge-4,-8),(edge,-2)),
                  gdstk.rectangle((edge-4,-6),(left1+300,-2))]
        difference=gdstk.boolean(added,bare,'xor',precision=.001)
        self.assertFalse(gdstk.boolean(difference,expected,'xor',precision=.001))

    def test_IDC_calibration_preserves_all_node_comb_and_wire_metal(self):
        p=self.parameters();p=replace(p,idc_pairs=30,idc_span=172,
                                      taper=replace(p.taper,l_left=504))
        full,meta,_=build(p,'idc');cal,cm,_=build(p,'idc_Lcal')
        left=meta['cell_centers'][1]-p.taper.w/2;input_left=left-p.idc_span
        mask=[gdstk.rectangle((left,-618),(left+p.taper.w,-8))]
        first=-8-p.taper.t_top-p.d_idc_top
        for row in range(1,2*p.idc_pairs,2):
            y=first-row*(p.idc_t+p.idc_d)
            mask.append(gdstk.rectangle((input_left+20+p.idc_d,y-p.idc_t),
                                        (left+p.taper.t_left,y)))
        self.assertFalse(gdstk.boolean(gdstk.boolean(full,mask,'and'),
                                      gdstk.boolean(cal,mask,'and'),'xor',precision=.001))
        self.assertTrue(meta['IDC_present'])
        self.assertFalse(cm['IDC_present'])
        self.assertTrue(cm['node_IDC_comb_present'])
        # Input comb removed, node comb retained in the calibration.
        self.assertEqual(gdstk.inside([(input_left+60,first-2),
                                      (input_left+60,first-10)],cal),(False,True))


if __name__=='__main__':unittest.main()
