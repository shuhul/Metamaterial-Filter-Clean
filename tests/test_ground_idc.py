"""Ground-IDC trials must preserve isolation, the grid and downstream metal."""
import unittest
import numpy as np
import gdstk
from geometry.meander import MeanderParams
from geometry.taper.layout import TaperParams
from geometry.taper.second_layout import SecondParams
from geometry.taper.locked_taper2 import T1
from geometry.taper.compact_stages import build
from geometry.taper.ground_idc import GroundIDCParams,apply


class GroundIDCTests(unittest.TestCase):
    def layout(self):
        p=SecondParams(taper=TaperParams(w=500,t_top=10,t_left=24,t_right=20,l_left=440,l_right=548,
            meander=MeanderParams(w=160,l=600,l_top=10,n=4,t=52,d=4)),taper1=T1,
            straight_inductor=True,idc_pairs=26,idc_span=174)
        original,meta,_=build(p,'idc');metal,result=apply(original,meta,GroundIDCParams())
        return original,metal,meta,result

    def test_ground_capacitor_does_not_short_T2_or_change_input_isolation(self):
        old,metal,meta,_=self.layout()
        count=lambda q:len(gdstk.boolean(q,[],'or',precision=.001))
        self.assertEqual(count(old),count(metal))
        x=meta['cell_centers'][1]
        wire=gdstk.rectangle((x-26,-618),(x+26,-18))
        self.assertEqual(count(gdstk.boolean(metal,[wire],'not',precision=.001)),count(metal)+1)

    def test_T1_input_IDC_and_grid_are_preserved(self):
        original,metal,meta,_=self.layout()
        left=meta['cell_centers'][1]-250;left1=meta['cell_centers'][2]-T1.w/2
        for window in (gdstk.rectangle((0,-1052),(left,448)),
                       gdstk.rectangle((left1,-1052),(meta['w_sim'],448))):
            self.assertFalse(gdstk.boolean(gdstk.boolean(original,[window],'and'),
                gdstk.boolean(metal,[window],'and'),'xor',precision=.001))
        # Existing downstream GDS leads are snapped by the standard native
        # importer. This feature must add/remove only exact 2 um grid edges.
        delta=gdstk.boolean(original,metal,'xor',precision=.001)
        for q in delta:np.testing.assert_allclose(q.points/2,np.round(q.points/2),atol=1e-9,rtol=0)
