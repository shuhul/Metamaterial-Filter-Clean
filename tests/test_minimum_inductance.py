"""Independent geometry and linewidth checks for the width-only controls."""
import unittest
import gdstk
import numpy as np
from geometry.taper.minimum_inductance import half_power_features
from geometry.taper.minimum_inductance_meander import build


class MinimumInductanceTests(unittest.TestCase):
    def test_absolute_half_power_width(self):
        f=np.linspace(4.8e9,5.2e9,40001)
        # Power Lorentzian has an independently specified 20 MHz FWHM.
        s=1/(1+2j*(f-5e9)/20e6)
        peaks=half_power_features(f,s)
        self.assertEqual(len(peaks),1)
        self.assertAlmostEqual(peaks[0]['power_FWHM_MHz'],20,places=6)

    def test_merged_half_power_region_has_no_individual_width(self):
        f=np.arange(9,dtype=float)*1e6
        power=np.array([0,.2,1,.7,.7,1,.2,0,0])
        peaks=half_power_features(f,np.sqrt(power))
        self.assertEqual(len(peaks),2)
        self.assertTrue(all(p['power_FWHM_MHz'] is None for p in peaks))

    def test_only_inductor_changes_with_width(self):
        a,ma,_=build(200,'direct');b,mb,_=build(80,'direct')
        delta=gdstk.boolean(a,b,'xor',precision=.001)
        x=ma['cell_centers'][1]
        allowed=gdstk.rectangle((x-100,-618),(x+100,-18))
        self.assertFalse(gdstk.boolean(delta,[allowed],'not',precision=.001))
        self.assertEqual(ma['second'],mb['second'] | {'taper':ma['second']['taper']})
        # Both endpoints remain connected to the cap and flat bottom ground.
        for polys in (a,b):
            self.assertTrue(all(gdstk.inside([(x,-19),(x,-617),(x,-619)],polys)))
            # A former 52 um straight wire must not remain between coil rows.
            self.assertFalse(gdstk.inside([(x+10,-35)],polys)[0])
        self.assertEqual(ma['meander']['n'],58)
        self.assertEqual(mb['meander']['n'],58)


if __name__=='__main__':unittest.main()
