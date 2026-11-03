"""Verify fixed-reference topology and fit gates used to bracket minimum L."""
import unittest
import numpy as np
from geometry.taper.minimum_inductance_review import PROFILES,predictor,assessment
from geometry.taper.staged_fitting import circuit,LT1,LR


class MinimumInductanceReviewTests(unittest.TestCase):
    def test_profiles_map_to_the_same_declared_circuit(self):
        f=np.linspace(4.65e9,5.4e9,101)
        for name,p in PROFILES.items():
            for l in (.74,5):
                expected=circuit(f,[l,LT1,LR,LR],[100,p['ct'],p['cr'],p['ce']],
                                 [None,56.7,p['gj'],p['go']],left_ps=6)[:,1,0]
                np.testing.assert_allclose(predictor(f,l,name)([100,56.7,6]),expected,rtol=1e-10,atol=1e-11)

    def test_identical_merged_peaks_do_not_have_invented_linewidths(self):
        f=np.arange(9,dtype=float)*1e6
        a=np.sqrt(np.array([0,.2,1,.7,.7,1,.2,0,0]))
        r=assessment(f,a,a,dict(complex_RMS=0,complex_max=0,active_bounds=[0,0,0]))
        self.assertTrue(r['pass_curve'])
        self.assertTrue(all(p['power_FWHM_MHz'] is None for p in r['observed_peaks']))

    def test_bad_error_or_bound_cannot_pass(self):
        f=np.arange(9,dtype=float)*1e6
        a=np.sqrt(np.array([0,.2,1,.7,.7,1,.2,0,0]))
        self.assertFalse(assessment(f,a,a,dict(complex_RMS=.02,complex_max=.04,active_bounds=[0,0,0]))['pass_curve'])
        self.assertFalse(assessment(f,a,a,dict(complex_RMS=0,complex_max=0,active_bounds=[-1,0,0]))['pass_curve'])


if __name__=='__main__':unittest.main()
