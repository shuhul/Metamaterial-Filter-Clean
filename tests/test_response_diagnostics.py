import unittest
import numpy as np
from geometry.taper.response_diagnostics import local_pole


class PoleDiagnosticsTests(unittest.TestCase):
    def test_recovers_shared_causal_pole_with_different_complex_backgrounds(self):
        # Known analytic pole with the exp(+jwt) convention, independent of
        # the fitting code. Different backgrounds mimic three S channels.
        f=np.linspace(4.83e9,4.86e9,601)
        frequency=4.84403e9;gamma=1.24e6
        residue=np.array([.8+.1j,.4-.2j,-.3+.7j])*1e6
        background=np.array([.4+.2j,.02-.01j,-.6+.1j])
        y=background+residue/(f[:,None]-frequency-1j*gamma)
        result=local_pole(f,y,4.8441,2.5,15)
        self.assertTrue(result['optimizer_success'])
        self.assertAlmostEqual(result['frequency_GHz'],frequency/1e9,places=8)
        self.assertAlmostEqual(result['loaded_linewidth_MHz'],2*gamma/1e6,places=6)
        self.assertLess(result['local_complex_RMS'],1e-8)


if __name__=='__main__':unittest.main()
