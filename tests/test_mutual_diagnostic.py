import unittest
import numpy as np
from geometry.taper.mutual_diagnostic import bare,loaded,full
from geometry.taper.staged_fitting import circuit,LT1,LR


class MutualDiagnosticTests(unittest.TestCase):
    def test_zero_mutual_reproduces_uncoupled_target_topologies(self):
        f=np.linspace(4e9,7e9,101)
        expected=circuit(f,[.74,LT1,LR],[0,100.69929340643347,190.4],[None,56.7,8.2865],left_ps=8)[:,1,0]
        np.testing.assert_allclose(bare(f,.74,56.7,0,8),expected,rtol=1e-10,atol=1e-11)
        expected=circuit(f,[.74,LT1],[900,100.69929340643347],[459,56.7],left_ps=8)[:,1,0]
        np.testing.assert_allclose(loaded(f,.74,100,459,56.7,0,8,add=800),expected,rtol=1e-10,atol=1e-11)
        for n in (1,4):
            expected=circuit(f,[.74,LT1]+[LR]*n,[100,100.69929340643347]+[190.4]*n,
                             [459,56.7]+[8.2865]*n,left_ps=8)[:,1,0]
            np.testing.assert_allclose(full(f,.74,100,459,56.7,0,8,n_regular=n),expected,rtol=1e-10,atol=1e-11)


if __name__=='__main__':unittest.main()
