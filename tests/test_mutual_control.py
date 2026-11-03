"""Check coupled nodal model against the independent uncoupled ABCD ladder."""
import unittest
import numpy as np
from geometry.taper.mutual_control import coupled
from geometry.taper.staged_fitting import circuit


class MutualControlTests(unittest.TestCase):
    def test_zero_mutual_matches_abcd_with_and_without_input_cap(self):
        f=np.linspace(4.6e9,5.4e9,201)
        ls=[.74,6.014,4.9227,4.9227];cs=[100,101.49,188.878,188.878];gs=[56.7,8.406,8.242]
        for g in (None,459):
            expected=circuit(f,ls,cs,[g]+gs,left_ps=4,right_ps=1)
            np.testing.assert_allclose(coupled(f,ls,cs,gs,input_cg=g,left_ps=4,right_ps=1),expected,atol=3e-12)

    def test_coupled_lossless_network_is_reciprocal_and_unitary(self):
        f=np.linspace(4.5e9,5.5e9,101)
        s=coupled(f,[.74,6,4.92,4.92],[100,101,189,189],[56.7,8.4,8.2],mutual_nh=.2)
        np.testing.assert_allclose(s[:,0,1],s[:,1,0],atol=1e-13)
        np.testing.assert_allclose(np.swapaxes(s.conj(),1,2)@s,np.broadcast_to(np.eye(2),s.shape),atol=1e-12)


if __name__=='__main__':unittest.main()
