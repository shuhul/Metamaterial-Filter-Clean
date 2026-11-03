"""Audit controls independently verify the exact loading and topology."""
import unittest
import numpy as np
from geometry.taper.ideal_inductor_audit import input_shunt,synthetic_operator
from geometry.taper.ideal_wire_control import load_inductor
from geometry.taper.ideal_topology_diagnostics import prediction
from geometry.taper.staged_fitting import circuit,LT1,CT1,LR,CR,GOUT
from geometry.taper.fitting import chain
from geometry.taper.locked_taper2 import FIXED
from geometry.taper.paired_coupling import abcd


class IdealInductorAuditTests(unittest.TestCase):
    def test_ideal_load_matches_independent_nodal_model(self):
        f=np.linspace(4.65e9,5.4e9,101)
        for l in (.74,2,5):
            actual=load_inductor(synthetic_operator(f),f,l)
            expected=circuit(f,[l,LT1,LR,LR],[130,CT1,CR,CR],[None,56.7,8.40642914296533,GOUT])
            np.testing.assert_allclose(actual,expected,rtol=1e-10,atol=1e-11)

    def test_external_shunt_formula_matches_abcd_multiplication(self):
        f=np.linspace(4.65e9,5.4e9,101)
        original=circuit(f,[LT1,LR],[CT1,CR],[56.7,8.4],left_ps=3,right_ps=2)
        m=abcd(original)
        y=1/(2j*np.pi*f*.74e-9)
        m[:,1,:]+=y[:,None]*m[:,0,:]
        a,b,c,d=m[:,0,0],m[:,0,1],m[:,1,0],m[:,1,1]
        den=a+b/50+50*c+d
        expected=np.empty_like(original)
        expected[:,0,0]=(a+b/50-50*c-d)/den
        expected[:,1,1]=(-a+b/50-50*c+d)/den
        expected[:,0,1]=expected[:,1,0]=2/den
        np.testing.assert_allclose(input_shunt(original,f,.74),expected,rtol=1e-10,atol=1e-11)

    def test_old_topology_reproduces_original_chain(self):
        f=np.linspace(4.65e9,5.4e9,101)
        actual=prediction(f,None,[0,56.7,0],'selected_original_topology')
        expected=chain(2*np.pi*f,[LT1,LR],[CT1,FIXED['C_regular_fF']],
                       [56.7,FIXED['Cg_junction_fF'],FIXED['Cg_output_fF']])[:,1,0]
        np.testing.assert_allclose(actual,expected,rtol=1e-10,atol=1e-11)


if __name__=='__main__':unittest.main()
