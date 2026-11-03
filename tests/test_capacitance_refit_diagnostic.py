import unittest
import numpy as np
from geometry.taper.capacitance_refit_diagnostic import predict,fit_data
from geometry.taper.staged_fitting import circuit,LT1,LR


class CapacitanceRefitDiagnosticTests(unittest.TestCase):
    def test_predictor_matches_cascade_with_zero_T2_capacitance(self):
        f=np.linspace(4.65e9,5.4e9,101)
        x=[101,189,188,57,8.2,8.3,6]
        expected=circuit(f,[2.52,LT1,LR,LR],[0,*x[:3]],[None,*x[3:6]],left_ps=x[6])[:,1,0]
        np.testing.assert_allclose(predict(f,2.52,x),expected,rtol=1e-10,atol=1e-11)


if __name__=='__main__':unittest.main()
