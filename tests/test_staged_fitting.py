"""Independent nodal check of the directly driven taper circuit."""
import unittest
import numpy as np
from geometry.taper.staged_fitting import circuit


class StagedFittingTests(unittest.TestCase):
    def test_direct_feed_and_output_end_match_nodal_elimination(self):
        f=np.linspace(4.6e9,5.4e9,51);w=2*np.pi*f
        ls=[6.014,4.9227,4.9227];cs=[101.49,188.45,188.45];gs=[None,8.23,8.31]
        y=np.zeros((len(f),3,3),complex)
        for i,(l,c) in enumerate(zip(ls,cs)):
            y[:,i,i]=1j*w*c*1e-15+1/(1j*w*l*1e-9)
        for a,b,c in [(0,1,gs[1]),(1,2,gs[2])]:
            v=1j*w*c*1e-15
            y[:,a,a]+=v;y[:,b,b]+=v;y[:,a,b]-=v;y[:,b,a]-=v
        external=[0,2]
        yy=y[:,external,:][:,:,external]-y[:,external,1:2]@np.linalg.solve(y[:,1:2,1:2],y[:,1:2,external])
        identity=np.eye(2)[None,:,:]
        expected=(identity-50*yy)@np.linalg.inv(identity+50*yy)
        np.testing.assert_allclose(circuit(f,ls,cs,gs),expected,atol=2e-13)


if __name__=='__main__':unittest.main()
