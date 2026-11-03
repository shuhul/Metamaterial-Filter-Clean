"""Independent circuit-boundary and internal-load inversion checks."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
import unittest
import numpy as np
from geometry.taper.optimize_response import theory,FIXED
from geometry.taper.required_load import required_gamma

class ResponseObjectiveTests(unittest.TestCase):
    def test_nominal_schematic_includes_the_port_loaded_end_LC(self):
        f=np.linspace(3e9,8e9,101);expected=[]
        for v in f:
            w=2*np.pi*v;m=np.eye(2,dtype=complex)
            for cg,ell,cap in zip([459,56.7,8.2865,8.2865],
                                 [.74,6.06,4.9194,4.9194],[100,100,190.4,190.4]):
                m=m@np.array([[1,1/(1j*w*cg*1e-15)],[0,1]])
                m=m@np.array([[1,0],[1j*w*cap*1e-15+1/(1j*w*ell*1e-9),1]])
            a,b,c,d=m.ravel();expected.append(2/(a+b/50+c*50+d))
        expected=np.array(expected)*np.exp(-2j*np.pi*f*FIXED['effective_fixture_delay_ps']*1e-12)
        np.testing.assert_allclose(theory(f)[:,1,0],expected,atol=1e-11,rtol=1e-10)

    def test_internal_required_load_recovers_independently_terminated_network(self):
        rng=np.random.default_rng(75);mat=rng.normal(size=(6,6))
        z=np.tile((mat@mat.T*10+np.eye(6)*20+1j*(mat+mat.T)*8)[None,:,:],(4,1,1))
        # Form a three-port independently by eliminating shorted gap voltages.
        reduced=z[:,:3,:3]-z[:,:3,3:]@np.linalg.solve(z[:,3:,3:],z[:,3:,:3])
        eye=np.tile(np.eye(3),(4,1,1));s=(reduced-50*eye)@np.linalg.inv(reduced+50*eye)
        gamma=np.array([-.9,.1+.2j,.3-.2j,-.4+.1j])
        target=s[:,1,0]+s[:,1,2]*gamma*s[:,2,0]/(1-s[:,2,2]*gamma)
        np.testing.assert_allclose(required_gamma(z,target),gamma,atol=1e-11)

if __name__=='__main__':unittest.main()
