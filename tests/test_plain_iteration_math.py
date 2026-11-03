"""Independent nodal checks of the plain-geometry diagnostic algebra."""
import unittest
import numpy as np
from geometry.taper.ideal_load_search import terminate
from geometry.taper.plain_series import change_caps,A
from geometry.taper.idc_series_control import change as change_idc
from geometry.taper.ideal_wire_control import load_inductor


class PlainIterationMathTests(unittest.TestCase):
    def test_floating_inductor_load_matches_direct_nodal_reduction(self):
        f=np.array([4.7,4.9,5.1,5.3])*1e9;eye=np.eye(4)[None]
        g=np.eye(4)*.001;c=np.eye(4)*30e-15
        for i,j,r in ((0,2,50),(1,3,80),(0,1,400),(2,3,600)):
            g[i,i]+=1/r;g[j,j]+=1/r;g[i,j]-=1/r;g[j,i]-=1/r
        y=g[None]+2j*np.pi*f[:,None,None]*c[None]
        s=(eye-50*y)@np.linalg.inv(eye+50*y);loaded=y.copy()
        v=1/(2j*np.pi*f*.74e-9)
        loaded[:,2,2]+=v;loaded[:,3,3]+=v;loaded[:,2,3]-=v;loaded[:,3,2]-=v
        yy=loaded[:,:2,:2]-loaded[:,:2,2:]@np.linalg.solve(loaded[:,2:,2:],loaded[:,2:,:2])
        small=np.eye(2)[None];expected=(small-50*yy)@np.linalg.inv(small+50*yy)
        np.testing.assert_allclose(load_inductor(s,f,.74),expected,rtol=1e-11,atol=1e-12)

    def test_five_conductor_idc_control_matches_independent_nodal_circuit(self):
        f=np.array([4.7,4.9,5.1,5.3])*1e9;w=2*np.pi*f
        cap=np.diag([85.,130.,152.,274.,300.])
        for i,j,c in ((0,1,403),(1,2,69),(2,3,18),(3,4,15),(0,2,.7)):
            cap[i,i]+=c;cap[j,j]+=c;cap[i,j]-=c;cap[j,i]-=c
        b=np.zeros((5,6))
        for port,node in enumerate([0,4,1,2,3,4]):b[node,port]=1
        ell=np.diag([.74,6,4.92,4.92])*1e-9
        ell[0,1]=ell[1,0]=.1e-9
        z=(b.T@np.linalg.inv(cap*1e-15)@b)[None]/(1j*w[:,None,None])
        z[:,2:,2:]+=1j*w[:,None,None]*ell[None]
        changed=change_idc(z,f,cap,[110,450,80])
        expected_c=cap.copy();expected_c[1,1]-=20
        for i,j,delta in ((0,1,47),(1,2,11)):
            expected_c[i,i]+=delta;expected_c[j,j]+=delta;expected_c[i,j]-=delta;expected_c[j,i]-=delta
        y=1j*w[:,None,None]*expected_c[None]*1e-15
        y[:,1:,1:]+=np.linalg.inv(ell)[None]/(1j*w[:,None,None])
        ports=[0,4];internal=[1,2,3]
        yy=y[:,ports,:][:,:,ports]-y[:,ports,:][:,:,internal]@np.linalg.solve(y[:,internal,:][:,:,internal],y[:,internal,:][:,:,ports])
        eye=np.eye(2)[None];expected=(eye-50*yy)@np.linalg.inv(eye+50*yy)
        np.testing.assert_allclose(A['short_internal_z'](changed),expected,rtol=1e-9,atol=1e-10)

    def test_multiple_floating_capacitors_match_direct_nodal_loading(self):
        f=np.array([4.65,4.9,5.1,5.4])*1e9;eye=np.eye(8)[None]
        conductance=np.eye(8)*.001
        capacitance=np.eye(8)*20e-15
        for i,j,r,c in ((0,2,300,25),(1,6,400,30),(2,4,800,40),(3,5,500,35),
                        (4,7,900,15),(5,6,700,20),(0,3,600,10),(1,7,450,15)):
            for matrix,value in ((conductance,1/r),(capacitance,c*1e-15)):
                matrix[i,i]+=value;matrix[j,j]+=value;matrix[i,j]-=value;matrix[j,i]-=value
        y=conductance[None]+2j*np.pi*f[:,None,None]*capacitance[None]
        s=(eye-50*y)@np.linalg.inv(eye+50*y)
        caps=[900,0,17]
        loaded=y.copy()
        for k,c in enumerate(caps):
            a,b=2+2*k,3+2*k;v=2j*np.pi*f*c*1e-15
            loaded[:,a,a]+=v;loaded[:,b,b]+=v;loaded[:,a,b]-=v;loaded[:,b,a]-=v
        effective=loaded[:,:2,:2]-loaded[:,:2,2:]@np.linalg.solve(loaded[:,2:,2:],loaded[:,2:,:2])
        small=np.eye(2)[None];expected=(small-50*effective)@np.linalg.inv(small+50*effective)
        np.testing.assert_allclose(terminate(s,f,caps),expected,rtol=1e-11,atol=1e-12)

    def test_plain_two_cap_changes_match_independent_coupled_nodal_circuit(self):
        f=np.array([4.7,4.9,5.1,5.3])*1e9;w=2*np.pi*f
        cap=np.diag([130.,152.,274.,300.])
        for i,j,c in ((0,1,69),(1,2,18),(2,3,15),(0,2,1.3)):
            cap[i,i]+=c;cap[j,j]+=c;cap[i,j]-=c;cap[j,i]-=c
        mapping=np.zeros((4,6))
        for port,node in enumerate([0,3,0,1,2,3]):mapping[node,port]=1
        inductance=np.diag([.74,6,4.92,4.92])*1e-9
        inductance[0,1]=inductance[1,0]=.1e-9
        inductance[0,2]=inductance[2,0]=.019e-9
        z=(mapping.T@np.linalg.inv(cap*1e-15)@mapping)[None]/(1j*w[:,None,None])
        z[:,2:,2:]+=1j*w[:,None,None]*inductance[None]
        changed=change_caps(z,f,cap,[110,80])
        expected_c=cap.copy();expected_c[0,0]-=20
        expected_c[0,0]+=11;expected_c[1,1]+=11;expected_c[0,1]-=11;expected_c[1,0]-=11
        y=1j*w[:,None,None]*expected_c[None]*1e-15+np.linalg.inv(inductance)[None]/(1j*w[:,None,None])
        ports=[0,3];internal=[1,2]
        effective=y[:,ports,:][:,:,ports]-y[:,ports,:][:,:,internal]@np.linalg.solve(y[:,internal,:][:,:,internal],y[:,internal,:][:,:,ports])
        eye=np.eye(2)[None];expected=(eye-50*effective)@np.linalg.inv(eye+50*effective)
        np.testing.assert_allclose(A['short_internal_z'](changed),expected,rtol=1e-9,atol=1e-10)


if __name__=='__main__':unittest.main()
