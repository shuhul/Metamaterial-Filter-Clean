"""Zero-EM ideal T2 embedding into the preserved measured T1 two-port."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks
from circuit_common import read_snp
from geometry.taper.locked_taper2 import VERIFIED
from geometry.taper.optimize_response import OUT,theory,THEORY_VERSION

def cascade(a,b):
    den=1-a[:,1,1]*b[:,0,0];v=np.empty_like(a)
    v[:,0,0]=a[:,0,0]+a[:,0,1]*b[:,0,0]*a[:,1,0]/den
    v[:,1,1]=b[:,1,1]+b[:,1,0]*a[:,1,1]*b[:,0,1]/den
    v[:,1,0]=a[:,1,0]*b[:,1,0]/den;v[:,0,1]=b[:,0,1]*a[:,0,1]/den
    return v

def embedding(f,downstream,l=.74,c=0,cg=459,advance_ps=0):
    downstream=downstream.copy();w=2*np.pi*f
    downstream[:,0,0]*=np.exp(2j*w*advance_ps*1e-12)
    downstream[:,0,1]*=np.exp(1j*w*advance_ps*1e-12);downstream[:,1,0]*=np.exp(1j*w*advance_ps*1e-12)
    z=1/(1j*w*cg*1e-15);y=1j*w*c*1e-15+1/(1j*w*l*1e-9)
    series=np.empty_like(downstream);shunt=np.empty_like(downstream)
    series[:,0,0]=series[:,1,1]=z/(100+z);series[:,0,1]=series[:,1,0]=100/(100+z)
    shunt[:,0,0]=shunt[:,1,1]=-50*y/(2+50*y);shunt[:,0,1]=shunt[:,1,0]=2/(2+50*y)
    return cascade(cascade(series,shunt),downstream)

def main():
    net=read_snp(VERIFIED['artifacts']['touchstone']);mask=(net.f>=4.5e9)&(net.f<=5.5e9)
    f=net.f[mask];s=net.s[mask];target=theory(f)[:,1,0]
    reports=[]
    for advance in (0.,.22,4.39575829091819):
        def residual(x):
            pred=embedding(f,s,l=x[0],c=x[1],cg=x[2],advance_ps=advance)[:,1,0]
            d=pred-target*np.exp(-2j*np.pi*f*x[3]*1e-12)
            return np.r_[d.real,d.imag]
        seeds=[[l,c,cg,0] for l,c,cg in ((.74,0,459),(.74,30,459),(.74,80,459),(1.2,30,459),(.74,30,800))]
        fits=[least_squares(residual,x,bounds=([.2,0,100,-30],[2,200,1200,30]),
            x_scale=[.74,100,459,10],max_nfev=300) for x in seeds]
        q=min(fits,key=lambda a:a.fun@a.fun)
        reports.append({'downstream_port_advance_ps':advance,'best_L_nH_Cextra_fF_Cg_fF_delay_ps':q.x.tolist(),
            'complex_RMSE':float(np.sqrt(np.mean(q.fun**2)*2)),
            'interpretation':'Cextra is additional to the native bare input-coupler capacitance; advance other than 0 or ~0.22 ps is a diagnostic reference-plane hypothesis, not a physical deembedding proof'})
    result={'theory_version':THEORY_VERSION,'source':VERIFIED['artifacts']['touchstone'],'new_EM_frequencies':0,
        'checks':reports,'geometry_acceptance':False}
    (OUT/'ideal_embedding_control.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)

if __name__=='__main__':main()
