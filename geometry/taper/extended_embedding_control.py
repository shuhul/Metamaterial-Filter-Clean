"""Test physical extra matching reactances against the preserved T1 EM block."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import numpy as np
from scipy.optimize import least_squares
from geometry.taper.ideal_embedding_control import cascade,read_snp,VERIFIED,OUT,theory

def series(z):
    s=np.empty((len(z),2,2),complex)
    s[:,0,0]=s[:,1,1]=z/(100+z);s[:,0,1]=s[:,1,0]=100/(100+z)
    return s

def shunt(y):
    s=np.empty((len(y),2,2),complex)
    s[:,0,0]=s[:,1,1]=-50*y/(2+50*y);s[:,0,1]=s[:,1,0]=2/(2+50*y)
    return s

def main():
    net=read_snp(VERIFIED['artifacts']['touchstone']);mask=(net.f>=4.5e9)&(net.f<=5.5e9)
    f=net.f[mask];w=2*np.pi*f;s=net.s[mask];target=theory(f)
    # Lnode, additional Cnode, Cgin, input series L, input shunt C,
    # right series L, total reference-plane delay. All physical reactances >=0.
    def response(x):
        a=cascade(shunt(1j*w*x[4]*1e-15),series(1/(1j*w*x[2]*1e-15)+1j*w*x[3]*1e-9))
        a=cascade(a,shunt(1j*w*x[1]*1e-15+1/(1j*w*x[0]*1e-9)))
        return cascade(cascade(a,series(1j*w*x[5]*1e-9)),s)
    def residual(x):
        d=response(x)[:,1,0]-target[:,1,0]*np.exp(-1j*w*x[6]*1e-12)
        return np.r_[d.real,d.imag]
    lower=[.2,0,100,0,0,0,-30];upper=[3,1000,1600,3,500,3,30]
    starts=[[.92,1,341,ls,ci,lr,-8] for ls,ci,lr in ((0,0,0),(.2,70,0),(.4,70,.2),(.1,150,.5),(1,50,.1),(.3,0,1))]
    fits=[least_squares(residual,x,bounds=(lower,upper),x_scale=[.74,100,459,.3,100,.2,10],max_nfev=600) for x in starts]
    q=min(fits,key=lambda a:a.fun@a.fun)
    r={'parameters':['Lnode_nH','additional_Cnode_fF','Cgin_fF','input_series_L_nH','input_shunt_C_fF','right_series_L_nH','delay_ps'],
        'best_x':q.x.tolist(),'complex_S21_RMSE':float(np.sqrt(np.mean(q.fun*q.fun)*2)),
        'active_bounds':q.active_mask.tolist(),'new_EM_frequencies':0,
        'all_starts':[{'x':v.x.tolist(),'cost':float(v.fun@v.fun)} for v in fits],
        'physical_acceptance':False}
    (OUT/'extended_embedding_control.json').write_text(json.dumps(r,indent=2)+'\n')
    print(json.dumps(r,indent=2),flush=True)

if __name__=='__main__':main()
