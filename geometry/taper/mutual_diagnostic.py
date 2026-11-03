"""Diagnostic coupled-inductor model; never silently replaces the target circuit.

All self-inductances stay at independently measured values. A fitted mutual
term tests one candidate explanation of the linewidth mismatch. A good fit
alone cannot uniquely identify a physical magnetic coupling mechanism.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import numpy as np
from scipy.optimize import least_squares
from circuit_common import read_snp
from geometry.taper.progressive_fixtures import ROOT,NATIVE,OUT
from geometry.taper.staged_fitting import LT1,LR
from geometry.taper.run_bridged_stage import digest


def bare(f,l,cg,mutual,delay,c0=0):
    w=2*np.pi*np.asarray(f)
    inductance=np.diag([l,LT1,LR]);inductance[0,1]=inductance[1,0]=mutual
    inv=np.linalg.inv(inductance*1e-9)
    cap=np.diag([c0,100.69929340643347,190.4])
    for i,j,g in ((0,1,cg),(1,2,8.2865)):
        cap[i,i]+=g;cap[j,j]+=g;cap[i,j]-=g;cap[j,i]-=g
    y=1j*w[:,None,None]*cap[None]*1e-15+inv[None]/(1j*w[:,None,None])
    ports=[0,2]
    yp=y[:,ports][:,:,ports]-y[:,ports,1,None]*y[:,None,1,ports]/y[:,1,1,None,None]
    den=(1+50*yp[:,0,0])*(1+50*yp[:,1,1])-2500*yp[:,0,1]*yp[:,1,0]
    return -100*yp[:,1,0]/den*np.exp(-1j*w*delay*1e-12)


def loaded(f,l,c0,gin,gj,mutual,delay,add=0):
    w=2*np.pi*np.asarray(f)
    inv=np.linalg.inv(np.array([[l,mutual],[mutual,LT1]])*1e-9)
    cap=np.diag([0.,c0+add,100.69929340643347])
    for i,j,g in ((0,1,gin),(1,2,gj)):
        cap[i,i]+=g;cap[j,j]+=g;cap[i,j]-=g;cap[j,i]-=g
    y=1j*w[:,None,None]*cap[None]*1e-15
    y[:,1:,1:]+=inv[None]/(1j*w[:,None,None])
    ports=[0,2]
    yp=y[:,ports][:,:,ports]-y[:,ports,1,None]*y[:,None,1,ports]/y[:,1,1,None,None]
    den=(1+50*yp[:,0,0])*(1+50*yp[:,1,1])-2500*yp[:,0,1]*yp[:,1,0]
    return -100*yp[:,1,0]/den*np.exp(-1j*w*delay*1e-12)


def full(f,l,c0,gin,gj,mutual,delay,n_regular=1):
    """Input series capacitor, coupled T2/T1 coils, then regular chain."""
    w=2*np.pi*np.asarray(f);ls=[l,LT1]+[LR]*n_regular
    inductance=np.diag(ls);inductance[0,1]=inductance[1,0]=mutual
    inv=np.linalg.inv(inductance*1e-9);size=len(ls)+1
    cap=np.diag([0,c0,100.69929340643347]+[190.4]*n_regular)
    for i,g in enumerate([gin,gj]+[8.2865]*n_regular):
        cap[i,i]+=g;cap[i+1,i+1]+=g;cap[i,i+1]-=g;cap[i+1,i]-=g
    y=1j*w[:,None,None]*cap[None]*1e-15
    y[:,1:,1:]+=inv[None]/(1j*w[:,None,None])
    ports=[0,size-1];inside=list(range(1,size-1))
    yp=y[:,ports][:,:,ports]-y[:,ports][:,:,inside]@np.linalg.solve(y[:,inside][:,:,inside],y[:,inside][:,:,ports])
    den=(1+50*yp[:,0,0])*(1+50*yp[:,1,1])-2500*yp[:,0,1]*yp[:,1,0]
    return -100*yp[:,1,0]/den*np.exp(-1j*w*delay*1e-12)


def synthesize(mutual=.1642629940183481):
    """Can positive physical C compensate M while retaining fixed T1 geometry?"""
    from geometry.taper.staged_fitting import circuit
    f=np.linspace(4.65e9,5.4e9,1001)
    target=circuit(f,[.74,6.06,4.9194],[100,100,190.4],[459,56.7,8.2865])[:,1,0]
    def residual(x):
        d=full(f,.74,*x[:3],mutual,x[3])-target;return np.r_[d.real,d.imag]
    bounds=([0,10,2,-100],[3000,5000,200,100])
    starts=[(c,g,j,0) for c,g,j in ((100,459,56.7),(100,800,80),(500,800,80),(1000,1200,100),(20,200,40),(500,200,80))]
    fits=[least_squares(residual,x,bounds=bounds,x_scale=[100,459,57,20],max_nfev=400,
                        ftol=1e-11,xtol=1e-11,gtol=1e-11) for x in starts]
    best=min(fits,key=lambda q:q.fun@q.fun);e=abs(full(f,.74,*best.x[:3],mutual,best.x[3])-target)
    wide=np.linspace(4e9,7e9,3001)
    target4=circuit(wide,[.74,6.06]+[4.9194]*4,[100,100]+[190.4]*4,[459,56.7]+[8.2865]*4)[:,1,0]
    actual4=full(wide,.74,*best.x[:3],mutual,best.x[3],n_regular=4)
    r=dict(M_nH=mutual,physical_C0_T2_fF=float(best.x[0]),physical_Cg_input_fF=float(best.x[1]),
           physical_Cg_T2_T1_fF=float(best.x[2]),delay_ps=float(best.x[3]),
           complex_RMS=float(np.sqrt(np.mean(e**2))),complex_max=float(e.max()),active_bounds=best.active_mask.tolist(),
           four_regular_wideband_RMS=float(np.sqrt(np.mean(abs(actual4-target4)**2))),
           diagnostic_only=True,scope='Numerical inverse circuit test, no geometry accepted')
    folder=OUT/'mutual_inverse_synthesis';folder.mkdir(parents=True,exist_ok=True)
    (folder/'results.json').write_text(json.dumps(r,indent=2)+'\n')
    print('MUTUAL SYNTHESIS',r,flush=True);return r


def fit_bare(name,l,*,source=None,output=None):
    source=source or NATIVE/name/'generated_ncell_3.snp'
    net=read_snp(source);ix=np.arange(0,len(net.f),5)
    def residual(x,indices=ix):
        d=bare(net.f[indices],l,*x)-net.s[indices,1,0];return np.r_[d.real,d.imag]
    fits=[least_squares(residual,[56.7,m,10],bounds=([10,-.8,-100],[150,.8,100]),
                        x_scale=[57,.1,20],max_nfev=300,ftol=1e-11,xtol=1e-11,gtol=1e-11)
          for m in (-.4,-.15,0,.15,.4)]
    best=min(fits,key=lambda q:q.fun@q.fun)
    best=least_squares(lambda x:residual(x,np.arange(len(net.f))),best.x,
                       bounds=([10,-.8,-100],[150,.8,100]),x_scale=[57,.1,20],max_nfev=200)
    model=bare(net.f,l,*best.x);error=abs(model-net.s[:,1,0])
    r=dict(source=str(source),source_sha256=digest(source),code_sha256=digest(__file__),
           L_T2_fixed_nH=l,Cg_T2_T1_fF=float(best.x[0]),mutual_nH=float(best.x[1]),delay_ps=float(best.x[2]),
           complex_RMS=float(np.sqrt(np.mean(error**2))),complex_max=float(error.max()),
           active_bounds=best.active_mask.tolist(),diagnostic_only=True,
           caveat='Extra coupled-inductor degree of freedom; success does not uniquely identify the physical mechanism.')
    out=output or OUT/name/'mutual_diagnostic';out.mkdir(parents=True,exist_ok=True)
    (out/'fit.json').write_text(json.dumps(r,indent=2)+'\n');np.savez(out/'response.npz',f=net.f,actual=net.s[:,1,0],model=model)
    print('MUTUAL BARE',name,r,flush=True);return r


def fit_loaded():
    from geometry.taper.measure_compact_inductance import AP
    from geometry.taper.response_cases import W
    folder=OUT/'03d_known_capacitor_operator';net=W['read_network'](folder/'operator.s4p',4)
    s,z0=AP['reduce_pairs'](net.s);y=AP['s_to_y'](s,z0)
    basis=lambda f:(np.asarray(f)[:,None]/5e9)**np.array([-1,1,3,5])
    held=np.isclose(net.f,5.5e9,rtol=0,atol=1)
    co=np.linalg.solve(basis(net.f[~held]),y[~held].reshape(4,-1))
    f=np.linspace(4e9,7e9,1001);yy=(basis(f)@co).reshape(len(f),3,3);adds=[0,400,800,1200]
    actual=np.array([AP['add_capacitors'](yy,f,[c])[:,1,0] for c in adds])
    def prediction(x):return np.array([loaded(f,.7341365873386546,*x,add=c) for c in adds])
    def residual(x):
        d=(prediction(x)-actual).ravel();return np.r_[d.real,d.imag]
    fits=[least_squares(residual,[100,459,56.7,m,10],bounds=([0,30,5,-.8,-100],[1000,3000,180,.8,100]),
                        x_scale=[100,459,57,.1,20],max_nfev=400,ftol=1e-10,xtol=1e-10,gtol=1e-10)
          for m in (-.2,0,.2)]
    best=min(fits,key=lambda q:q.fun@q.fun);error=abs(prediction(best.x)-actual)
    r=dict(parameter_names=['C0_T2_fF','Cg_input_fF','Cg_T2_T1_fF','mutual_nH','delay_ps'],
           values=best.x.tolist(),complex_RMS=float(np.sqrt(np.mean(error**2))),complex_max=float(error.max()),
           active_bounds=best.active_mask.tolist(),known_Cadd_fF=adds,diagnostic_only=True)
    (folder/'mutual_diagnostic.json').write_text(json.dumps(r,indent=2)+'\n')
    print('MUTUAL LOADED',r,flush=True);return r


if __name__=='__main__':
    fit_bare('07_plain_T2_direct_regular',.7327695472953403)
    fit_bare('06_loop_return_T2_direct_regular',.7447382669006548)
    fit_loaded()
