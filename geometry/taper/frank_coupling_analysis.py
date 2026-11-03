"""Compare coupling-only, bridged T1 and full T1 with coherent fixed references.

No EM calls. Every fitted electrical parameter and fixture nuisance parameter
is recorded. Fit agreement is not a guarantee of unique physical capacitance.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import numpy as np
from scipy.optimize import least_squares,minimize_scalar
from circuit_common import read_snp
from geometry.taper.frank_coupling_controls import ROOT,OUT,NATIVE,CASES
from geometry.taper.staged_fitting import circuit,LT1,CT1,LR,CR,GOUT
from geometry.taper.locked_taper2 import FIXED,VERIFIED
from geometry.taper.run_bridged_stage import digest
from geometry.taper.minimum_inductance import half_power_features

PROFILES=dict(N4_independent=dict(cr=CR,go=GOUT),
              prior_T1=dict(cr=FIXED['C_regular_fF'],go=FIXED['Cg_output_fF']))


def error_report(actual,model):
    error=abs(actual-model)
    return dict(S21_complex_RMS=float(np.sqrt(np.mean(error[:,1,0]**2))),
                S21_complex_max=float(error[:,1,0].max()),
                S21_magnitude_RMS=float(np.sqrt(np.mean((abs(actual[:,1,0])-abs(model[:,1,0]))**2))),
                unfitted_S11_RMS=float(np.sqrt(np.mean(error[:,0,0]**2))),
                unfitted_S22_RMS=float(np.sqrt(np.mean(error[:,1,1]**2))))


def plot(f,actual,model,r,path,title):
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(12,7))
    ax.plot(f/1e9,abs(actual[:,1,0]),color='#b22222',lw=2.8,label='Sonnet simulation')
    ax.plot(f/1e9,abs(model[:,1,0]),'k--',lw=2.3,label='Complex S21 circuit fit')
    box=(f"Cg (T1 to regular) = {r['Cg_junction_fF']:.5f} fF\n"
         f"Target difference = {100*(r['Cg_junction_fF']/8.2865-1):+.2f}%\n"
         f"Complex S21 RMS = {r['S21_complex_RMS']:.5f}\n"
         f"Fixed regular reference: {r['reference_profile']}\nDiagnostic extraction")
    ax.text(.03,.96,box,transform=ax.transAxes,va='top',fontsize=14,
            bbox=dict(boxstyle='round,pad=.5',fc='white',ec='gray'))
    ax.set_ylim(0,max(abs(actual[:,1,0]).max(),abs(model[:,1,0]).max())*1.6)
    ax.set_xlabel('Frequency (GHz)',fontsize=18);ax.set_ylabel(r'$|S_{21}|$',fontsize=22)
    ax.set_title(title,fontsize=21,pad=14);ax.tick_params(labelsize=15)
    ax.legend(loc='upper right',fontsize=13);ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(path,dpi=190);plt.close(fig)


def fit(snp,out,kind,profile,*,source_cap=False,junction_fixed=None,regular_C_override=None,output_Cg_override=None):
    """Cg only for direct controls; full T1 uses C1/Gin and optionally Cg."""
    out.mkdir(parents=True,exist_ok=True)
    n=read_snp(snp);mask=(n.f>=4.75e9-1)&(n.f<=5.3e9+1)
    f=n.f[mask];actual=n.s[mask]
    ix=np.unique(np.linspace(0,len(f)-1,min(1501,len(f))).astype(int))
    cr,go=PROFILES[profile]['cr'],PROFILES[profile]['go']
    if output_Cg_override is not None:go=output_Cg_override
    main_cr=cr if regular_C_override is None else regular_C_override
    if kind=='full_T1':
        names=['C0_T1_fF','Cg_input_T1_fF']+([] if junction_fixed is not None else ['Cg_junction_fF'])+['delay_ps']
        lower=[80,35]+([] if junction_fixed is not None else [5])+[0]
        upper=[120,85]+([] if junction_fixed is not None else [12])+[25]
        seeds=[[c,g]+([] if junction_fixed is not None else [j])+[5]
               for c in (98,102) for g in (50,60) for j in ((junction_fixed,) if junction_fixed is not None else (8,8.6))]
    else:
        names=['Cg_junction_fF','delay_ps']+(['source_C_fF'] if source_cap else [])
        lower=[3,0]+([0] if source_cap else []);upper=[15,30]+([1000] if source_cap else [])
        seeds=[[g,7]+([100] if source_cap else []) for g in (6,8.2865,11)]
    def model(x,frequency=f,split=.8):
        if kind=='full_T1':
            ct,gin=x[:2];gj=junction_fixed if junction_fixed is not None else x[2];tau=x[-1]
            ls=[LT1,LR,LR];cs=[ct,main_cr,cr];gs=[gin,gj,go]
        else:
            gj,tau=x[:2]
            if kind=='bridged_T1':ls=[LT1,LR,LR];cs=[CT1,main_cr,cr];gs=[None,gj,go]
            elif source_cap:ls=[1e12,LR,LR];cs=[x[2],main_cr,cr];gs=[None,gj,go]
            else:ls=[LR,LR];cs=[main_cr,cr];gs=[gj,go]
        return circuit(frequency,ls,cs,gs,left_ps=tau*split,right_ps=tau*(1-split))
    def residual(x,indices=ix):
        d=model(x,f[indices])[:,1,0]-actual[indices,1,0]
        return np.r_[d.real,d.imag]
    candidates=[least_squares(residual,x,bounds=(lower,upper),x_scale='jac',max_nfev=300,
                             ftol=1e-11,xtol=1e-11,gtol=1e-11) for x in seeds]
    best=min(candidates,key=lambda q:q.fun@q.fun)
    best=least_squares(lambda x:residual(x,np.arange(len(f))),best.x,bounds=(lower,upper),
                       x_scale='jac',max_nfev=200,ftol=1e-11,xtol=1e-11,gtol=1e-11)
    # S21 determines only total matched delay. Reflection diagnostics use the
    # best split at that fixed total; no electrical parameter is changed here.
    def reflection_cost(split):
        s=model(best.x,split=split)
        return np.mean(abs(s[:,0,0]-actual[:,0,0])**2+abs(s[:,1,1]-actual[:,1,1])**2)
    split=minimize_scalar(reflection_cost,bounds=(0,1),method='bounded').x
    prediction=model(best.x,split=split)
    params=dict(zip(names,map(float,best.x)))
    if junction_fixed is not None:params['Cg_junction_fF']=float(junction_fixed)
    r=dict(source=str(snp.resolve()),source_sha256=digest(snp),kind=kind,reference_profile=profile,
           **params,**error_report(actual,prediction),active_bounds=best.active_mask.tolist(),
           free_electrical_parameters=[q for q in names if q!='delay_ps'],
           fixed_junction_fF=junction_fixed,fixture_delay_split_left=float(split),
           fixed_regular=dict(L_nH=LR,C0_fF=main_cr,C0_end_fF=cr,Cg_output_fF=go,output_end_LC_included=True),
           samples=len(f),interval_GHz=[float(f[0]/1e9),float(f[-1]/1e9)],
           observed_peaks=half_power_features(f,actual[:,1,0]),
           model_peaks=half_power_features(f,prediction[:,1,0]),
           extraction_scope='Effective circuit coupling, fixed 2 um mesh; no absolute/unique physical-C guarantee')
    (out/'fit.json').write_text(json.dumps(r,indent=2)+'\n')
    np.savez(out/'response.npz',f=f,actual=actual,model=prediction)
    plot(f,actual,prediction,r,out/'s21_fit.png',kind.replace('_',' ').title()+': Coupling Check')
    print(json.dumps(r),flush=True)
    return r


def run(*,include_source_C_diagnostic=False):
    rows=[]
    for case in CASES:
        snp=NATIVE/case/'generated_ncell_2.snp'
        if not snp.exists():continue
        status=OUT/case/'run_status.json'
        if not status.exists() or json.loads(status.read_text()).get('status')!='complete':continue
        for profile in PROFILES:
            rows.append(dict(case=case,fit=fit(snp,OUT/case/profile,'coupling_only',profile)))
            if include_source_C_diagnostic:
                rows.append(dict(case=case+'_sourceC_diagnostic',fit=fit(snp,OUT/case/(profile+'_sourceC'),'coupling_only',profile,source_cap=True)))
    direct=ROOT/'sonnet_files/taper_exploration/bridged_feed_tests/01_t1_input_bridged/generated_ncell_2.snp'
    full=ROOT/VERIFIED['artifacts']['touchstone']
    for profile in PROFILES:
        rows.append(dict(case='bridged_T1_existing',fit=fit(direct,OUT/'bridged_T1_existing'/profile,'bridged_T1',profile)))
        rows.append(dict(case='full_T1_free_junction',fit=fit(full,OUT/'full_T1_free_junction'/profile,'full_T1',profile)))
        primary=next((r['fit'] for r in rows if r['case']=='minimal_feed10' and r['fit']['reference_profile']==profile),None)
        if primary is not None:
            rows.append(dict(case='full_T1_control_junction_fixed',fit=fit(full,OUT/'full_T1_control_junction_fixed'/profile,
                                        'full_T1',profile,junction_fixed=primary['Cg_junction_fF'])))
    (OUT/'comparison.json').write_text(json.dumps(rows,indent=2)+'\n')
    return rows


def regular_sensitivity(case='minimal_feed10'):
    """Diagnostic only: can a small regular-C error explain the one-C fit?

    Keep L and output coupling/end-cell C fixed. Release only the resonating
    regular cell C within +/-2%; never update a production reference file.
    """
    snp=NATIVE/case/'generated_ncell_2.snp';n=read_snp(snp)
    mask=(n.f>=4.75e9-1)&(n.f<=5.3e9+1);f=n.f[mask];actual=n.s[mask]
    ix=np.unique(np.linspace(0,len(f)-1,1501).astype(int));rows=[]
    for profile in PROFILES:
        cr,go=PROFILES[profile]['cr'],PROFILES[profile]['go']
        def predict(x,frequency=f,split=.8):
            return circuit(frequency,[LR,LR],[x[1],cr],[x[0],go],left_ps=x[2]*split,right_ps=x[2]*(1-split))
        def residual(x,indices=ix):
            d=predict(x,f[indices])[:,1,0]-actual[indices,1,0];return np.r_[d.real,d.imag]
        bounds=([3,cr*.98,0],[15,cr*1.02,30])
        fits=[least_squares(residual,[g,c,6],bounds=bounds,x_scale=[8,190,6],max_nfev=300,
                              ftol=1e-11,xtol=1e-11,gtol=1e-11)
              for g in (6,9) for c in (cr-1,cr+1)]
        q=min(fits,key=lambda p:p.fun@p.fun)
        q=least_squares(lambda x:residual(x,np.arange(len(f))),q.x,bounds=bounds,x_scale=[8,190,6],max_nfev=200)
        split=minimize_scalar(lambda s:np.mean(abs(predict(q.x,split=s)[:,0,0]-actual[:,0,0])**2+
                                             abs(predict(q.x,split=s)[:,1,1]-actual[:,1,1])**2),
                              bounds=(0,1),method='bounded').x
        model=predict(q.x,split=split)
        r=dict(case=case,reference_profile=profile,diagnostic_only=True,production_parameters_changed=False,
               Cg_junction_fF=float(q.x[0]),C0_regular_fF=float(q.x[1]),regular_C_change_percent=float(100*(q.x[1]/cr-1)),
               delay_ps=float(q.x[2]),active_bounds=q.active_mask.tolist(),**error_report(actual,model),
               original_regular_C_fF=cr,fixed_end_regular_C_fF=cr,fixed_output_Cg_fF=go,
               observed_peaks=half_power_features(f,actual[:,1,0]),model_peaks=half_power_features(f,model[:,1,0]))
        out=OUT/case/(profile+'_regular_C_sensitivity');out.mkdir(exist_ok=True)
        (out/'fit.json').write_text(json.dumps(r,indent=2)+'\n');np.savez(out/'response.npz',f=f,actual=actual,model=model)
        plot(f,actual,model,r,out/'s21_fit.png','Coupling Only: Regular-C Sensitivity Diagnostic')
        rows.append(r);print(json.dumps(r),flush=True)
    return rows


def joint_port_fit(case='minimal_feed10'):
    """Independent Cin/Cout/C0 extraction using S11, S21 and S22 together.

    Fix both inductances and output-end C. Electrical changes are local
    diagnostics; the production geometry and references are never modified.
    """
    if case=='bridged_T1_existing':
        snp=ROOT/'sonnet_files/taper_exploration/bridged_feed_tests/01_t1_input_bridged/generated_ncell_2.snp'
    else:snp=NATIVE/case/'generated_ncell_2.snp'
    n=read_snp(snp)
    mask=(n.f>=4.75e9-1)&(n.f<=5.3e9+1);f=n.f[mask];actual=n.s[mask]
    ix=np.unique(np.linspace(0,len(f)-1,1501).astype(int))
    def predict(x,frequency=f):
        if case=='bridged_T1_existing':
            return circuit(frequency,[LT1,LR,LR],[CT1,x[1],CR],[None,x[0],x[2]],left_ps=x[3],right_ps=x[4])
        return circuit(frequency,[LR,LR],[x[1],CR],[x[0],x[2]],left_ps=x[3],right_ps=x[4])
    def residual(x,indices=ix):
        q=predict(x,f[indices])-actual[indices]
        d=q[:,[0,1,1],[0,0,1]].ravel()
        return np.r_[d.real,d.imag]
    bounds=([3,185,6,0,0],[15,195,10,15,8])
    starts=[[g,c,8.25,4.8,1.2] for g in (5.8,8.3) for c in (189,191)]
    fits=[least_squares(residual,x,bounds=bounds,x_scale=[8,190,8,5,1],max_nfev=300,
                        ftol=1e-11,xtol=1e-11,gtol=1e-11) for x in starts]
    q=min(fits,key=lambda p:p.fun@p.fun)
    q=least_squares(lambda x:residual(x,np.arange(len(f))),q.x,bounds=bounds,x_scale=[8,190,8,5,1],max_nfev=200)
    model=predict(q.x)
    jac=q.jac*np.array([8,190,8,5,1])[None,:]
    condition=float(np.linalg.cond(jac))
    r=dict(case=case,source=str(snp),source_sha256=digest(snp),reference_profile='joint_S11_S21_S22',diagnostic_only=True,production_parameters_changed=False,
           objective='Equal-weight complex S11, S21, S22; S12 held out',
           Cg_junction_fF=float(q.x[0]),C0_regular_fF=float(q.x[1]),Cg_output_fF=float(q.x[2]),
           delay_left_ps=float(q.x[3]),delay_right_ps=float(q.x[4]),active_bounds=q.active_mask.tolist(),
           multistart_solutions=[dict(x=p.x.tolist(),success=bool(p.success),cost=float(p.fun@p.fun)) for p in fits],
           multistart_parameter_spread=np.ptp(np.array([p.x for p in fits]),axis=0).tolist(),
           scaled_Jacobian_condition=condition,**error_report(actual,model),
           heldout_S12_RMS=float(np.sqrt(np.mean(abs(actual[:,0,1]-model[:,0,1])**2))),
           observed_peaks=half_power_features(f,actual[:,1,0]),model_peaks=half_power_features(f,model[:,1,0]))
    r['S11_complex_RMS']=r.pop('unfitted_S11_RMS');r['S22_complex_RMS']=r.pop('unfitted_S22_RMS')
    r['S12_scope']='Reciprocity consistency check, not an independent physical validation'
    out=OUT/case/'joint_port_fit';out.mkdir(parents=True,exist_ok=True)
    (out/'fit.json').write_text(json.dumps(r,indent=2)+'\n');np.savez(out/'response.npz',f=f,actual=actual,model=model)
    plot(f,actual,model,r,out/'s21_fit.png','Coupling Only: Joint Transmission and Reflection Fit')
    print(json.dumps(r),flush=True)
    return r


def transfer_joint_control():
    """Freeze the complete isolated-reference tuple, then restore T1."""
    measured=json.loads((OUT/'minimal_feed10/joint_port_fit/fit.json').read_text())
    kwargs=dict(regular_C_override=measured['C0_regular_fF'],output_Cg_override=measured['Cg_output_fF'])
    full=ROOT/VERIFIED['artifacts']['touchstone']
    direct=ROOT/'sonnet_files/taper_exploration/bridged_feed_tests/01_t1_input_bridged/generated_ncell_2.snp'
    rows=[]
    for kind,snp,fixed,name in [('full_T1',full,measured['Cg_junction_fF'],'full_T1_fixed_control'),
                              ('full_T1',full,None,'full_T1_free_junction'),
                              ('bridged_T1',direct,None,'bridged_T1')]:
        r=fit(snp,OUT/'joint_reference_transfer'/name,kind,'N4_independent',junction_fixed=fixed,**kwargs)
        rows.append(dict(case=name,fit=r))
    (OUT/'joint_reference_transfer/comparison.json').write_text(json.dumps(rows,indent=2)+'\n')
    return rows


if __name__=='__main__':run()
