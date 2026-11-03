"""Two-electrical-parameter check of restored IDC fixtures, without running EM.

T1 and regular values remain fixed. A separately labelled total fixture delay
is a nuisance parameter. Convergence alone never establishes design acceptance.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks,peak_widths
from circuit_common import read_snp
from geometry.taper.staged_fitting import circuit,LT1,CT1,LR,CR,GOUT
from geometry.taper.run_bridged_stage import digest
from geometry.taper.optimize_response import theory
import matplotlib.pyplot as plt


def features(f,s):
    v=abs(s[:,1,0]);peaks,_=find_peaks(v,prominence=.025)
    width=peak_widths(v*v,peaks)[0]*(f[1]-f[0])/1e6
    return [dict(GHz=float(f[i]/1e9),magnitude=float(v[i]),power_FWHM_MHz=float(w))
            for i,w in zip(peaks,width)]


def fit(snp,out,l2,*,downstream_g=57.09412572268222,data_label='Sonnet'):
    snp=Path(snp);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    net=read_snp(snp);f=net.f;observed=net.s
    ix=np.unique(np.linspace(0,len(f)-1,min(1501,len(f))).astype(int))
    def predict(x,indices):
        return circuit(f[indices],[l2,LT1,LR,LR],[x[0],CT1,CR,CR],
            [x[1],downstream_g,8.40642914296533,GOUT],left_ps=x[2])
    def residual(x,indices=ix,denominator=False):
        a=predict(x,indices)[:,1,0];b=observed[indices,1,0]
        d=(1/a-1/b)*np.maximum(abs(b),.03)**2 if denominator else a-b
        return np.r_[d.real,d.imag]
    bounds=([0,30,-30],[2000,3000,30])
    starts=[(c,g,5) for c in (50,100,500) for g in (150,459,1200)]
    smooth=[least_squares(lambda x:residual(x,denominator=True),x,bounds=bounds,
                         x_scale=[100,459,5],max_nfev=150) for x in starts]
    seeds=sorted(smooth,key=lambda q:np.linalg.norm(residual(q.x)))[:3]
    candidates=[least_squares(residual,x,bounds=bounds,x_scale=[100,459,5],max_nfev=200)
                for x in [q.x for q in seeds]+[[100,459,5]]]
    best=min(candidates,key=lambda q:q.fun@q.fun);full=np.arange(len(f))
    best=least_squares(lambda x:residual(x,full),best.x,bounds=bounds,
                       x_scale=[100,459,5],max_nfev=150)
    model=predict(best.x,full);target=theory(f)
    report=dict(source=str(snp.resolve()),source_sha256=digest(snp),data_label=data_label,
        free_electrical_parameters=['C0_T2_fF','Cg_input_T2_fF'],
        L_T2_nH=l2,C0_T2_fF=float(best.x[0]),Cg_input_T2_fF=float(best.x[1]),
        fixture_delay_ps=float(best.x[2]),active_bounds=best.active_mask.tolist(),
        errors_percent=dict(L_T2=100*(l2/.74-1),C0_T2=100*(best.x[0]/100-1),
                            Cg_input_T2=100*(best.x[1]/459-1)),
        fixed=dict(L_T1_nH=LT1,C0_T1_fF=CT1,Cg_T2_T1_fF=downstream_g,
            L_regular_nH=LR,C0_regular_fF=CR,Cg_T1_regular_fF=8.40642914296533,Cg_output_fF=GOUT),
        S21_complex_RMSE=float(np.sqrt(np.mean(abs(model[:,1,0]-observed[:,1,0])**2))),
        S21_max_complex_error=float(np.max(abs(model[:,1,0]-observed[:,1,0]))),
        S21_magnitude_RMSE=float(np.sqrt(np.mean((abs(model[:,1,0])-abs(observed[:,1,0]))**2))),
        nominal_S21_magnitude_RMSE=float(np.sqrt(np.mean((abs(target[:,1,0])-abs(observed[:,1,0]))**2))),
        peaks=features(f,observed),fitted_peaks=features(f,model),target_peaks=features(f,target),
        accepted=False,
        limitation='Downstream Cg comes from the wire-removed control with provisional C0. Its transfer must be validated. No automatic final acceptance.')
    (out/'fit.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez(out/'response.npz',f=f,observed=observed,model=model,target=target)
    fig,ax=plt.subplots(figsize=(12,7))
    ax.plot(f/1e9,abs(observed[:,1,0]),color='#b22222',lw=2.8,label=data_label)
    ax.plot(f/1e9,abs(model[:,1,0]),'k--',lw=2.2,label='Complex circuit fit')
    ax.plot(f/1e9,abs(target[:,1,0]),color='#246cb5',ls=':',lw=2.2,label='Target circuit')
    text=(f"Exploratory — not accepted\nLₜ₂ = {l2:.5f} nH ({report['errors_percent']['L_T2']:+.2f}%)\n"
          f"C₀,ₜ₂ = {best.x[0]:.2f} fF ({report['errors_percent']['C0_T2']:+.2f}%)\n"
          f"Cg,ₜ₂ = {best.x[1]:.2f} fF ({report['errors_percent']['Cg_input_T2']:+.2f}%)\n"
          f"Cg,ₜ₁ = {downstream_g:.3f} fF (fixed control)\nComplex S₂₁ RMS = {report['S21_complex_RMSE']:.4f}")
    ax.text(.035,.96,text,transform=ax.transAxes,va='top',fontsize=13,
            bbox=dict(boxstyle='round,pad=.5',fc='white',ec='gray'))
    ax.set_ylim(0,max(abs(observed[:,1,0]).max(),abs(model[:,1,0]).max(),.7)*1.65)
    ax.set_title('Taper 2 Transmission — Design Validation',fontsize=20,pad=14)
    ax.set_xlabel('Frequency (GHz)',fontsize=17);ax.set_ylabel(r'$|S_{21}|$',fontsize=20)
    ax.tick_params(labelsize=15);ax.legend(fontsize=13,loc='upper right');ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(out/'s21_fit_linear.png',dpi=220);plt.close(fig)
    print(json.dumps(report,indent=2),flush=True)
    return report
