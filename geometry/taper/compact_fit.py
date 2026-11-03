"""Fit the no-IDC stage: bare T2 parasitic C and T2-to-T1 coupling only.

Bare parasitic C is not the eventual 100 fF IDC-equipped target. T1 and the
regular circuit remain fixed. Fixture delay is labelled separately.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from circuit_common import read_snp
from geometry.taper.staged_fitting import circuit,LT1,CT1,LR,CR,GOUT
from geometry.taper.paired_coupling import abcd
from geometry.taper.fit_restored_t2 import features
from geometry.taper.run_bridged_stage import digest
import matplotlib.pyplot as plt

GJ=8.40642914296533


def direct_predictor(f,l2):
    f=np.asarray(f);w=2*np.pi*f
    n=abcd(circuit(f,[LT1,LR,LR],[CT1,CR,CR],[None,GJ,GOUT]))
    def predict(x):
        c,g,tau=x;z=1/(1j*w*g*1e-15)
        a=n[:,0,0]+z*n[:,1,0];b=n[:,0,1]+z*n[:,1,1]
        y=1j*w*c*1e-15+1/(1j*w*l2*1e-9)
        denominator=(a+b/50)*(1+50*y)+50*n[:,1,0]+n[:,1,1]
        return 2/denominator*np.exp(-1j*w*tau*1e-12)
    return predict


def fit_direct(snp,out,l2,*,known_cg=None,junction_idc=False,data_label='Sonnet, no IDC'):
    snp=Path(snp);out=Path(out);out.mkdir(parents=True,exist_ok=True)
    n=read_snp(snp);f=n.f;actual=n.s
    ix=np.unique(np.linspace(0,len(f)-1,min(1501,len(f))).astype(int))
    predict=direct_predictor(f[ix],l2);obs=actual[ix,1,0]
    expand=lambda x: [x[0],known_cg,x[1]] if known_cg is not None else x
    def residual(x):
        d=predict(expand(x))-obs;return np.r_[d.real,d.imag]
    bounds=([0,-20],[2000,30]) if known_cg is not None else ([0,5,-20],[2000,180,30])
    scale=[100,5] if known_cg is not None else [100,57,5]
    seeds=[]
    for c in (0,30,60,120,300,800):
        # Narrow peaks move by many linewidths on a coarse Cg grid. Resolve
        # the target neighbourhood before selecting least-squares starts.
        for g in ([known_cg] if known_cg is not None else np.unique(np.r_[np.linspace(5,180,36),np.arange(15,150.001,.125)])):
            v=predict([c,g,0])
            for tau in (0,5,15):
                cost=np.mean(abs(v*np.exp(-2j*np.pi*f[ix]*tau*1e-12)-obs)**2)
                seeds.append((cost,[c,tau] if known_cg is not None else [c,g,tau]))
    seeds.sort(key=lambda q:q[0])
    fits=[least_squares(residual,x,bounds=bounds,x_scale=scale,max_nfev=200,
                        ftol=1e-11,xtol=1e-11,gtol=1e-11) for _,x in seeds[:8]]
    best=min(fits,key=lambda q:q.fun@q.fun)
    full_predict=direct_predictor(f,l2)
    def full_res(x):
        d=full_predict(expand(x))-actual[:,1,0];return np.r_[d.real,d.imag]
    best=least_squares(full_res,best.x,bounds=bounds,x_scale=scale,max_nfev=100)
    c,g,tau=expand(best.x)
    model=circuit(f,[l2,LT1,LR,LR],[c,CT1,CR,CR],[None,g,GJ,GOUT],left_ps=tau)
    error=abs(model[:,1,0]-actual[:,1,0])
    report=dict(source=str(snp.resolve()),source_sha256=digest(snp),data_label=data_label,IDC_present=junction_idc,
        input_IDC_present=False,junction_IDC_present=junction_idc,
        free_electrical_parameters=['Cparasitic_bare_T2_fF'] if known_cg is not None else ['Cparasitic_bare_T2_fF','Cg_T2_T1_fF'],
        known_ideal_junction_Cg_fF=known_cg,
        L_T2_nH=float(l2),Cparasitic_bare_T2_fF=float(c),Cg_T2_T1_fF=float(g),fixture_delay_ps=float(tau),
        L_error_percent=float(100*(l2/.74-1)),Cg_error_percent=float(100*(g/56.7-1)),
        active_bounds=best.active_mask.tolist(),S21_complex_RMSE=float(np.sqrt(np.mean(error**2))),
        S21_max_complex_error=float(error.max()),
        S21_magnitude_RMSE=float(np.sqrt(np.mean((abs(model[:,1,0])-abs(actual[:,1,0]))**2))),
        peaks=features(f,actual),model_peaks=features(f,model),
        fixed=dict(L_T1_nH=LT1,C0_T1_fF=CT1,L_regular_nH=LR,C0_regular_fF=CR,
                   Cg_T1_regular_fF=GJ,Cg_output_fF=GOUT),
        complete_T2_accepted=False,
        limitation='No input IDC. Bare parasitic capacitance is not constrained to the eventual 100 fF target. Independent coupling control and assembly checks remain required.')
    (out/'fit.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez(out/'response.npz',f=f,actual=actual,model=model)
    fig,ax=plt.subplots(figsize=(12,7))
    ax.plot(f/1e9,abs(actual[:,1,0]),color='#b22222',lw=2.6,label=data_label)
    ax.plot(f/1e9,abs(model[:,1,0]),'k--',lw=2.2,label='Fixed downstream circuit + T2 fit')
    box=(("Known ideal junction control\n" if known_cg is not None else "No-input-IDC stage\n")+f"L T2 = {l2:.5f} nH ({report['L_error_percent']:+.2f}%)\n"
         f"Cg T2 to T1 = {g:.3f} fF ({report['Cg_error_percent']:+.2f}%)\n"
         f"Bare parasitic C = {c:.2f} fF\nComplex S21 RMS = {report['S21_complex_RMSE']:.4f}")
    ax.text(.03,.96,box,transform=ax.transAxes,va='top',fontsize=13,
            bbox=dict(boxstyle='round,pad=.5',fc='white',ec='gray'))
    ax.set_ylim(0,max(abs(actual[:,1,0]).max(),abs(model[:,1,0]).max())*1.6)
    ax.set_title('Compact Taper 2: Direct-Feed Validation',fontsize=20,pad=14)
    ax.set_xlabel('Frequency (GHz)',fontsize=17);ax.set_ylabel(r'$|S_{21}|$',fontsize=20)
    ax.tick_params(labelsize=15);ax.legend(fontsize=12,loc='upper right');ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(out/'s21_fit_linear.png',dpi=210);plt.close(fig)
    print(json.dumps(report,indent=2),flush=True)
    return report
