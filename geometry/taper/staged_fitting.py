"""One-electrical-parameter fits for the reviewed feed-driven fixtures.

All retained resonators, including directly fed ones, remain in the model.
Reference-plane delay is a separately labelled nuisance parameter, not a C/L.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json
import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from circuit_common import read_snp
from geometry.taper.locked_taper2 import FIXED

ROOT=Path(__file__).resolve().parents[2]
UNIT_RESULT=ROOT/'output_data/generated_ncell_4_w186_d56/fit_parameters.json'
UNIT=json.loads(UNIT_RESULT.read_text())

LR=FIXED['L_regular_nH']
CR=UNIT['C0_fF']
LT1=FIXED['L_taper_nH']
CT1=FIXED['C_taper_fF']
GOUT=UNIT['Cg_fF']


def circuit(f,ls,cs,gs,left_ps=0,right_ps=0,line_z=50):
    """Series coupling then shunt LC for every branch; None means direct feed.

The last branch is the feed-connected output end. No capacitance is inserted
after it. Equal port reference impedances are 50 ohms.
"""
    if not len(ls)==len(cs)==len(gs):raise ValueError('Every branch needs an input coupling or None')
    w=2*np.pi*np.asarray(f)
    m=np.tile(np.eye(2,dtype=complex),(len(w),1,1))
    def line(tau):
        a=np.zeros_like(m);theta=w*tau*1e-12
        a[:,0,0]=a[:,1,1]=np.cos(theta)
        a[:,0,1]=1j*line_z*np.sin(theta);a[:,1,0]=1j*np.sin(theta)/line_z
        return a
    m=m@line(left_ps)
    for l,c,g in zip(ls,cs,gs):
        if g is not None:
            a=np.tile(np.eye(2,dtype=complex),(len(w),1,1));a[:,0,1]=1/(1j*w*g*1e-15);m=m@a
        a=np.tile(np.eye(2,dtype=complex),(len(w),1,1))
        a[:,1,0]=1j*w*c*1e-15+1/(1j*w*l*1e-9);m=m@a
    m=m@line(right_ps)
    a,b,c,d=m[:,0,0],m[:,0,1],m[:,1,0],m[:,1,1]
    denominator=a+b/50+c*50+d
    s=np.empty_like(m)
    s[:,0,0]=(a+b/50-c*50-d)/denominator
    s[:,1,1]=(-a+b/50-c*50+d)/denominator
    s[:,0,1]=s[:,1,0]=2/denominator
    return s


def save_plot(f,actual,model,report,out,title):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    fig,ax=plt.subplots(figsize=(11,6))
    ax.plot(f/1e9,abs(actual[:,1,0]),color='#b22222',lw=2.5,label='Sonnet')
    ax.plot(f/1e9,abs(model[:,1,0]),'k--',lw=2.2,label='One-parameter circuit fit')
    key=report['free_electrical_parameters'][0]
    target_line=(f"Target error = {report['target_error_percent']:+.2f}%" if report.get('target_error_percent') is not None
                 else 'Interaction diagnostic')
    text=f"{key} = {report['value']:.6g}\n{target_line}\nComplex S21 RMS = {report['S21_complex_RMSE']:.4g}"
    ax.text(.035,.95,text,transform=ax.transAxes,va='top',fontsize=13,
            bbox=dict(boxstyle='round,pad=.5',fc='white',ec='gray'))
    ax.set_ylim(0,max(abs(actual[:,1,0]).max(),abs(model[:,1,0]).max())*1.45)
    ax.set_xlabel('Frequency (GHz)',fontsize=16);ax.set_ylabel(r'$|S_{21}|$',fontsize=18)
    ax.set_title(title,fontsize=17,pad=12);ax.tick_params(labelsize=13);ax.legend(fontsize=12);ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(out/'s21_fit_linear.png',dpi=200);plt.close(fig)
    fig,axes=plt.subplots(2,2,figsize=(13,8),sharex=True)
    for ax,(i,j) in zip(axes.flat,[(1,0),(0,0),(1,1),(0,1)]):
        ax.plot(f/1e9,actual[:,i,j].real,'r',label='EM real')
        ax.plot(f/1e9,model[:,i,j].real,'r--',label='Model real')
        ax.plot(f/1e9,actual[:,i,j].imag,'b',label='EM imag')
        ax.plot(f/1e9,model[:,i,j].imag,'b--',label='Model imag')
        ax.set_title(f'S{i+1}{j+1}');ax.grid(alpha=.25);ax.legend(fontsize=8)
    fig.tight_layout();fig.savefig(out/'complex_diagnostics.png',dpi=170);plt.close(fig)


def fit_t1_junction(snp,out,*,reference='unit_cell_N4'):
    if reference=='unit_cell_N4':rc,go=CR,GOUT
    elif reference=='prior_T1_fit':rc,go=FIXED['C_regular_fF'],FIXED['Cg_output_fF']
    else:raise ValueError('Unknown fixed regular reference')
    net=read_snp(snp);f=net.f;actual=net.s
    peaks,_=find_peaks(abs(actual[:,1,0]),prominence=.02)
    if len(peaks)!=1:raise ValueError(f'Directly fed T1 should expose one regular-cell peak; found {len(peaks)}')
    wp=2*np.pi*f[peaks[0]]
    seed=1/(wp*wp*LR*1e-9)*1e15-rc-go
    fraction=FIXED['fixture_delay_left_ps']/(FIXED['fixture_delay_left_ps']+FIXED['fixture_delay_right_ps'])
    def predict(x,frequency=f):
        return circuit(frequency,[LT1,LR,LR],[CT1,rc,rc],[None,x[0],go],
                       left_ps=x[1]*fraction,right_ps=x[1]*(1-fraction))
    def residual(x):
        delta=predict(x)[:,1,0]-actual[:,1,0]
        return np.r_[delta.real,delta.imag]
    candidates=[]
    for g in (np.clip(seed,2.01,19.99),8.2865,6.,11.):
        q=least_squares(residual,[g,5.],bounds=([2.,0.],[20.,20.]),
                        x_scale=[8.,5.],max_nfev=250,ftol=1e-11,xtol=1e-11,gtol=1e-11)
        candidates.append(q)
    best=min(candidates,key=lambda q:q.fun@q.fun);prediction=predict(best.x)
    report=dict(free_electrical_parameters=['Cg_T1_to_regular_fF'],value=float(best.x[0]),target=8.2865,
        target_error_percent=float(100*(best.x[0]/8.2865-1)),free_fixture_parameters=['total_reference_delay_ps'],
        total_reference_delay_ps=float(best.x[1]),delay_split_fixed_fraction=float(fraction),
        fixed=dict(L_T1_nH=LT1,C0_T1_fF=CT1,L_regular_nH=LR,C0_regular_fF=rc,Cg_output_fF=go,
                   output_end_LC_included=True),
        fixed_regular_reference=reference,
        fixed_parameter_source=str(UNIT_RESULT) if reference=='unit_cell_N4' else 'geometry/taper/selected_design.json primary_fit',
        S21_complex_RMSE=float(np.sqrt(np.mean(abs(prediction[:,1,0]-actual[:,1,0])**2))),
        S21_magnitude_RMSE=float(np.sqrt(np.mean((abs(prediction[:,1,0])-abs(actual[:,1,0]))**2))),
        S11_complex_RMSE=float(np.sqrt(np.mean(abs(prediction[:,0,0]-actual[:,0,0])**2))),
        S22_complex_RMSE=float(np.sqrt(np.mean(abs(prediction[:,1,1]-actual[:,1,1])**2))),
        peak_GHz=float(f[peaks[0]]/1e9),success=bool(best.success),active_bounds=best.active_mask.tolist(),
        source=str(Path(snp).resolve()),samples=len(f),method='Complex S21; one electrical parameter plus bounded fixture delay',
        acceptance_pending=True)
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    (out/'junction_fit.json').write_text(json.dumps(report,indent=2)+'\n')
    save_plot(f,actual,prediction,report,out,'Directly Fed T1: Downstream Coupling Extraction')
    print(json.dumps(report,indent=2),flush=True)
    return report


def fit_t1_restored(snp,out,junction_ff,c0_fixed=CT1):
    """Reuse unchanged T1 C0/L; measure input Cg with downstream Cg frozen."""
    net=read_snp(snp)
    mask=(net.f>=4.75e9)&(net.f<=5.3e9)
    f,actual=net.f[mask],net.s[mask]
    fraction=FIXED['fixture_delay_left_ps']/(FIXED['fixture_delay_left_ps']+FIXED['fixture_delay_right_ps'])
    def predict(x):
        return circuit(f,[LT1,LR,LR],[c0_fixed,CR,CR],[x[0],junction_ff,GOUT],
                       left_ps=x[1]*fraction,right_ps=x[1]*(1-fraction))
    def residual(x):
        delta=predict(x)[:,1,0]-actual[:,1,0]
        return np.r_[delta.real,delta.imag]
    fits=[least_squares(residual,[g,5],bounds=([20,0],[100,20]),x_scale=[56,5],
        max_nfev=300,ftol=1e-11,xtol=1e-11,gtol=1e-11) for g in (45,56.7,65)]
    best=min(fits,key=lambda q:q.fun@q.fun);prediction=predict(best.x)
    report=dict(free_electrical_parameters=['Cg_input_T1_fF'],value=float(best.x[0]),target=56.7,
        target_error_percent=float(100*(best.x[0]/56.7-1)),
        free_fixture_parameters=['total_reference_delay_ps'],total_reference_delay_ps=float(best.x[1]),
        fixed=dict(L_T1_nH=LT1,C0_T1_fF=c0_fixed,Cg_T1_to_regular_fF=junction_ff,
                   L_regular_nH=LR,C0_regular_fF=CR,Cg_output_fF=GOUT,output_end_LC_included=True),
        fixed_C0_source='Prior selected T1 response; same complete unit geometry, reused as known value',
        source=str(Path(snp).resolve()),samples=len(f),
        S21_complex_RMSE=float(np.sqrt(np.mean(abs(prediction[:,1,0]-actual[:,1,0])**2))),
        S21_magnitude_RMSE=float(np.sqrt(np.mean((abs(prediction[:,1,0])-abs(actual[:,1,0]))**2))),
        S11_complex_RMSE=float(np.sqrt(np.mean(abs(prediction[:,0,0]-actual[:,0,0])**2))),
        S22_complex_RMSE=float(np.sqrt(np.mean(abs(prediction[:,1,1]-actual[:,1,1])**2))),
        active_bounds=best.active_mask.tolist(),acceptance_pending=True)
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    (out/'input_coupling_fit.json').write_text(json.dumps(report,indent=2)+'\n')
    save_plot(f,actual,prediction,report,out,'T1 Input Coupling: Measured Downstream Coupling Fixed')
    print(json.dumps(report,indent=2),flush=True)
    return report


if __name__=='__main__':
    import sys
    fit_t1_junction(sys.argv[1],sys.argv[2])
