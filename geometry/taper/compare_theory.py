"""Compare an ordinary physical EM sweep with the fixed theoretical filter.

No gain, frequency shift or circuit value is optimized. One total port delay
may be reported to account for reference planes; magnitude is always raw.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import find_peaks,peak_widths
from scipy.optimize import minimize_scalar
from circuit_common import read_snp
from geometry.taper.optimize_response import theory,THEORY_VERSION

def compare(snp,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    n=read_snp(snp);f=n.f;actual=n.s[:,1,0];target_network=theory(f);target=target_network[:,1,0]
    def features(s):
        ix,_=find_peaks(abs(s),prominence=.02)
        widths=peak_widths(abs(s)**2,ix,rel_height=.5)[0]*(f[1]-f[0])/1e6
        return [{'GHz':float(f[i]/1e9),'height':float(abs(s[i])),
                 'FWHM_MHz':float(w)} for i,w in zip(ix,widths)]
    def cost(tau):return np.mean(abs(actual-target*np.exp(-2j*np.pi*f*tau*1e-12))**2)
    delay=minimize_scalar(cost,bounds=(-30,30),method='bounded')
    aligned=target*np.exp(-2j*np.pi*f*delay.x*1e-12)
    # Transmission fixes the sum of the two reference-plane delays. Check
    # reflections using the same sum, with only its left/right split varied.
    # This diagnostic matters when embedding the fixture in a longer filter.
    def reflection_error(left_ps):
        right_ps=delay.x-left_ps
        s11=target_network[:,0,0]*np.exp(-4j*np.pi*f*left_ps*1e-12)
        s22=target_network[:,1,1]*np.exp(-4j*np.pi*f*right_ps*1e-12)
        return np.r_[n.s[:,0,0]-s11,n.s[:,1,1]-s22]
    split=minimize_scalar(lambda t:np.mean(abs(reflection_error(t))**2),
                          bounds=(-50,50),method='bounded')
    reflection_residual=reflection_error(split.x).reshape(2,-1)
    ap,tp=features(actual),features(target)
    errors=[]
    if len(ap)==len(tp):
        for a,t in zip(ap,tp):
            errors.append({k:100*(a[k]/t[k]-1) for k in ('GHz','height','FWHM_MHz')})
    result={'snp':str(Path(snp).resolve()),'comparison':'Fixed nominal circuit from centralized specifications',
            'theory_version':THEORY_VERSION,
            'fit_variables':['additional_reference_plane_delay_ps'],
            'additional_reference_plane_delay_ps':float(delay.x),
            'circuit_values_fitted':False,'gain_fitted':False,'frequency_shift_fitted':False,
            'S21_complex_RMSE_raw':float(np.sqrt(cost(0))),
            'S21_complex_RMSE_delay_aligned':float(np.sqrt(delay.fun)),
            'S21_complex_max_error_delay_aligned':float(abs(actual-aligned).max()),
            'S21_magnitude_RMSE':float(np.sqrt(np.mean((abs(actual)-abs(target))**2))),
            'S21_magnitude_max_error':float(abs(abs(actual)-abs(target)).max()),
            'S11_magnitude_RMSE_unfitted':float(np.sqrt(np.mean((abs(n.s[:,0,0])-abs(target_network[:,0,0]))**2))),
            'S22_magnitude_RMSE_unfitted':float(np.sqrt(np.mean((abs(n.s[:,1,1])-abs(target_network[:,1,1]))**2))),
            'diagnostic_extra_left_reference_delay_ps':float(split.x),
            'diagnostic_extra_right_reference_delay_ps':float(delay.x-split.x),
            'S11_complex_RMSE_same_total_delay':float(np.sqrt(np.mean(abs(reflection_residual[0])**2))),
            'S22_complex_RMSE_same_total_delay':float(np.sqrt(np.mean(abs(reflection_residual[1])**2))),
            'physical_EM_peaks':ap,'theory_peaks':tp,'peak_errors_percent':errors,
            'acceptance_pending':True}
    (out/'theory_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    fig,ax=plt.subplots(figsize=(11,6.7))
    ax.plot(f/1e9,abs(actual),color='#b22222',lw=2.7,label='Sonnet simulation')
    ax.plot(f/1e9,abs(target),color='#173f5f',ls='--',lw=2.5,label='Specified theoretical filter')
    ax.set_xlabel('Frequency (GHz)',fontsize=18);ax.set_ylabel(r'$|S_{21}|$',fontsize=20)
    ax.set_title('Second Taper: Transmission and Theoretical Target',fontsize=20,pad=14)
    box=('Theoretical T2 design\n'+r'$L_{t2}=0.74$ nH'+'\n'+r'$C_{0,t2}=100$ fF'+'\n'+r'$C_{g,t2}=459$ fF'+
         f"\nMagnitude RMS = {result['S21_magnitude_RMSE']:.3g}")
    ax.text(.035,.95,box,transform=ax.transAxes,va='top',fontsize=13,
            bbox=dict(boxstyle='round,pad=.6',fc='white',ec='#555',alpha=.95))
    ax.tick_params(labelsize=15);ax.legend(loc='upper right',fontsize=12);ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(out/'taper2_s21_theory_linear.png',dpi=230);plt.close(fig)
    fig,axes=plt.subplots(2,1,figsize=(11,7),sharex=True)
    for ax,transform,label in zip(axes,[np.real,np.imag],['Real S21','Imaginary S21']):
        ax.plot(f/1e9,transform(actual),'r',lw=2,label='Sonnet')
        ax.plot(f/1e9,transform(aligned),'k--',lw=2,label='Theory with reference-plane delay')
        ax.set_ylabel(label);ax.grid(alpha=.25);ax.legend()
    axes[1].set_xlabel('Frequency (GHz)');fig.tight_layout();fig.savefig(out/'taper2_s21_theory_complex.png',dpi=200);plt.close(fig)
    print(json.dumps(result,indent=2),flush=True)
    return result

if __name__=='__main__':
    import sys
    compare(sys.argv[1],sys.argv[2])
