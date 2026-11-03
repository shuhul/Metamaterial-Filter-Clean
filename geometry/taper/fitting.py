"""Nonuniform two-resonator circuit fitting, with independent couplings."""
from pathlib import Path
import json
import numpy as np
from scipy.optimize import least_squares, minimize_scalar
from scipy.signal import find_peaks
from circuit_common import read_snp, magnitude_db
import matplotlib.pyplot as plt
from geometry.specifications import TAPER_1, UNIT_CELL


def plot_linear_s21(net, model, result, path):
    """Show linear transmission and signed errors against the taper-1 targets."""
    targets = {'L0_nH': TAPER_1.L0_nH, 'C0_fF': TAPER_1.C0_fF,
               'Cg_input_fF': TAPER_1.Cg_input_fF, 'Cg_junction_fF': UNIT_CELL.Cg_input_fF}
    parameters = [
        (r'$L_{t1}$', 'L_taper_nH', 'L0_nH', 'nH', 4),
        (r'$C_{0,t1}$', 'C_taper_fF', 'C0_fF', 'fF', 3),
        (r'$C_{g,t1}$', 'Cg_input_fF', 'Cg_input_fF', 'fF', 3),
        (r'$C_g$', 'Cg_junction_fF', 'Cg_junction_fF', 'fF', 4),
    ]
    lines = []
    for label, key, target_key, unit, digits in parameters:
        value = result[key]
        error = 100 * (value / targets[target_key] - 1)
        lines.append(f'{label} = {value:.{digits}f} {unit} ({error:+.2f}%)')
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.plot(net.f / 1e9, abs(net.s[:, 1, 0]), 'r-', lw=2.3, label='Sonnet simulation')
    ax.plot(net.f / 1e9, abs(model[:, 1, 0]), 'k--', lw=2.3, label='Circuit fit')
    # Reserve space above the peaks so the large parameter box never hides them.
    peak = max(abs(net.s[:, 1, 0]).max(), abs(model[:, 1, 0]).max())
    ax.set_ylim(0, peak * 1.65)
    ax.text(.04, .95, '\n'.join(lines), transform=ax.transAxes,
            ha='left', va='top', fontsize=13, linespacing=1.35,
            bbox=dict(boxstyle='round,pad=0.45', facecolor='white', edgecolor='gray', alpha=.95))
    ax.set_xlabel('Frequency (GHz)', fontsize=15)
    ax.set_ylabel(r'$|S_{21}|$', fontsize=16)
    ax.set_title('First Taper Cell: Sonnet Simulation and Circuit Fit', fontsize=17, pad=12)
    ax.tick_params(axis='both', labelsize=13)
    ax.legend(loc='upper right', fontsize=13, framealpha=.95)
    ax.grid(linestyle=':', alpha=.65)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def chain(omega, inductances_nh, caps_ff, couplings_ff, z0=50):
    if len(inductances_nh) != len(caps_ff) or len(couplings_ff) != len(caps_ff)+1:
        raise ValueError('N resonators require N+1 couplings')
    m=np.tile(np.eye(2,dtype=complex),(len(omega),1,1))
    def series(c):
        a=np.tile(np.eye(2,dtype=complex),(len(omega),1,1))
        a[:,0,1]=1/(1j*omega*c*1e-15)
        return a
    m=m@series(couplings_ff[0])
    for l,c,g in zip(inductances_nh,caps_ff,couplings_ff[1:]):
        a=np.tile(np.eye(2,dtype=complex),(len(omega),1,1))
        a[:,1,0]=1j*omega*c*1e-15+1/(1j*omega*l*1e-9)
        m=m@a@series(g)
    a,b,c,d=m[:,0,0],m[:,0,1],m[:,1,0],m[:,1,1]
    den=a+b/z0+c*z0+d
    s=np.empty_like(m)
    s[:,0,0]=(a+b/z0-c*z0-d)/den
    s[:,1,1]=(-a+b/z0-c*z0+d)/den
    s[:,1,0]=s[:,0,1]=2/den
    return s


def fit(snp,output,l_taper=6.048402335624231,l_regular=5.014942955,c_regular=191.192,free_regular=False):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    net=read_snp(snp);idx,_=find_peaks(abs(net.s[:,1,0]),prominence=.005)
    if len(idx)!=2:raise ValueError(f'Expected two resolved peaks; found {len(idx)}')
    mask=(net.f>=net.f[idx[0]]-.35e9)&(net.f<=net.f[idx[-1]]+.35e9)
    indices=np.flatnonzero(mask)
    if len(indices)>901:indices=indices[np.unique(np.linspace(0,len(indices)-1,901).astype(int))]
    omega=2*np.pi*net.f[indices];observed=net.s[indices,1,0]
    def response(x,w,delay_fraction=.5):
        s=chain(w,[l_taper,l_regular],[x[0],x[4] if free_regular else c_regular],x[1:4])
        if len(x)>(5 if free_regular else 4):
            # Equal matched delay sections at both ideal circuit ports.
            # This represents smooth residual fixture phase, not a geometric measurement.
            tau=x[-1]*1e-12
            s[:,0,0]*=np.exp(-2j*w*tau*delay_fraction)
            s[:,1,1]*=np.exp(-2j*w*tau*(1-delay_fraction))
            s[:,1,0]*=np.exp(-1j*w*tau)
            s[:,0,1]*=np.exp(-1j*w*tau)
        return s
    lower=[40,10,2,2]+([160] if free_regular else [])
    upper=[220,130,25,25]+([220] if free_regular else [])
    seed=[100,56.7,8.2865,8.2865]+([c_regular] if free_regular else [])
    rng=np.random.default_rng(42)
    fits={}
    for kind in ['amplitude','complex','complex_with_delay']:
        def residual(x):
            prediction=response(x,omega)[:,1,0]
            if kind=='amplitude':return abs(prediction)-abs(observed)
            delta=prediction-observed
            return np.r_[delta.real,delta.imag]
        starts=([seed]+[rng.uniform(lower,upper) for _ in range(3)]) if kind=='amplitude' else [amplitude_best,seed]
        fit_lower,fit_upper=lower,upper
        if kind=='complex_with_delay':
            starts=[list(complex_best)+[3.2],list(amplitude_best)+[6.0]]
            fit_lower,fit_upper=lower+[0],upper+[15]
        results=[least_squares(residual,x,bounds=(fit_lower,fit_upper),max_nfev=600,x_scale='jac',ftol=1e-10) for x in starts]
        best=min(results,key=lambda r:r.fun@r.fun);x=best.x
        if kind=='amplitude':amplitude_best=x.copy()
        if kind=='complex':complex_best=x.copy()
        delay_fraction=.5
        if kind=='complex_with_delay' and x[-1]>1e-6:
            def reflection_error(fraction):
                v=response(x,omega,fraction)
                return np.mean(abs(v[:,0,0]-net.s[indices,0,0])**2+abs(v[:,1,1]-net.s[indices,1,1])**2)
            split=minimize_scalar(reflection_error,bounds=(0,1),method='bounded')
            delay_fraction=min([0.,1.,float(split.x)],key=reflection_error)
        model=response(x,2*np.pi*net.f,delay_fraction)
        fit_mask=model[mask];data=net.s[mask]
        result={'L_taper_nH':l_taper,'L_regular_nH':l_regular,'C_taper_fF':float(x[0]),
                'Cg_input_fF':float(x[1]),'Cg_junction_fF':float(x[2]),'Cg_output_fF':float(x[3]),
                'C_regular_fF':float(x[4] if free_regular else c_regular),'regular_capacitance_fitted':free_regular,
                'effective_fixture_delay_ps':float(x[-1]) if kind=='complex_with_delay' else 0.0,
                'fixture_delay_left_ps':float(x[-1]*delay_fraction) if kind=='complex_with_delay' else 0.0,
                'fixture_delay_right_ps':float(x[-1]*(1-delay_fraction)) if kind=='complex_with_delay' else 0.0,
                'delay_split_from_reflections':kind=='complex_with_delay',
                'delay_model':'Matched port delays; total from S21, split from reflections; effective fixture phase' if kind=='complex_with_delay' else 'none',
                'success':bool(best.success),'starts':len(starts),'sample_count':len(indices),
                'peaks_GHz':(net.f[idx]/1e9).tolist(),'objective':kind,
                's21_amplitude_rmse':float(np.sqrt(np.mean((abs(fit_mask[:,1,0])-abs(data[:,1,0]))**2))),
                's21_complex_rmse':float(np.sqrt(np.mean(abs(fit_mask[:,1,0]-data[:,1,0])**2))),
                's11_complex_rmse':float(np.sqrt(np.mean(abs(fit_mask[:,0,0]-data[:,0,0])**2))),
                's22_complex_rmse':float(np.sqrt(np.mean(abs(fit_mask[:,1,1]-data[:,1,1])**2)))}
        sv=np.linalg.svd(best.jac*x[None,:],compute_uv=False)
        result['relative_jacobian_condition']=float(sv[0]/sv[-1]);fits[kind]=result
        print(f'{kind}: {result}',flush=True)
        (output/'fit.json').write_text(json.dumps(fits,indent=2)+'\n')
        fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True)
        axes[0].plot(net.f/1e9,magnitude_db(net.s[:,1,0]),'k',label='Sonnet')
        axes[0].plot(net.f/1e9,magnitude_db(model[:,1,0]),'r--',label='Independent-coupling fit')
        axes[0].text(.02,.08,f'L={l_taper:.4f} nH; C0={x[0]:.3f} fF; input={x[1]:.3f} fF; junction={x[2]:.3f} fF',transform=axes[0].transAxes,fontsize=8,bbox=dict(facecolor='white',alpha=.85))
        axes[0].set(ylabel='S21 (dB)',ylim=(-65,2),title=f'{Path(snp).stem}: {kind} fit');axes[0].legend()
        axes[1].plot(net.f/1e9,abs(net.s[:,1,0]),'k');axes[1].plot(net.f/1e9,abs(model[:,1,0]),'r--')
        axes[1].set(xlabel='Frequency (GHz)',ylabel='|S21|')
        for ax in axes:ax.grid(alpha=.3)
        fig.tight_layout();fig.savefig(output/f's21_fit_{kind}.png',dpi=160);plt.close(fig)
        plot_linear_s21(net, model, result, output/f's21_fit_{kind}_linear.png')
    (output/'fit.json').write_text(json.dumps(fits,indent=2)+'\n')
    print(json.dumps(fits,indent=2),flush=True)
    return fits
