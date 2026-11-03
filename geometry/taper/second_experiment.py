"""Named taper-2 trials: matched L calibration, full EM fixture, complex fit."""
from pathlib import Path
from dataclasses import asdict
from hashlib import sha256
import json
import numpy as np
from scipy.optimize import least_squares, minimize_scalar
from scipy.ndimage import gaussian_filter1d
import matplotlib.pyplot as plt
from geometry.specifications import TAPER_1,TAPER_2,UNIT_CELL,specification
from geometry.taper.second_layout import write,validate,from_dict,require_same_isolation_geometry
from geometry.taper.sonnet import add_capacitor,configure,simulate,solver_environment
from geometry.prepare_sonnet import main as prepare
from extract_generated_inductance import main as extract_l
from geometry.taper.fitting import chain
from circuit_common import read_snp,magnitude_db

ROOT=Path(__file__).resolve().parents[2]
L1=6.014105932449249
LR=4.922723945068658


def plot_final(snp,result,path,*,title='Second Taper Cell: Sonnet Simulation and Circuit Fit'):
    """Linear |S21| with stage-specific notation and signed target errors."""
    net=read_snp(snp);model=model_response(net.f,result['x'],result['L_t2_nH'],l1=result.get('L_t1_nH',L1))
    params=[(r'$L_{t2}$',result['L_t2_nH'],TAPER_2.L0_nH,'nH',4),
            (r'$C_{0,t2}$',result['C0_t2_fF'],TAPER_2.C0_fF,'fF',3),
            (r'$C_{g,t2}$',result['Cg_t2_fF'],TAPER_2.Cg_input_fF,'fF',3),
            (r'$C_{g,t1}$',result['Cg_t1_fF'],TAPER_1.Cg_input_fF,'fF',3),
            (r'$C_g$',result['Cg_fF'],UNIT_CELL.Cg_input_fF,'fF',4)]
    text='\n'.join(f'{label} = {v:.{digits}f} {unit} ({100*(v/t-1):+.2f}%)' for label,v,t,unit,digits in params)
    fig,ax=plt.subplots(figsize=(10,6))
    ax.plot(net.f/1e9,abs(net.s[:,1,0]),'r-',lw=2.3,label='Sonnet simulation')
    ax.plot(net.f/1e9,abs(model[:,1,0]),'k--',lw=2.3,label='Circuit fit')
    peak=max(abs(net.s[:,1,0]).max(),abs(model[:,1,0]).max())
    ax.set_ylim(0,peak*1.7)
    ax.text(.035,.95,text,transform=ax.transAxes,va='top',fontsize=13,linespacing=1.3,
            bbox=dict(boxstyle='round,pad=0.45',facecolor='white',edgecolor='gray',alpha=.95))
    ax.set_xlabel('Frequency (GHz)',fontsize=15);ax.set_ylabel(r'$|S_{21}|$',fontsize=16)
    ax.set_title(title,fontsize=17,pad=12)
    ax.tick_params(labelsize=13);ax.grid(linestyle=':',alpha=.65);ax.legend(loc='upper right',fontsize=12)
    fig.tight_layout();fig.savefig(path,dpi=220);plt.close(fig)


def measure(name,p,*,c_add_sweep_pf=(.1,.2,.3,.4),calibration_band=(3,14),calibration_target=2200):
    case=name+'_isolation';output=ROOT/'output_data/taper_exploration'/case
    if (output/'results.json').exists():
        meta=json.loads((ROOT/'geometry/taper/output'/f'{case}.json').read_text())
        # External t1 settings do not change this calibration geometry.
        if from_dict(meta['second'])!=p:raise ValueError('Existing calibration differs')
        return output/'results.json',json.loads((output/'results.json').read_text())
    dest=ROOT/'sonnet_files/taper_exploration'/case
    project=dest/'generated_ncell_1.sonx'
    gds=ROOT/'geometry/taper/output'/f'{case}.gds'
    if project.exists() and (output/'sweep_points.json').exists():
        meta=json.loads(gds.with_suffix('.json').read_text())
        if from_dict(meta['second'])!=p:raise ValueError('Interrupted calibration geometry differs')
        from geometry.taper.audit import check_geometry
        check_geometry(project,gds,meta,geometry_validator=validate)
        print('Resuming validated calibration with existing EM cache:',case,flush=True)
    else:
        gds=write(p,ROOT/'geometry/taper/output'/case,isolation=True)
        prepare(1,gds_path=gds,output_dir=dest,geometry_validator=validate)
        add_capacitor(project,json.loads(gds.with_suffix('.json').read_text()))
        configure(project,*calibration_band,calibration_target)
    result=extract_l(project_file=project,output_dir=output,fit_capacitance=False,
                     c_add_sweep_pf=c_add_sweep_pf,target_l0_nh=TAPER_2.L0_nH,solver_env=solver_environment())
    return output/'results.json',result


def model_response(f,x,l2,*,l1=L1):
    # x = C02, C01, C0r, Cgin2, Cg21, Cg1r, Cgout, total delay ps.
    s=chain(2*np.pi*f,[l2,l1,LR],x[:3],x[3:7])
    total=x[7]+x[8] if len(x)==9 else x[7]
    s[:,1,0]*=np.exp(-2j*np.pi*f*total*1e-12)
    s[:,0,1]=s[:,1,0]
    if len(x)==9:
        s[:,0,0]*=np.exp(-4j*np.pi*f*x[7]*1e-12)
        s[:,1,1]*=np.exp(-4j*np.pi*f*x[8]*1e-12)
    return s


def fit_network(snp,output,l2,*,l1=L1):
    """Diagnostic joint complex S21/S11 fit: the strongly loaded t2 is visible in reflection."""
    net=read_snp(snp);output=Path(output)
    indices=np.flatnonzero(net.f<=9e9)[::2]
    f=net.f[indices];obs=net.s[indices]
    lower=np.array([20,60,150,100,10,2,2,0,0])
    upper=np.array([300,160,230,1600,140,25,25,20,20])
    seed=np.array([100,99,186,459,56.7,8.2865,8.2865,5,.6])
    def residual(x):
        pred=model_response(f,x,l2,l1=l1)
        transmission=3*(pred[:,1,0]-obs[:,1,0])
        reflection=pred[:,0,0]-obs[:,0,0]
        return np.r_[transmission.real,transmission.imag,reflection.real,reflection.imag]
    starts=[seed]
    # First align broadened peaks. This avoids narrow-resonance local minima.
    # Only this initialization holds the input seed; every final parameter is free.
    band=np.flatnonzero(net.f<=6e9)
    fb=net.f[band];amplitude=abs(net.s[band,1,0]);spacing=float(np.median(np.diff(fb)))
    for coupling in (25,56.7,90):
        v=seed.copy();v[4]=coupling
        for width_hz in (60e6,20e6,6e6):
            sigma=max(.5,width_hz/spacing)
            observed=gaussian_filter1d(amplitude,sigma)
            def smooth_residual(z):
                trial=v.copy();trial[[1,2,4,5,6]]=z
                prediction=abs(model_response(fb,trial,l2,l1=l1)[:,1,0])
                return gaussian_filter1d(prediction,sigma)-observed
            keys=[1,2,4,5,6]
            init=least_squares(smooth_residual,v[keys],bounds=(lower[keys],upper[keys]),
                               x_scale='jac',max_nfev=150)
            v[keys]=init.x
        starts.append(v)
    # Sharp downstream resonances create local minima if a seed misses a peak.
    # Vary the strong junction and both downstream capacitances independently.
    for coupling in (25,40,70):
        for c1,cr in ((99,186),(101.49,188.45)):
            v=seed.copy();v[1:3]=[c1,cr];v[4]=coupling;starts.append(v)
    previous=output/'fit.json'
    if previous.exists():
        old=json.loads(previous.read_text())['complex_with_delay']['x']
        starts.append(np.clip(np.r_[old[:7],1,1],lower+1e-6,upper-1e-6))
    starts=sorted(starts,key=lambda x:float(residual(x)@residual(x)))[:3]
    results=[least_squares(residual,x,bounds=(lower,upper),x_scale='jac',max_nfev=500,ftol=1e-10) for x in starts]
    r=min(results,key=lambda r:r.fun@r.fun);x=r.x;pred=model_response(f,x,l2,l1=l1)
    sv=np.linalg.svd(r.jac*np.maximum(abs(x),1)[None,:],compute_uv=False)
    result={'x':x.tolist(),'L_t2_nH':l2,'L_t1_nH':l1,'L_regular_nH':LR,
            'C0_t2_fF':float(x[0]),'C0_t1_fF':float(x[1]),'C0_regular_fF':float(x[2]),
            'Cg_t2_fF':float(x[3]),'Cg_t1_fF':float(x[4]),'Cg_fF':float(x[5]),'Cg_output_fF':float(x[6]),
            'delay_left_ps':float(x[7]),'delay_right_ps':float(x[8]),'delay_ps':float(x[7]+x[8]),
            'objective':'joint complex S21/S11, weights 3:1, independently fitted matched port delays',
            'fit_range_GHz':[float(f[0]/1e9),float(f[-1]/1e9)],'success':bool(r.success),
            's21_complex_rmse':float(np.sqrt(np.mean(abs(pred[:,1,0]-obs[:,1,0])**2))),
            's11_complex_rmse':float(np.sqrt(np.mean(abs(pred[:,0,0]-obs[:,0,0])**2))),
            'active_bounds':r.active_mask.tolist(),'relative_jacobian_condition':float(sv[0]/sv[-1])}
    result['starts_optimized']=len(starts)
    (output/'fit_network.json').write_text(json.dumps(result,indent=2)+'\n')
    print('NETWORK FIT',json.dumps(result),flush=True)
    return result


def fit_response(snp,output,l2,*,l1=L1,seed_result=None):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    net=read_snp(snp)
    # Fit the broad t2 feature and both lower-frequency downstream modes.
    low=np.flatnonzero(net.f<=6e9)
    high=np.flatnonzero((net.f>6e9)&(net.f<=9e9))
    # Preserve the narrow downstream peaks; thin only the broad upper feature.
    indices=np.unique(np.r_[low,high[::max(1,len(high)//601)]])
    f=net.f[indices];obs=net.s[indices,1,0]
    lower=[20,60,150,100,15,2,2,0]
    upper=[250,160,230,1200,140,25,25,40]
    seed=[100,101.492,188.451,459,56.7,8.2865,8.2865,5]
    if seed_result is not None:
        seed=np.r_[np.clip(seed_result['x'][:7],np.array(lower[:7])+1e-6,np.array(upper[:7])-1e-6),seed_result['delay_ps']].tolist()
    rng=np.random.default_rng(42)
    fits={};best_x=None
    for objective in ('amplitude','complex_with_delay'):
        def residual(x):
            pred=model_response(f,x,l2,l1=l1)[:,1,0]
            if objective=='amplitude':return abs(pred)-abs(obs)
            delta=pred-obs
            return np.r_[delta.real,delta.imag]
        if objective=='amplitude':
            # Delay is invisible to magnitude: hold it at zero, fit seven variables.
            def amp_res(v):return residual(np.r_[v,0])
            starts=[seed[:7]]+[np.clip(np.array(seed[:7])*rng.uniform(.7,1.3,7),
                                      np.array(lower[:7])+1e-6,np.array(upper[:7])-1e-6) for _ in range(1)]
            if seed_result is not None:starts=starts[:1]
            results=[least_squares(amp_res,x,bounds=(lower[:7],upper[:7]),x_scale='jac',max_nfev=300,ftol=1e-10) for x in starts]
            chosen=min(results,key=lambda r:r.fun@r.fun);best_x=np.r_[chosen.x,0]
        else:
            starts=[seed,list(best_x[:7])+[5],list(best_x[:7])+[15]]
            if seed_result is not None:starts=[seed,list(best_x[:7])+[seed_result['delay_ps']]]
            results=[least_squares(residual,x,bounds=(lower,upper),x_scale='jac',max_nfev=500,ftol=1e-10) for x in starts]
            chosen=min(results,key=lambda r:r.fun@r.fun);best_x=chosen.x
        x=best_x;pred=model_response(net.f,x,l2,l1=l1)
        # The split does not alter S21 or refit any capacitance.
        base=chain(2*np.pi*f,[l2,l1,LR],x[:3],x[3:7])
        def reflection_error(fraction):
            l=np.exp(-4j*np.pi*f*x[7]*1e-12*fraction)
            r=np.exp(-4j*np.pi*f*x[7]*1e-12*(1-fraction))
            return np.mean(abs(base[:,0,0]*l-net.s[indices,0,0])**2+abs(base[:,1,1]*r-net.s[indices,1,1])**2)
        frac=float(minimize_scalar(reflection_error,bounds=(0,1),method='bounded').x) if x[7] else .5
        result={'L_t2_nH':l2,'L_t1_nH':l1,'L_regular_nH':LR,
                'C0_t2_fF':float(x[0]),'C0_t1_fF':float(x[1]),'C0_regular_fF':float(x[2]),
                'Cg_t2_fF':float(x[3]),'Cg_t1_fF':float(x[4]),'Cg_fF':float(x[5]),'Cg_output_fF':float(x[6]),
                'delay_ps':float(x[7]),'delay_left_ps':float(x[7]*frac),'delay_right_ps':float(x[7]*(1-frac)),
                'objective':objective,'success':bool(chosen.success),'cost':float(chosen.fun@chosen.fun),
                's21_complex_rmse':float(np.sqrt(np.mean(abs(pred[:,1,0]-net.s[:,1,0])**2))),
                's21_amplitude_rmse':float(np.sqrt(np.mean((abs(pred[:,1,0])-abs(net.s[:,1,0]))**2))),
                'x':list(map(float,x)),'samples':len(indices),'bounds':[lower,upper],
                'fit_range_GHz':[float(f[0]/1e9),float(f[-1]/1e9)]}
        singular=np.linalg.svd(chosen.jac*chosen.x[None,:],compute_uv=False)
        result['relative_jacobian_condition']=float(singular[0]/singular[-1])
        result['active_bounds']=chosen.active_mask.tolist()
        fits[objective]=result
        fig,axes=plt.subplots(2,1,figsize=(10,7),sharex=True)
        axes[0].plot(net.f/1e9,magnitude_db(net.s[:,1,0]),'r-',label='Sonnet')
        axes[0].plot(net.f/1e9,magnitude_db(pred[:,1,0]),'k--',label='Circuit fit')
        axes[0].set(ylabel='S21 (dB)',ylim=(-70,2),title='Taper 2 + Taper 1 + Regular Cell');axes[0].legend()
        axes[1].plot(net.f/1e9,abs(net.s[:,1,0]),'r-');axes[1].plot(net.f/1e9,abs(pred[:,1,0]),'k--')
        axes[1].set(xlabel='Frequency (GHz)',ylabel='|S21|')
        for ax in axes:ax.grid(alpha=.3)
        fig.tight_layout();fig.savefig(output/f's21_fit_{objective}.png',dpi=160);plt.close(fig)
        print(objective,json.dumps(result),flush=True)
    (output/'fit.json').write_text(json.dumps(fits,indent=2)+'\n')
    plot_final(snp,fits['complex_with_delay'],output/'s21_fit_linear.png')
    return fits


def run(name,p,inductance_reference=None,*,inductance_t1_reference=None):
    output=ROOT/'output_data/taper_exploration'/name
    if (output/'trial.json').exists():raise ValueError('Choose a new case name')
    l1=L1
    if p.taper1 is not None:
        if inductance_t1_reference is None:raise ValueError('Changed t1 requires its own matching L calibration')
        inductance_t1_reference=Path(inductance_t1_reference).resolve()
        t1_meta=json.loads((ROOT/'geometry/taper/output'/f'{inductance_t1_reference.parent.name}.json').read_text())
        if t1_meta['taper']!=asdict(p.taper1):raise ValueError('T1 L calibration geometry differs')
        if t1_meta.get('w_top_slot',0)!=p.w_t1_ground_slot or t1_meta.get('l_top_slot',100)!=p.l_t1_ground_slot:
            raise ValueError('T1 calibration ground slot differs')
        l1=json.loads(inductance_t1_reference.read_text())['L0_nH']
    if inductance_reference is None:inductance_reference,lr=measure(name,p)
    else:
        inductance_reference=Path(inductance_reference)
        # Changing neighbor geometry is allowed; taper and its comb must be identical.
        meta=asdict(from_dict(json.loads((ROOT/'geometry/taper/output'/f'{inductance_reference.parent.name}.json').read_text())['second']))
        expected=asdict(p)
        for key in ('taper','idc_pairs','idc_t','idc_d','idc_span'):
            if expected[key]!=meta[key]:raise ValueError('L reference has different t2 geometry')
        for key in ('d_idc_top','x_meander','t_ground_tail'):
            if expected[key]!=meta[key]:raise ValueError('L reference has different IDC/meander placement')
        for key in ('t_coupling_finger','l_coupling_finger','n_coupling_teeth','w_coupling_teeth',
                    'w_upper_plate','l_coupling_bridge','l_cal_output','w_t1_ground_slot','l_t1_ground_slot'):
            if expected[key]!=meta[key]:raise ValueError('L reference has a different coupling finger')
        if p.t_coupling_finger and expected['d_t2_t1']!=meta['d_t2_t1']:
            raise ValueError('Coupling-finger position differs from L calibration')
        if p.t_coupling_finger and expected['taper1']!=meta['taper1']:
            require_same_isolation_geometry(p,ROOT/'geometry/taper/output'/f'{inductance_reference.parent.name}.gds')
        lr=json.loads(inductance_reference.read_text())
    gds=write(p,ROOT/'geometry/taper/output'/name)
    dest=ROOT/'sonnet_files/taper_exploration'/name
    prepare(3,gds_path=gds,output_dir=dest,geometry_validator=validate)
    project=configure(dest/'generated_ncell_3.sonx',3.8,9,2600)
    snp=simulate(project,output)
    network=fit_network(snp,output,lr['L0_nH'],l1=l1)
    fits=fit_response(snp,output,lr['L0_nH'],l1=l1,seed_result=network)
    primary=network
    plot_final(snp,primary,output/'s21_fit_linear.png')
    comparisons={'L0_t2':(lr['L0_nH'],TAPER_2.L0_nH),'C0_t2':(primary['C0_t2_fF'],TAPER_2.C0_fF),
                 'Cg_t2':(primary['Cg_t2_fF'],TAPER_2.Cg_input_fF),'Cg_t1':(primary['Cg_t1_fF'],TAPER_1.Cg_input_fF),
                 'Cg':(primary['Cg_fF'],UNIT_CELL.Cg_input_fF),
                 'L0_t1':(l1,TAPER_1.L0_nH),
                 'C0_t1':(primary['C0_t1_fF'],TAPER_1.C0_fF),'C0_regular':(primary['C0_regular_fF'],UNIT_CELL.C0_fF)}
    errors={k:100*(v/t-1) for k,(v,t) in comparisons.items()}
    result={'name':name,'parameters':asdict(p),'specification':specification(),
            'inductance_reference':str(inductance_reference),'gds':str(gds),'project':str(project),'snp':str(snp),
            'inductance_t1_reference':str(inductance_t1_reference) if inductance_t1_reference else None,
            'fits':fits,'network_fit':network,'primary_fit':'joint complex S21/S11 with matched port delays',
            'errors_percent':errors,'passes_2_percent':primary['success'] and not any(primary['active_bounds'][:7]) and all(abs(v)<=2 for v in errors.values()),
            'hashes':{str(f):sha256(Path(f).read_bytes()).hexdigest() for f in [gds,project,snp,inductance_reference]+([inductance_t1_reference] if inductance_t1_reference else [])}}
    (output/'trial.json').write_text(json.dumps(result,indent=2)+'\n')
    print('ERRORS',errors,flush=True)
    return result
