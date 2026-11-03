"""Offline audit of ideal-inductor loading, topology and fit comparisons.

Does not launch Sonnet or change production parameters. In particular, an
S21-only fit is not a complete equivalent two-port under a new termination.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
from pathlib import Path
import numpy as np
from scipy.optimize import least_squares
from circuit_common import read_snp
from geometry.taper.compact_stages import ROOT
from geometry.taper.ideal_wire_control import load_inductor
from geometry.taper.compact_fit import direct_predictor
from geometry.taper.fitting import chain
from geometry.taper.locked_taper2 import VERIFIED
from geometry.taper.staged_fitting import circuit, LT1,CT1,LR,CR,GOUT
from geometry.taper.response_cases import W
from geometry.taper.run_bridged_stage import digest

OUT=ROOT/'output_data/taper_exploration/ideal_inductor_audit_2026-10-23'
VALUES=(None,5.,3.,2.,1.,.74)


def metrics(actual,model):
    d=model-actual
    return dict(S21_complex_RMS=float(np.sqrt(np.mean(abs(d[:,1,0])**2))),
                S21_complex_max=float(abs(d[:,1,0]).max()),
                S21_magnitude_RMS=float(np.sqrt(np.mean((abs(model[:,1,0])-abs(actual[:,1,0]))**2))),
                S11_complex_RMS=float(np.sqrt(np.mean(abs(d[:,0,0])**2))),
                S22_complex_RMS=float(np.sqrt(np.mean(abs(d[:,1,1])**2))))


def input_shunt(s,f,l):
    """Independent ABCD shunt operation written directly in S parameters."""
    if l is None:return s.copy()
    y=50/(2j*np.pi*f*l*1e-9)
    den=2+y*(1+s[:,0,0])
    q=np.empty_like(s)
    q[:,1,0]=2*s[:,1,0]/den;q[:,0,1]=2*s[:,0,1]/den
    q[:,0,0]=(2*s[:,0,0]-y*(1+s[:,0,0]))/den
    q[:,1,1]=s[:,1,1]-y*s[:,1,0]*s[:,0,1]/den
    return q


def synthetic_operator(f):
    """Independent nodal circuit with two external ports and load terminals.

    Nodes: T2, T1, regular, output end. Port 3 shares the T2 node; port 4
    connects to ideal reference ground. Solve nodal equations with each wave
    port terminated in 50 ohms, rather than building a cascade of ABCD blocks.
    """
    c=np.diag([130,CT1,CR,CR])*1e-15
    for a,b,g in ((0,1,56.7),(1,2,8.40642914296533),(2,3,GOUT)):
        value=g*1e-15;c[a,a]+=value;c[b,b]+=value;c[a,b]-=value;c[b,a]-=value
    w=2*np.pi*f
    y=1j*w[:,None,None]*c[None]
    for i,l in enumerate((None,LT1,LR,LR)):
        if l is not None:y[:,i,i]+=1/(1j*w*l*1e-9)
    b=np.zeros((4,4));b[0,0]=b[3,1]=b[0,2]=1
    return .04*b.T[None]@np.linalg.solve(y+(b@b.T)[None]/50,b[None])-np.eye(4)[None]


def validate_math():
    f=np.linspace(4.65e9,5.4e9,1501)
    s=synthetic_operator(f)
    rows=[]
    for l in (.74,2.,5.):
        actual=load_inductor(s,f,l)
        expected=circuit(f,[l,LT1,LR,LR],[130,CT1,CR,CR],[None,56.7,8.40642914296533,GOUT])
        error=float(abs(actual-expected).max())
        predict=direct_predictor(f,l)
        observed=actual[:,1,0]*np.exp(-2j*np.pi*f*6e-12)
        def residual(x):
            d=predict(x)-observed;return np.r_[d.real,d.imag]
        fit=least_squares(residual,[100,56,5],bounds=([0,5,-20],[2000,180,30]),
                          x_scale=[100,57,5],max_nfev=400,ftol=1e-13,xtol=1e-13,gtol=1e-13)
        rows.append(dict(L_nH=l,all_S_loading_error=error,recovered=fit.x.tolist(),
                         S21_recovery_RMS=float(np.sqrt(np.mean(abs(predict(fit.x)-observed)**2)))))
    if max(r['all_S_loading_error'] for r in rows)>1e-10 or max(r['S21_recovery_RMS'] for r in rows)>1e-8:
        raise ValueError('Ideal loading / synthetic recovery failed')
    return rows


def earlier_t1_loading():
    """Same data, same grid, same port: load the measured and fitted two-ports.

    This is deliberately an input-reference-plane shunt, not a claim that the
    new T2 geometry is identical to the old fixture.
    """
    v=VERIFIED['primary_fit'];source=Path(VERIFIED['artifacts']['touchstone'])
    n=read_snp(source);ix=(n.f>=4.65e9)&(n.f<=5.4e9)
    f=n.f[ix];actual=n.s[ix];w=2*np.pi*f
    model=chain(w,[v['L_taper_nH'],v['L_regular_nH']],
                [v['C_taper_fF'],v['C_regular_fF']],
                [v['Cg_input_fF'],v['Cg_junction_fF'],v['Cg_output_fF']])
    model[:,1,0]*=np.exp(-1j*w*v['effective_fixture_delay_ps']*1e-12)
    model[:,0,1]=model[:,1,0]
    model[:,0,0]*=np.exp(-2j*w*v['fixture_delay_left_ps']*1e-12)
    model[:,1,1]*=np.exp(-2j*w*v['fixture_delay_right_ps']*1e-12)
    rows=[];arrays=dict(f=f,baseline_actual=actual,baseline_model=model)
    for l in VALUES:
        a=input_shunt(actual,f,l);m=input_shunt(model,f,l)
        row=dict(L_nH=l,**metrics(a,m))
        if l is not None:
            # Replace only S11 to isolate the old input-reflection error.
            hybrid=model.copy();hybrid[:,1,0]=actual[:,1,0]
            q=input_shunt(hybrid,f,l)
            row['reflection_error_alone_S21_RMS']=float(np.sqrt(np.mean(abs(q[:,1,0]-a[:,1,0])**2)))
        rows.append(row)
        label='open' if l is None else f'L{l:g}'
        arrays[f'{label}_actual']=a;arrays[f'{label}_model']=m
    np.savez(OUT/'earlier_t1_loading.npz',**arrays)
    return dict(source=str(source),source_sha256=digest(source),samples=len(f),
                interval_GHz=[float(f[0]/1e9),float(f[-1]/1e9)],
                production_parameters_changed=False,rows=rows)


def same_geometry_baseline():
    """Compare the actual wire-absent case and ideal loaded cases equally."""
    source=ROOT/'sonnet_files/taper_exploration/ideal_wire_diagnostic/wire_terminal_control.snp'
    n=W['read_network'](source,4);eye=np.eye(2)[None]
    opened=n.s[:,:2,:2]+n.s[:,:2,2:]@np.linalg.solve(eye-n.s[:,2:,2:],n.s[:,2:,:2])
    regular_refs=dict(L_T1_nH=LT1,C_T1_fF=CT1,L_regular_nH=LR,C_regular_fF=CR,
                      junction_fF=8.40642914296533,output_fF=GOUT)
    data=[]
    for l in VALUES:
        label='open' if l is None else f'L{l:g}'
        actual=opened if l is None else load_inductor(n.s,n.f,l)
        folder=OUT/'same_geometry'/label;folder.mkdir(parents=True,exist_ok=True)
        # Full established fit routine, including global seed search.
        import skrf
        from geometry.taper.compact_fit import fit_direct
        path=folder/'diagnostic.s2p'
        if not path.exists():
            skrf.Network(frequency=skrf.Frequency.from_f(n.f,unit='hz'),s=actual,z0=50).write_touchstone(str(path))
        fitfile=folder/'fit/fit.json'
        r=json.loads(fitfile.read_text()) if fitfile.exists() else fit_direct(path,folder/'fit',1e12 if l is None else l,
                           data_label='Same EM geometry: '+('no inductor' if l is None else f'ideal {l:g} nH'))
        data.append(dict(L_nH=l,S21_complex_RMS=r['S21_complex_RMSE'],S21_complex_max=r['S21_max_complex_error'],
                         bare_C_fF=r['Cparasitic_bare_T2_fF'],side_Cg_fF=r['Cg_T2_T1_fF'],
                         delay_ps=r['fixture_delay_ps'],active_bounds=r['active_bounds']))
    return dict(operator=str(source),operator_sha256=digest(source),
                fixed_downstream=regular_refs,samples=len(n.f),rows=data)


def ground_reference_control():
    """Diagnostic: local bottom-metal terminal is not ideal reference ground.

    Compare the floating load between ports 3/4 with a one-port ideal load at
    port 3 relative to the network's reference, leaving port 4 open. The latter
    is an artificial grounding control, not a manufacturable replacement.
    """
    import skrf
    from geometry.taper.compact_fit import fit_direct
    source=ROOT/'sonnet_files/taper_exploration/ideal_wire_diagnostic/wire_terminal_control.snp'
    n=W['read_network'](source,4);eye=np.eye(4)[None]
    z=50*(eye+n.s)@np.linalg.inv(eye-n.s)
    summary=[]
    for frequency in (4.65,4.9,5.1,5.4):
        i=int(np.argmin(abs(n.f-frequency*1e9)))
        summary.append(dict(GHz=float(n.f[i]/1e9),Z44_real_ohm=float(z[i,3,3].real),
                            Z44_imag_ohm=float(z[i,3,3].imag),
                            apparent_L44_nH=float(z[i,3,3].imag/(2*np.pi*n.f[i])*1e9)))
    rows=[]
    for l in (5.,2.,.74):
        folder=OUT/'global_ground_control'/f'L{l:g}';folder.mkdir(parents=True,exist_ok=True)
        load=np.zeros((len(n.f),2,2),complex)
        zz=2j*np.pi*n.f*l*1e-9
        load[:,0,0]=(zz-50)/(zz+50);load[:,1,1]=1
        actual=n.s[:,:2,:2]+n.s[:,:2,2:]@load@np.linalg.solve(np.eye(2)[None]-n.s[:,2:,2:]@load,n.s[:,2:,:2])
        path=folder/'diagnostic_global_ground.s2p'
        if not path.exists():skrf.Network(frequency=skrf.Frequency.from_f(n.f,unit='hz'),s=actual,z0=50).write_touchstone(str(path))
        fitfile=folder/'fit/fit.json'
        r=json.loads(fitfile.read_text()) if fitfile.exists() else fit_direct(path,folder/'fit',l,
            data_label=f'Ideal {l:g} nH to network reference (diagnostic)')
        rows.append(dict(L_nH=l,S21_complex_RMS=r['S21_complex_RMSE'],S21_complex_max=r['S21_max_complex_error'],
                         C_fF=r['Cparasitic_bare_T2_fF'],Cg_fF=r['Cg_T2_T1_fF'],active_bounds=r['active_bounds']))
    return dict(Z44_condition='Other three ports open; driving-point impedance, not a uniquely isolated lumped ground inductance',
                ground_terminal_impedance=summary,global_ground_fit=rows,
                physical_design=False,production_parameters_changed=False)


def plots(report):
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(17,5.8))
    for ax,l in zip(axes,(None,5,.74)):
        label='open' if l is None else f'L{l:g}'
        a=np.load(OUT/'same_geometry'/label/'fit/response.npz')
        r=next(q for q in report['same_geometry']['rows'] if q['L_nH']==l)
        ax.plot(a['f']/1e9,abs(a['actual'][:,1,0]),color='#b22222',lw=2.6,label='EM layout + load')
        ax.plot(a['f']/1e9,abs(a['model'][:,1,0]),'k--',lw=2.1,label='Frozen downstream fit')
        title='No inductor (open)' if l is None else f'Ideal {l:g} nH'
        ax.set_title(f'{title}\nComplex RMS = {r["S21_complex_RMS"]:.4f}',fontsize=17)
        ax.set_xlabel('Frequency (GHz)',fontsize=15);ax.set_ylabel(r'$|S_{21}|$',fontsize=18)
        ax.tick_params(labelsize=12);ax.grid(alpha=.25);ax.legend(fontsize=10)
    fig.suptitle('Fair Comparison: Same EM Layout, Same Fixed Circuit Parameters',fontsize=21)
    fig.tight_layout();fig.savefig(OUT/'same_geometry_comparison.png',dpi=190);plt.close(fig)
    a=np.load(OUT/'earlier_t1_loading.npz')
    fig,axes=plt.subplots(1,2,figsize=(13,6))
    for ax,label,title in zip(axes,('open','L0.74'),('Earlier T1: Unloaded','Same Two-Port + Ideal 0.74 nH at Input')):
        ax.plot(a['f']/1e9,abs(a[f'{label}_actual'][:,1,0]),color='#b22222',lw=2.6,label='Measured EM two-port')
        ax.plot(a['f']/1e9,abs(a[f'{label}_model'][:,1,0]),'k--',lw=2.1,label='Earlier fitted two-port')
        ax.set_title(title,fontsize=17);ax.set_xlabel('Frequency (GHz)',fontsize=15);ax.set_ylabel(r'$|S_{21}|$',fontsize=18)
        ax.tick_params(labelsize=13);ax.grid(alpha=.25);ax.legend(fontsize=11)
    fig.suptitle('Identical Ideal Loading Exposes Existing Two-Port Model Error',fontsize=20)
    fig.tight_layout();fig.savefig(OUT/'reflection_sensitivity_control.png',dpi=190);plt.close(fig)


def run():
    OUT.mkdir(parents=True,exist_ok=True)
    report=dict(no_new_EM=True,production_parameters_changed=False,
                synthetic_loading_and_recovery=validate_math(),
                earlier_T1=earlier_t1_loading(),same_geometry=same_geometry_baseline(),
                ground_reference=ground_reference_control())
    (OUT/'audit.json').write_text(json.dumps(report,indent=2)+'\n')
    plots(report)
    print('AUDIT COMPLETE',OUT,flush=True)
    return report


if __name__=='__main__':run()
