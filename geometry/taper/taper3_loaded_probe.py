"""Known native shunt capacitors identify T3 under unchanged EM geometry.

One EM calculation with native CAP terminals; later C_add values reuse cache.
All resonance/load responses come directly from Sonnet, never a Python Y fit.
The production project stays unchanged. Run serially after production EM.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import shutil
import xml.etree.ElementTree as ET
import numpy as np
from scipy.optimize import least_squares
from geometry.taper import taper3 as work
from geometry.taper.sonnet import add_capacitor,configure,simulate
from geometry.taper.progressive_fitting import transmission
from circuit_common import read_snp,write_sonx
from geometry.prepare_sonnet import tag,write_xml,BIN
from geometry.taper.solver_guard import run as guarded_solver
from geometry.taper.solver_resources import solver_environment


CAPS_FF=[.1,100.,250.,400.]


def simulate_loads(case):
    source=work.NATIVE/case/'response/generated_ncell_4.sonx'
    meta=json.loads((source.parent/'geometry.json').read_text())
    native=work.NATIVE/case/'known_Cadd';out=work.OUT/case/'known_Cadd'
    native.mkdir(exist_ok=True);out.mkdir(exist_ok=True)
    project=native/source.name
    if not project.exists():
        shutil.copy2(source,project);add_capacitor(project,meta);configure(project,4.65,5.4,7500)
    manifest=out/'manifest.json'
    expected=dict(source_sha256=work.digest(source),C_add_fF=CAPS_FF,
                  method='Native explicit CAP; one new EM, cache-only later loads',
                  source_project=str(source),mesh_um=[2,2])
    if manifest.exists() and json.loads(manifest.read_text())!=expected:raise ValueError('Probe source/config changed')
    manifest.write_text(json.dumps(expected,indent=2)+'\n')
    for i,cap in enumerate(CAPS_FF):
        folder=out/f'cadd_{cap:g}ff';folder.mkdir(exist_ok=True)
        marker=folder/'completed.json'
        if marker.exists():
            record=json.loads(marker.read_text())
            for name,key in [('project.sonx','project_sha256'),('response.snp','snp_sha256')]:
                if work.digest(folder/name)!=record[key]:raise ValueError('Probe archive changed')
            continue
        write_sonx({'C_add':str(cap/1000)},project)
        print('KNOWN LOAD',cap,'fF',flush=True)
        simulate(project,folder,cache_only=i>0)
        record=dict(C_add_fF=cap,project_sha256=work.digest(folder/'project.sonx'),
                    snp_sha256=work.digest(folder/'response.snp'))
        marker.write_text(json.dumps(record,indent=2)+'\n')
    return out


def fit_loads(case, *, free_L=False, separate_delays=False, older_fraction=.01):
    root=work.OUT/case/'known_Cadd'
    name='shared_fit'+('_freeL_DIAGNOSTIC' if free_L else '')+('_separate_delays_DIAGNOSTIC' if separate_delays else '')
    if older_fraction!=.01:name+=f'_older{100*older_fraction:g}percent_DIAGNOSTIC'
    out=root/name;out.mkdir(exist_ok=True)
    paths=[root/f'cadd_{c:g}ff/response.snp' for c in CAPS_FF]
    nets=[read_snp(p) for p in paths]
    refs=work.REFERENCE_FIT
    ls=[json.loads((work.OUT/case/'inductance.json').read_text())['L0_nH']]+refs['L_nH']
    old=np.array(refs['C_fF']+refs['G_fF'][1:])
    lo=np.r_[0,.1,.1,(1-older_fraction)*old,-100];hi=np.r_[1000,3000,3000,(1+older_fraction)*old,100]
    scale=np.r_[251,261,153,old,20]
    if separate_delays:
        lo=np.r_[lo,[-100]*3];hi=np.r_[hi,[100]*3];scale=np.r_[scale,[20]*3]
    if free_L:
        lo=np.r_[lo,.5];hi=np.r_[hi,6.];scale=np.r_[scale,2.77]
    ix=[np.linspace(0,len(n.f)-1,1501).astype(int) for n in nets]
    def predict(f,x,cap,index=0):
        inductances=[x[-1]]+ls[1:] if free_L else ls
        return transmission(f,inductances,[x[0]+cap]+x[3:7].tolist(),[x[1],x[2]]+x[7:10].tolist(),x[10+index if separate_delays else 10])
    def residual(x,full=False,smooth=False):
        blocks=[]
        for i,(n,cap) in enumerate(zip(nets,CAPS_FF)):
            j=np.arange(len(n.f)) if full else ix[i]
            model=predict(n.f[j],x,cap,i);actual=n.s[j,1,0]
            q=((1/model-1/actual)*np.maximum(abs(actual),.03)**2 if smooth else model-actual)
            blocks.extend([q.real,q.imag])
        return np.concatenate(blocks)
    starts=[]
    for c in (50,176,251,500,900):
        for g in (100,261,500):
            for gj in (100,153,300):
                x=np.r_[c,g,gj,old,5].astype(float)
                if separate_delays:x=np.r_[x,[5.]*3]
                if free_L:x=np.r_[x,ls[0]]
                n=nets[0];peak=np.argmax(abs(n.s[:,1,0]))
                zero=predict(n.f[peak:peak+1],x,CAPS_FF[0])[0]
                correction=-np.angle(n.s[peak,1,0]/zero)/(2*np.pi*n.f[peak])*1e12
                if separate_delays:x[10:14]+=correction
                else:x[10]+=correction
                starts.append(np.clip(x,lo+1e-8,hi-1e-8))
    starts.sort(key=lambda x:residual(x)@residual(x))
    chosen=[]
    for x in starts:
        if not chosen or all(np.linalg.norm((x[:3]-q[:3])/scale[:3])>.25 for q in chosen):chosen.append(x)
        if len(chosen)==12:break
    smooth_fits=[least_squares(lambda x:residual(x,smooth=True),x,bounds=(lo,hi),x_scale=scale,max_nfev=500)
                 for x in chosen]
    seeds=chosen+[q.x for q in smooth_fits]
    seeds.sort(key=lambda x:residual(x)@residual(x))
    candidates=[least_squares(residual,x,bounds=(lo,hi),x_scale=scale,max_nfev=800,
        ftol=1e-11,xtol=1e-11,gtol=1e-11) for x in seeds[:8]]
    best=min(candidates,key=lambda q:q.fun@q.fun)
    best=least_squares(lambda x:residual(x,True),best.x,bounds=(lo,hi),x_scale=scale,max_nfev=800,
        ftol=1e-11,xtol=1e-11,gtol=1e-11)
    report=dict(case=case,method='Shared C0 and couplings across known native Sonnet CAP loads',
        L_nH=([float(best.x[-1])]+ls[1:] if free_L else ls),independent_calibrated_L3_nH=ls[0],
        free_L_diagnostic=free_L,separate_delays_diagnostic=separate_delays,
        older_fractional_bounds=older_fraction,
        C_fF=[float(best.x[0])]+best.x[3:7].tolist(),G_fF=best.x[1:3].tolist()+best.x[7:10].tolist(),
        delay_ps=(best.x[10:14].tolist() if separate_delays else float(best.x[10])),active_bounds=best.active_mask.tolist(),
        older_fractional_changes=(best.x[3:10]/old-1).tolist(),
        optimizer_success=bool(best.success),sources={str(p):work.digest(p) for p in paths},loads=[])
    sensitivity_scale=np.r_[251,261,153,.01*old,20]
    if separate_delays:sensitivity_scale=np.r_[sensitivity_scale,[20]*3]
    if free_L:sensitivity_scale=np.r_[sensitivity_scale,2.77]
    jac=best.jac*sensitivity_scale[None,:]
    singular=np.linalg.svd(jac,compute_uv=False)
    inverse=np.linalg.pinv(jac.T@jac,rcond=1e-12)
    standard=np.sqrt(np.maximum(np.diag(inverse),1e-30))
    correlation=inverse/(standard[:,None]*standard[None,:])
    report['local_sensitivity']=dict(singular_values=singular.tolist(),
        condition_number=float(singular[0]/singular[-1]),
        new_capacitance_correlations=correlation[:3,:3].tolist(),
        note='Scaled local Jacobian diagnostic; not a noise-based confidence interval or global uniqueness proof.')
    for index,(n,cap) in enumerate(zip(nets,CAPS_FF)):
        m=predict(n.f,best.x,cap,index);err=abs(m-n.s[:,1,0])
        report['loads'].append(dict(C_add_fF=cap,RMS=float(np.sqrt(np.mean(err**2))),maximum=float(err.max())))
        np.savez(out/f'cadd_{cap:g}ff.npz',f=n.f,actual=n.s,model=m)
    baseline=read_snp(work.NATIVE/case/'response/generated_ncell_4.snp')
    report['unloaded_2port_vs_0p1fF_probe_RMS']=float(np.sqrt(np.mean(abs(baseline.s[:,1,0]-nets[0].s[:,1,0])**2)))
    report['accepted']=False
    (out/'fit.json').write_text(json.dumps(report,indent=2)+'\n')
    print('LOADED FIT',json.dumps(report),flush=True);return report


def plot_loads(case, fit_name='shared_fit'):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    root=work.OUT/case/'known_Cadd'/fit_name
    report=json.loads((root/'fit.json').read_text())
    fig,axes=plt.subplots(2,2,figsize=(16,11))
    for ax,row in zip(axes.flat,report['loads']):
        cap=row['C_add_fF'];data=np.load(root/f'cadd_{cap:g}ff.npz')
        ax.plot(data['f']/1e9,abs(data['actual'][:,1,0]),color='#2768aa',lw=2.3,label='Native Sonnet CAP')
        ax.plot(data['f']/1e9,abs(data['model']),'k--',lw=2,label='Shared-component fit')
        ax.set_title(f'Added C = {cap:g} fF; Fit RMS = {row["RMS"]:.4f}',fontsize=16)
        ax.set(xlabel='Frequency (GHz)',ylabel=r'$|S_{21}|$');ax.tick_params(labelsize=12)
        ax.xaxis.label.set_size(14);ax.yaxis.label.set_size(18);ax.grid(alpha=.2);ax.legend(fontsize=11)
    fig.suptitle('T3 Known-Capacitor Diagnostic — Not an Accepted Design',fontsize=21)
    fig.tight_layout(rect=[0,0,1,.95]);fig.savefig(root/'loaded_comparison.png',dpi=180);plt.close(fig)
    return root/'loaded_comparison.png'


def direct_check(case, frequency_GHz=None):
    """One fresh frequency validates cached loading at its strongest new peak."""
    source=work.OUT/case/'known_Cadd/cadd_100ff'
    net=read_snp(source/'response.snp')
    index=(int(np.argmax(abs(net.s[:,1,0]))) if frequency_GHz is None else
           int(np.argmin(abs(net.f-frequency_GHz*1e9))))
    frequency=float(net.f[index]/1e9)
    suffix='' if frequency_GHz is None else '_f'+str(frequency).replace('.','p')
    native=work.NATIVE/case/('known_Cadd_direct100'+suffix);native.mkdir(exist_ok=True)
    out=work.OUT/case/'known_Cadd'/('direct100_validation'+suffix);out.mkdir(exist_ok=True)
    project=native/'direct100.sonx';marker=out/'comparison.json'
    if marker.exists():
        record=json.loads(marker.read_text())
        if work.digest(project)!=record['project_sha256'] or work.digest(project.with_suffix('.snp'))!=record['snp_sha256']:
            raise ValueError('Direct validation artifacts changed')
        return record
    tree=ET.parse(source/'project.sonx');root=tree.getroot();sweeps=root.find(tag('Sweeps'))
    for q in list(sweeps):sweeps.remove(q)
    fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
    ET.SubElement(fs,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(frequency),Stop=str(frequency),Step='1')
    for q in root.iter(tag('File')):
        if q.get('Type')=='TOUCH':q.set('Name',project.with_suffix('.snp').name);q.set('UseInterpolatedData','FALSE')
    existing_hash=work.digest(project) if project.exists() else None
    write_xml(tree,project)
    log_path=out/'em.log'
    complete=(project.with_suffix('.snp').exists() and existing_hash==work.digest(project)
              and log_path.exists() and 'Em simulation completed' in log_path.read_text()
              and 'Em errors detected: 0' in log_path.read_text())
    if not complete:
        with log_path.open('w') as log:
            process=guarded_solver([str(BIN/'em.exe'),'-v',str(project)],cwd=native,
                                  stdout=log,env=solver_environment())
        process.check_returncode()
    # The project-wide reader intentionally rejects single-frequency fits.
    # Use scikit-rf only to read this independent one-point comparison.
    import skrf
    direct_path=out/'direct_response.s2p'
    shutil.copy2(project.with_suffix('.snp'),direct_path)
    direct=skrf.Network(str(direct_path))
    if len(direct.f)!=1 or abs(direct.f[0]-net.f[index])>1:raise ValueError('Wrong direct-check frequency')
    error=abs(direct.s[0]-net.s[index])
    record=dict(C_add_fF=100,frequency_GHz=frequency,cache_source_sha256=work.digest(source/'response.snp'),
        project_sha256=work.digest(project),snp_sha256=work.digest(project.with_suffix('.snp')),
        complex_S21_difference=float(error[1,0]),all_S_max_difference=float(error.max()),
        cached_magnitude=float(abs(net.s[index,1,0])),direct_magnitude=float(abs(direct.s[0,1,0])))
    marker.write_text(json.dumps(record,indent=2)+'\n');print('DIRECT LOAD CHECK',record,flush=True);return record


def derivative_check(case):
    """Separate coupling transfer from T3 C0 using known-load differences.

    In the chain, d(1/S21)/dC_add is independent of T3 L and C0. This checks
    whether changing C0/L could possibly fix the observed load dependence.
    """
    root=work.OUT/case/'known_Cadd';out=root/'derivative_diagnostic';out.mkdir(exist_ok=True)
    nets=[read_snp(root/f'cadd_{c:g}ff/response.snp') for c in CAPS_FF]
    obs=np.array([n.s[:,1,0] for n in nets]);f=nets[0].f
    intercept,slope=np.linalg.lstsq(np.c_[np.ones(4),CAPS_FF],1/obs,rcond=None)[0]
    old=np.array(work.REFERENCE_FIT['C_fF']+work.REFERENCE_FIT['G_fF'][1:])
    ls=[json.loads((work.OUT/case/'inductance.json').read_text())['L0_nH']]+work.REFERENCE_FIT['L_nH']
    ix=np.linspace(0,len(f)-1,1501).astype(int)
    weight=np.maximum(np.median(abs(obs),axis=0),.01)**2
    def values(freq,x):
        cs=[0.]+x[2:6].tolist();gs=x[:2].tolist()+x[6:9].tolist()
        zero=1/transmission(freq,ls,cs,gs,x[-1]);cs[0]=1.
        return zero,1/transmission(freq,ls,cs,gs,x[-1])-zero
    def residual(x):
        _,derivative=values(f[ix],x);d=400*(derivative-slope[ix])*weight[ix]
        return np.r_[d.real,d.imag]
    lo=np.r_[.1,.1,.99*old,-100];hi=np.r_[3000,3000,1.01*old,100]
    scale=np.r_[261,153,old,20]
    candidates=[]
    for gin in (100,261,500,1000):
        for gj in (80,153,300,1000):
            for delay in (-20,5,20):
                candidates.append(least_squares(residual,np.r_[gin,gj,old,delay],bounds=(lo,hi),
                    x_scale=scale,max_nfev=400))
    best=min(candidates,key=lambda q:q.fun@q.fun);zero,derivative=values(f,best.x)
    inferred=(intercept-zero)/derivative
    good=np.isfinite(inferred)&(weight>.0001)
    report=dict(method='Coupling-only fit to d(1/S21)/dC_add; T3 L/C0 cancel algebraically',
        Gin_fF=float(best.x[0]),Gjunction_fF=float(best.x[1]),
        older_relative_changes=(best.x[2:9]/old-1).tolist(),active_bounds=best.active_mask.tolist(),
        normalized_derivative_RMS=float(np.sqrt(np.mean(best.fun**2)*2)),
        inferred_C0_real_percentiles_fF=np.percentile(inferred.real[good],[5,50,95]).tolist(),
        inferred_C0_imag_percentiles_fF=np.percentile(inferred.imag[good],[5,50,95]).tolist(),
        scope='Diagnostic only; inferred C0 should be a real constant for an adequate chain model, not an accepted extraction.')
    (out/'result.json').write_text(json.dumps(report,indent=2)+'\n')
    np.savez(out/'diagnostic.npz',f=f,observed_slope=slope,model_slope=derivative,inferred_C0=inferred)
    print('DERIVATIVE CHECK',report,flush=True);return report


if __name__=='__main__':
    from geometry.taper.run_taper3 import CASE
    simulate_loads(CASE);fit_loads(CASE)
