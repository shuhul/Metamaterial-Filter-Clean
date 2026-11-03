"""Native known-CAP cross-check of a directly fed T3; production stays unchanged.

All inductances fixed. Shared C0/Gjunction across loads. Older Cs either fixed
or bounded once to1% of the selected T2 reference, never cumulatively changed.
One four-thread EM job; later capacitances use Sonnet's existing cache only.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import shutil
import numpy as np
from scipy.optimize import least_squares
from geometry.taper import taper3 as work
from geometry.taper.sonnet import add_capacitor, configure, simulate
from geometry.taper.progressive_fitting import transmission
from circuit_common import read_snp, write_sonx

CAPS_FF=[.1,100.,200.]


def simulate_loads(case):
    source=work.NATIVE/case/'direct/generated_ncell_4.sonx'
    meta=json.loads((source.parent/'geometry.json').read_text())
    native=work.NATIVE/case/'direct_known_Cadd'
    out=work.OUT/case/'direct_known_Cadd'
    native.mkdir(exist_ok=True);out.mkdir(exist_ok=True)
    project=native/source.name
    expected=dict(source_project=str(source),source_sha256=work.digest(source),
        C_add_fF=CAPS_FF,method='Native explicit CAP; one EM then cache-only changes')
    manifest=out/'manifest.json'
    if manifest.exists() and json.loads(manifest.read_text())!=expected:
        raise ValueError('Known-load source/settings changed')
    manifest.write_text(json.dumps(expected,indent=2)+'\n')
    if not project.exists():
        shutil.copy2(source,project)
        add_capacitor(project,meta)
        configure(project,4.65,5.4,7500)
    for i,cap in enumerate(CAPS_FF):
        folder=out/f'cadd_{cap:g}ff';folder.mkdir(exist_ok=True)
        marker=folder/'completed.json'
        if marker.exists():
            record=json.loads(marker.read_text())
            for name in ('project.sonx','response.snp'):
                if work.digest(folder/name)!=record['files_sha256'][name]:
                    raise ValueError('Archived load changed')
            continue
        write_sonx({'C_add':str(cap/1000)},project)
        print('DIRECT KNOWN LOAD',cap,'fF',flush=True)
        simulate(project,folder,cache_only=i>0)
        marker.write_text(json.dumps(dict(C_add_fF=cap,files_sha256={
            n:work.digest(folder/n) for n in ('project.sonx','response.snp')}),indent=2)+'\n')
    return out


def fit_loads(case,older_fraction=0.):
    if older_fraction not in (0.,.01):raise ValueError('Only fixed or one1% reference adjustment')
    root=work.OUT/case/'direct_known_Cadd'
    out=root/('fit_fixed_older' if not older_fraction else 'fit_older_1percent')
    out.mkdir(exist_ok=True)
    refs=work.REFERENCE_FIT
    ls=[json.loads((work.OUT/case/'inductance.json').read_text())['L0_nH']]+refs['L_nH']
    old=np.array(refs['C_fF']+refs['G_fF'][1:])
    paths=[root/f'cadd_{c:g}ff/response.snp' for c in CAPS_FF]
    nets=[read_snp(p) for p in paths]
    indices=[np.unique(np.linspace(0,len(n.f)-1,1501).astype(int)) for n in nets]
    lo=np.array([0.,.1,-100.]);hi=np.array([1000.,3000.,100.]);scale=np.array([251.,153.,20.])
    if older_fraction:
        lo=np.r_[lo,old*(1-older_fraction)];hi=np.r_[hi,old*(1+older_fraction)];scale=np.r_[scale,old]
    def expand(x,cap):
        older=x[3:] if older_fraction else old
        return [x[0]+cap]+older[:4].tolist(),[None,x[1]]+older[4:].tolist()
    def residual(x,full=False):
        result=[]
        for n,ix,cap in zip(nets,indices,CAPS_FF):
            j=np.arange(len(n.f)) if full else ix
            c,g=expand(x,cap)
            d=transmission(n.f[j],ls,c,g,x[2])-n.s[j,1,0]
            result.extend([d.real,d.imag])
        return np.concatenate(result)
    candidates=[]
    for c in (50.,150.,251.,450.,800.):
        for g in (80.,153.,300.,800.):
            x=np.array([c,g,0.])
            if older_fraction:x=np.r_[x,old]
            n=nets[0];peak=int(np.argmax(abs(n.s[:,1,0])))
            cs,gs=expand(x,CAPS_FF[0])
            zero=transmission(n.f[peak:peak+1],ls,cs,gs,0)[0]
            x[2]=np.clip(-np.angle(n.s[peak,1,0]/zero)/(2*np.pi*n.f[peak])*1e12,-99,99)
            candidates.append(least_squares(residual,x,bounds=(lo,hi),x_scale=scale,max_nfev=600))
    best=min(candidates,key=lambda q:q.fun@q.fun)
    best=least_squares(lambda x:residual(x,True),best.x,bounds=(lo,hi),x_scale=scale,
        max_nfev=800,ftol=1e-11,xtol=1e-11,gtol=1e-11)
    cs,gs=expand(best.x,0)
    report=dict(case=case,accepted=False,method='Shared native known-CAP direct-feed fit',
        L_nH=ls,C_fF=cs,G_fF=gs,delay_ps=float(best.x[2]),older_fractional_bounds=older_fraction,
        active_bounds=best.active_mask.tolist(),sources={str(p):work.digest(p) for p in paths},loads=[])
    if older_fraction:report['older_fractional_changes']=(best.x[3:]/old-1).tolist()
    jac=best.jac*scale[None,:];sv=np.linalg.svd(jac,compute_uv=False)
    inv=np.linalg.pinv(jac.T@jac);std=np.sqrt(np.maximum(np.diag(inv),1e-30))
    report['local_sensitivity']=dict(singular_values=sv.tolist(),condition_number=float(sv[0]/sv[-1]),
        correlation=(inv/(std[:,None]*std[None,:])).tolist(),
        scope='Local numerical diagnostic, not a confidence interval or proof of uniqueness')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,3,figsize=(19,6))
    for ax,n,cap in zip(axes,nets,CAPS_FF):
        cs,gs=expand(best.x,cap);model=transmission(n.f,ls,cs,gs,best.x[2])
        error=abs(model-n.s[:,1,0]);rms=float(np.sqrt(np.mean(error**2)))
        report['loads'].append(dict(C_add_fF=cap,complex_RMS=rms,complex_max=float(error.max())))
        np.savez(out/f'cadd_{cap:g}ff.npz',f=n.f,actual=n.s,model=model)
        ax.plot(n.f/1e9,abs(n.s[:,1,0]),lw=2,label='Native Sonnet CAP')
        ax.plot(n.f/1e9,abs(model),'k--',lw=2,label='Shared fit')
        ax.set(xlabel='Frequency (GHz)',ylabel=r'$|S_{21}|$',title=f'Added C={cap:g} fF; RMS={rms:.4f}')
        ax.grid(alpha=.2);ax.legend()
    fig.suptitle(f'Direct T3 Known-Capacitor Check: C0={best.x[0]:.2f} fF, Cg={best.x[1]:.2f} fF',fontsize=18)
    fig.tight_layout();fig.savefig(out/'shared_fit.png',dpi=170);plt.close(fig)
    unloaded=read_snp(work.NATIVE/case/'direct/generated_ncell_4.snp')
    report['unloaded_vs_native_0p1fF_complex_RMS']=float(np.sqrt(np.mean(abs(unloaded.s[:,1,0]-nets[0].s[:,1,0])**2)))
    (out/'fit.json').write_text(json.dumps(report,indent=2)+'\n')
    print('DIRECT LOAD FIT',json.dumps(report),flush=True)
    return report
