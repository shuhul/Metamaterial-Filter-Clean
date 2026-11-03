"""Joint complex-S21 extraction from matching bare and IDC-equipped fixtures.

All inductances and regular values remain fixed. The unchanged T1 electrode
geometry supplies three shared capacitor parameters, which are tested against
both responses rather than pinned to arbitrary earlier numbers.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import numpy as np
from scipy.optimize import least_squares
from geometry.taper.clearance_extraction import OUT
from geometry.taper.corrected_iteration import REGULAR_C,REGULAR_G
from geometry.taper.progressive_fitting import transmission
from geometry.taper.run_bridged_stage import digest
from geometry.specifications import S21_COMPLEX_RMS_CUTOFF
from circuit_common import read_snp
from geometry.taper.staged_fitting import circuit


def run(case='frame700_wire88_left214',*,profiles=False,known_loads=False,reflections=False,window_GHz=None):
    root=OUT/case
    name='staged_joint_with_reflections' if reflections else 'staged_joint_with_loads' if known_loads else 'staged_joint'
    if window_GHz is not None:name+='_window_'+str(window_GHz[0]).replace('.','p')+'_'+str(window_GHz[1]).replace('.','p')
    out=root/name;out.mkdir(exist_ok=True)
    records=[json.loads((root/s/'fit/fit.json').read_text()) for s in ('direct','idc')]
    arrays=[]
    for s in ('direct','idc'):
        with np.load(root/s/'fit/response.npz') as data:
            arrays.append({k:data[k] for k in data.files})
    for r,a in zip(records,arrays):
        source=r['provenance']['source']
        if digest(source)!=r['provenance']['sha256']:raise ValueError('Fit source changed')
        original=read_snp(source)
        if (not np.allclose(original.f,a['f'],rtol=0,atol=1)
            or not np.allclose(original.s,a['actual'],rtol=1e-11,atol=1e-11)):
            raise ValueError('Fit snapshot differs from its original source')
    original_arrays=[dict(f=a['f'].copy(),actual=a['actual'].copy()) for a in arrays]
    if window_GHz is not None:
        for a in arrays:
            mask=(a['f']>=window_GHz[0]*1e9)&(a['f']<=window_GHz[1]*1e9)
            a['f']=a['f'][mask];a['actual']=a['actual'][mask]
    ls=[r['L_nH'] for r in records]
    additions=[0,0];names=['bare','idc'];labels=['No IDC','Uniform IDC; 214 µm Left Clearance']
    loaded_sources=[];heldout=None
    if known_loads:
        control=root/(known_loads if isinstance(known_loads,str) else 'known_capacitor_control')
        if not json.loads((control/'open_validation.json').read_text())['passed']:
            raise ValueError('Known-load reference did not pass open validation')
        loads=json.loads((control/'loads.json').read_text())
        for cap in (400,1200,800):
            source=next(q for q in loads if q['known_C_fF']==cap)
            if digest(source['path'])!=source['sha256']:raise ValueError('Known-load data changed')
            net=read_snp(source['path'])
            if cap==800:
                heldout=(net,source)
                continue
            arrays.append(dict(f=net.f,actual=net.s));ls.append(list(ls[1]))
            additions.append(cap);names.append(f'idc_plus_{cap}fF');labels.append(f'IDC + Known {cap} fF')
            loaded_sources.append(source)
    freq=[a['f'] for a in arrays]
    observed=[a['actual'][:,1,0] for a in arrays]
    indices=[np.unique(np.linspace(0,len(f)-1,min(2001,len(f))).astype(int)) for f in freq]
    all_indices=[np.arange(len(f)) for f in freq]
    scale=np.array([100,459,100,56.7,8.2865,20,20.])
    lo=np.array([0,.1,0,.1,.1,-100,-100.])/scale
    hi=np.array([1000,3000,1000,3000,3000,100,100.])/scale
    if reflections:
        scale=np.r_[scale,20.];lo=np.r_[lo,-5.];hi=np.r_[hi,5.]
    def network(x,j,ix):
        v=x*scale
        return circuit(freq[j][ix],ls[j],[v[0]+additions[j] if j else 0,v[2],REGULAR_C,REGULAR_C],
             [v[1] if j else None,v[3],v[4],REGULAR_G],left_ps=v[5 if j==0 else 6]-v[7],right_ps=v[7])
    def predict(x,j,ix):
        v=x*scale
        return transmission(freq[j][ix],ls[j],[v[0]+additions[j] if j else 0,v[2],REGULAR_C,REGULAR_C],
             [v[1] if j else None,v[3],v[4],REGULAR_G],v[5 if j==0 else 6])
    def residual(x,ixs=indices,smooth=False):
        deltas=[]
        for j,ix in enumerate(ixs):
            s=predict(x,j,ix);a=observed[j][ix]
            d=(1/s-1/a)*np.maximum(abs(a),.03)**2 if smooth else s-a
            deltas.extend([d.real,d.imag])
            if reflections and not smooth:
                model=network(x,j,ix)
                for port in (0,1):
                    d=model[:,port,port]-arrays[j]['actual'][ix,port,port]
                    deltas.extend([d.real,d.imag])
        return np.concatenate(deltas)
    bare,idc=records
    seeds=[]
    for c,g in [(max(idc['C_fF'][0],1e-5),idc['G_fF'][0]),(100,459),(300,800),(600,1800)]:
        for shared in (bare,idc):
            v=[c,g,shared['C_fF'][1],shared['G_fF'][1],shared['G_fF'][2],bare['delay_ps'],idc['delay_ps']]
            if reflections:v.append(.74)
            seeds.append(np.clip(np.array(v)/scale,lo+1e-10,hi-1e-10))
    smooth=[least_squares(lambda x:residual(x,smooth=True),x,bounds=(lo,hi),max_nfev=400) for x in seeds]
    starts=sorted(seeds+[q.x for q in smooth],key=lambda x:np.linalg.norm(residual(x)))[:5]
    fits=[least_squares(residual,x,bounds=(lo,hi),max_nfev=400,
              ftol=1e-11,xtol=1e-11,gtol=1e-11) for x in starts]
    best=min(fits,key=lambda q:q.fun@q.fun)
    best=least_squares(lambda x:residual(x,all_indices),best.x,bounds=(lo,hi),max_nfev=400,
                       ftol=1e-11,xtol=1e-11,gtol=1e-11)
    v=best.x*scale
    # Remove both delay directions before testing the five capacitor directions.
    jac=best.jac/np.sqrt(sum(map(len,freq)))
    q,_=np.linalg.qr(jac[:,5:]);cj=jac[:,:5]-q@(q.T@jac[:,:5])
    _,sv,directions=np.linalg.svd(cj,full_matrices=False)
    errors=[abs(predict(best.x,j,all_indices[j])-observed[j]) for j in range(len(freq))]
    result=dict(fitted=dict(C0_T2_fF=v[0],Cg_input_T2_fF=v[1],C0_T1_fF=v[2],
                    Cg_T2_T1_fF=v[3],Cg_T1_regular_fF=v[4],delay_bare_ps=v[5],delay_IDC_ps=v[6]),
        fixed_inductances_nH=ls,fixed_regular_C_fF=REGULAR_C,fixed_output_Cg_fF=REGULAR_G,
        complex_RMS=[float(np.sqrt(np.mean(e**2))) for e in errors],
        complex_max=[float(e.max()) for e in errors],active_bounds=best.active_mask.tolist(),
        optimizer_success=bool(best.success),curve_RMS_cutoff=S21_COMPLEX_RMS_CUTOFF,
        capacitor_singular_values_after_delay_projection=sv.tolist(),
        capacitor_condition=float(sv[0]/sv[-1]),weakest_capacitor_direction=directions[-1].tolist(),
        parameter_scales=scale.tolist(),
        source_fit_sha256=[digest(root/s/'fit/fit.json') for s in ('direct','idc')],
        source_data_sha256=[r['provenance']['sha256'] for r in records],
        stage_labels=labels,known_added_capacitances_fF=additions,known_load_sources=loaded_sources,
        scope='Shared T1 capacitors tested jointly against bare and IDC responses. Fixed L and regular values; unweighted real/imaginary S21 residuals. No target optimization.')
    if reflections:
        result['scope']='Shared capacitors fitted to complex S11, S21 and S22 of both stages. Fixed L/regular values; three fixture delays (two total delays and one common output delay). S21 acceptance remains RMS <= 0.05.'
        result['common_right_delay_ps']=float(v[7])
        result['reflection_complex_RMS']=[{f'S{p+1}{p+1}':float(np.sqrt(np.mean(abs(network(best.x,j,all_indices[j])[:,p,p]-arrays[j]['actual'][:,p,p])**2))) for p in (0,1)} for j in range(len(freq))]
    result['curves_pass']=bool(best.success and all(e<=S21_COMPLEX_RMS_CUTOFF for e in result['complex_RMS']))
    result['fit_window_GHz']=window_GHz
    if window_GHz is not None:
        result['frequency_holdout']=[]
        for j,a in enumerate(original_arrays):
            prediction=circuit(a['f'],ls[j],[v[0] if j else 0,v[2],REGULAR_C,REGULAR_C],
                [v[1] if j else None,v[3],v[4],REGULAR_G],left_ps=v[5+j]-(v[7] if reflections else 0),
                right_ps=v[7] if reflections else 0)
            mask=(a['f']<window_GHz[0]*1e9)|(a['f']>window_GHz[1]*1e9)
            checks={f'S{p+1}{q+1}':float(np.sqrt(np.mean(abs(prediction[mask,p,q]-a['actual'][mask,p,q])**2))) for p,q in ((0,0),(1,0),(1,1))}
            result['frequency_holdout'].append(checks)
    result['interior_capacitors']=not any(best.active_mask[:5])
    result['extraction_accepted']=False  # Requires profiles/independent validation, not just a rank test.
    if heldout is not None:
        net,source=heldout
        prediction=transmission(net.f,ls[1],[v[0]+800,v[2],REGULAR_C,REGULAR_C],
            [v[1],v[3],v[4],REGULAR_G],v[6])
        err=abs(prediction-net.s[:,1,0])
        result['held_out_800fF']=dict(complex_RMS=float(np.sqrt(np.mean(err**2))),complex_max=float(err.max()),
            curve_pass=bool(np.sqrt(np.mean(err**2))<=S21_COMPLEX_RMS_CUTOFF),source=source,
            scope='Different termination of the same validated EM operator; not a separate EM simulation.')
        np.savez(out/'held_out_800fF.npz',f=net.f,actual=net.s[:,1,0],model=prediction)
    result['multistart_candidates']=[dict(values=(q.x*scale).tolist(),
        pooled_complex_RMS=float(np.linalg.norm(residual(q.x,all_indices))/np.sqrt(sum(map(len,freq))))) for q in fits]
    if profiles:
        result['profiles']=[]
        for parameter in (0,1):
            others=[j for j in range(len(scale)) if j!=parameter]
            center=float(v[parameter])
            branches=([50,100,200,400,800],) if center<10 else ([center*.9,center*.75,center*.5],[center*1.1,center*1.25,center*1.5])
            rows=[]
            for targets in branches:
                previous=best.x.copy();start=center
                for target in targets:
                    target=float(np.clip(target,lo[parameter]*scale[parameter]+1e-9,hi[parameter]*scale[parameter]-1e-9))
                    steps=max(1,int(np.ceil(abs(target-start)/max(center*.01,5))))
                    for fixed_value in np.linspace(start,target,steps+1)[1:]:
                        def expand(q):
                            x=previous.copy();x[parameter]=fixed_value/scale[parameter];x[others]=q
                            return x
                        q=least_squares(lambda y:residual(expand(y)),previous[others],
                            bounds=(lo[others],hi[others]),max_nfev=150,
                            ftol=1e-10,xtol=1e-10,gtol=1e-10)
                        previous=expand(q.x)
                    start=target
                    per_stage=[float(np.sqrt(np.mean(abs(predict(previous,j,all_indices[j])-observed[j])**2))) for j in range(len(freq))]
                    rows.append(dict(fixed_value_fF=target,all_values=(previous*scale).tolist(),stage_RMS=per_stage,
                        pooled_RMS=float(np.linalg.norm(residual(previous,all_indices))/np.sqrt(sum(map(len,freq)))),converged=bool(q.success)))
            result['profiles'].append(dict(parameter='C0_T2' if parameter==0 else 'Cg_input_T2',optimum_fF=center,rows=rows))
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    for j,name in enumerate(names):
        np.savez(out/(name+'_response.npz'),f=freq[j],actual=observed[j],model=predict(best.x,j,all_indices[j]))
    import matplotlib.pyplot as plt
    rows=int(np.ceil(len(freq)/2))
    fig,axs=plt.subplots(rows,2,figsize=(17,7*rows),sharey=True);axs=np.ravel(axs)
    for j,ax in enumerate(axs):
        ax.plot(freq[j]/1e9,abs(observed[j]),color='#b22222',lw=2.7,label='Sonnet EM')
        ax.plot(freq[j]/1e9,abs(predict(best.x,j,all_indices[j])),'k--',lw=2.1,
                label='Complex S-parameter fit' if reflections else 'Joint complex S21 fit')
        ax.set(title=labels[j],xlabel='Frequency (GHz)',ylim=(0,1.3))
        ax.title.set_size(19);ax.xaxis.label.set_size(18);ax.tick_params(labelsize=14);ax.grid(alpha=.22);ax.legend(fontsize=12,loc='upper right')
        text=f"Complex RMS = {result['complex_RMS'][j]:.5f}\nFixed L = {ls[j][0]:.5f} nH\nShared T1 C0 = {v[2]:.2f} fF\nShared T2–T1 Cg = {v[3]:.2f} fF\nShared T1–unit Cg = {v[4]:.3f} fF"
        if j:text+=f"\nT2 C0 = {v[0]:.2f} fF\nInput Cg = {v[1]:.2f} fF"
        ax.text(.025,.97,text,transform=ax.transAxes,va='top',fontsize=13,bbox=dict(boxstyle='round,pad=.5',fc='white',ec='gray'))
    axs[0].set_ylabel(r'$|S_{21}|$',fontsize=23)
    fig.suptitle('Staged Capacitance Extraction: All Inductances and Regular Values Fixed',fontsize=21)
    fig.tight_layout();fig.savefig(out/'s21_staged_joint.png',dpi=170);plt.close(fig)
    print('STAGED JOINT',json.dumps(dict(fitted=result['fitted'],RMS=result['complex_RMS'],
          active_bounds=result['active_bounds'],capacitor_condition=result['capacitor_condition'])),flush=True)
    return result


if __name__=='__main__':run()
