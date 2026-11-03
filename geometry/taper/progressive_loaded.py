"""Known-capacitor perturbations constrain low-L T2 using one EM operator.

The numeric loads are controls, not physical capacitors in the final design.
Five direct fields include an independent held-out interpolation check.
"""
import json
from shutil import copytree
import xml.etree.ElementTree as ET
import numpy as np
from scipy.optimize import least_squares
from geometry.taper.progressive_fixtures import ROOT,NATIVE,OUT
from geometry.taper.sonnet import add_capacitor
from geometry.prepare_sonnet import tag,write_xml
from geometry.taper.response_cases import W
from geometry.taper.measure_compact_inductance import AP
from geometry.taper.progressive_fitting import transmission
from geometry.taper.staged_fitting import LT1
from geometry.taper.run_bridged_stage import digest
from circuit_common import read_snp
from geometry.taper.native_records import four_port_direct


def prepare():
    source=NATIVE/'03_T2_IDC_T1';dest=NATIVE/'03d_known_capacitor_operator'
    out=OUT/'03d_known_capacitor_operator';out.mkdir(parents=True,exist_ok=True)
    if not dest.exists():copytree(source,dest)
    project=dest/'generated_ncell_2.sonx';marker=out/'prepared.json'
    if not marker.exists():
        meta=json.loads((ROOT/'geometry/taper/output/progressive_extraction_2026-10-23/03_T2_IDC_T1/geometry.json').read_text())
        add_capacitor(project,meta)
        tree=ET.parse(project);comp=next(tree.iter(tag('Component')))
        comp.set('Type','PORTS');comp.attrib.pop('Units',None)
        for child in list(comp):
            if child.tag==tag('Value'):comp.remove(child)
        sweeps=tree.find(tag('Sweeps'))
        for child in list(sweeps):sweeps.remove(child)
        block=ET.SubElement(sweeps,tag('Set'),On='TRUE');freq=ET.SubElement(block,tag('Frequencies'))
        for f in (4,5,5.5,6,7):ET.SubElement(freq,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(f),Stop=str(f),Step='1')
        for q in tree.iter(tag('File')):
            if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','FALSE')
        write_xml(tree,project)
        marker.write_text(json.dumps(dict(source_sha256=digest(source/project.name),project_sha256=digest(project)),indent=2)+'\n')
    elif json.loads(marker.read_text())['project_sha256']!=digest(project):raise ValueError('Changed load operator')
    return project,out


def analyze(net,out):
    s,z0=AP['reduce_pairs'](net.s);y=AP['s_to_y'](s,z0)
    held=np.isclose(net.f,5.5e9,rtol=0,atol=1);train=~held
    basis=lambda f:(np.asarray(f)[:,None]/5e9)**np.array([-1,1,3,5])
    coeff=np.linalg.solve(basis(net.f[train]),y[train].reshape(4,-1))
    operator=lambda f:(basis(f)@coeff).reshape(len(f),3,3)
    err=float(np.max(abs(operator(net.f[held])-y[held])))
    if err>1e-5:raise ValueError(f'Held-out operator interpolation failed: {err} S')
    f=np.linspace(4e9,7e9,3001);yy=operator(f)
    loads=np.array([0,400,800,1200])
    actual=np.array([AP['add_capacitors'](yy,f,[float(c)])[:,1,0] for c in loads])
    baseline_net=read_snp(NATIVE/'03c_T2_IDC_T1_broad/generated_ncell_2.snp')
    right=np.minimum(np.searchsorted(baseline_net.f,f),len(baseline_net.f)-1)
    left=np.maximum(right-1,0)
    indices=np.where(abs(baseline_net.f[left]-f)<abs(baseline_net.f[right]-f),left,right)
    if not np.allclose(baseline_net.f[indices],f,rtol=0,atol=1):raise ValueError('No common zero-load grid')
    difference=abs(actual[0]-baseline_net.s[indices,1,0])
    contraction=float(np.max(abs(AP['add_capacitors'](y,net.f,[800])-W['terminate_capacitor'](net.s,net.f,.8))))
    check=dict(zero_load_S21_RMS=float(np.sqrt(np.mean(difference**2))),zero_load_S21_max=float(difference.max()),
               held_out_Y_error_S=err,capacitor_termination_identity_error=contraction)
    check['passed']=check['zero_load_S21_RMS']<.005 and check['zero_load_S21_max']<.05 and contraction<1e-7
    (out/'operator_validation.json').write_text(json.dumps(check,indent=2)+'\n')
    if not check['passed']:raise ValueError('Loaded operator failed its zero-load/termination validation')
    baseline=json.loads((OUT/'01_T1_regular/twoC/fit.json').read_text())
    ct=baseline['C_fF'][0];g=baseline['G_fF'][0];l=.7341365873386546
    rows={}
    import matplotlib.pyplot as plt
    for mode in ('two_C','free_junction','free_L_diagnostic'):
        seed=[100,459,5]+([g] if mode!='two_C' else [])+([l] if mode=='free_L_diagnostic' else [])
        lo=[0,30,-30]+([1] if mode!='two_C' else [])+([.1] if mode=='free_L_diagnostic' else [])
        hi=[1000,2000,40]+([180] if mode!='two_C' else [])+([2] if mode=='free_L_diagnostic' else [])
        scales=[100,459,5]+([57] if mode!='two_C' else [])+([.74] if mode=='free_L_diagnostic' else [])
        def model(x):
            gj=x[3] if mode!='two_C' else g;lx=x[4] if mode=='free_L_diagnostic' else l
            return np.array([transmission(f,[lx,LT1],[x[0]+c,ct],[x[1],gj],x[2]) for c in loads])
        def residual(x):
            d=(model(x)-actual).ravel();return np.r_[d.real,d.imag]
        fits=[least_squares(residual,[c,gin,5]+seed[3:],bounds=(lo,hi),x_scale=scales,max_nfev=500,
                            ftol=1e-11,xtol=1e-11,gtol=1e-11) for c,gin in ((100,459),(60,270),(300,600))]
        best=min(fits,key=lambda q:q.fun@q.fun);pred=model(best.x);delta=abs(pred-actual)
        r=dict(parameters=best.x.tolist(),parameter_names=['C0_T2_fF','Cg_input_T2_fF','delay_ps']+
               (['Cg_T2_T1_fF'] if mode!='two_C' else [])+(['L_T2_nH'] if mode=='free_L_diagnostic' else []),
               fixed_L_T2_nH=l if mode!='free_L_diagnostic' else None,fixed_C0_T1_fF=ct,
               fixed_Cg_T2_T1_fF=g if mode=='two_C' else None,
               complex_RMS=float(np.sqrt(np.mean(delta**2))),complex_max=float(delta.max()),
               per_load_RMS=np.sqrt(np.mean(delta**2,axis=1)).tolist(),active_bounds=best.active_mask.tolist(),
               success=bool(best.success),known_added_C_fF=loads.tolist(),held_out_Y_error_S=err,
               diagnostic_only=True,physical_design_accepted=False)
        rows[mode]=r;(out/f'{mode}.json').write_text(json.dumps(r,indent=2)+'\n')
        fig,axes=plt.subplots(2,2,figsize=(13,9))
        for ax,c,a,m in zip(axes.flat,loads,actual,pred):
            ax.plot(f/1e9,abs(a),'r-',lw=2,label='Loaded EM operator')
            ax.plot(f/1e9,abs(m),'k--',lw=2,label='Shared circuit fit')
            ax.set_title(f'Known added C = {c} fF');ax.set_xlabel('Frequency (GHz)');ax.set_ylabel('|S21|');ax.grid(alpha=.25);ax.legend()
        fig.suptitle(f'Known-C Controls: {mode}; RMS {r["complex_RMS"]:.5f}',fontsize=18)
        fig.tight_layout();fig.savefig(out/f'{mode}.png',dpi=170);plt.close(fig)
        print('LOADED',mode,r,flush=True)
    (out/'summary.json').write_text(json.dumps(rows,indent=2)+'\n')
    return rows


def main():
    project,out=prepare();marker=out/'completed.json'
    if marker.exists():
        r=json.loads(marker.read_text())
        if r['project_sha256']!=digest(project) or r['snp_sha256']!=digest(project.with_suffix('.snp')):raise ValueError('Changed loaded fields')
        net=W['read_network'](out/'operator.s4p',4)
    else:
        log=project.parent/'em.log'
        # A completed native run can have a mixed old-2-port/new-4-port export.
        # Recover its five complete direct matrices without repeating any fields.
        if not (log.exists() and 'Em simulation completed' in log.read_text() and 'Em errors detected: 0' in log.read_text()):
            try:W['em'](project,4)
            except ValueError:
                if not (log.exists() and 'Em simulation completed' in log.read_text() and 'Em errors detected: 0' in log.read_text()):raise
        net=four_port_direct(project.with_suffix('.snp'),[4,5,5.5,6,7])
        net.write_touchstone(str(out/'operator.s4p'))
        marker.write_text(json.dumps(dict(project_sha256=digest(project),snp_sha256=digest(project.with_suffix('.snp')),
            clean_operator_sha256=digest(out/'operator.s4p'),mixed_legacy_rows_excluded=True),indent=2)+'\n')
    return analyze(net,out)


if __name__=='__main__':main()
