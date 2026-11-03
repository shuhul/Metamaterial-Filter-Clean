"""Fixed-reference fit/assessment of the width-only no-IDC inductance sweep.

Primary reference is the complete-T1 calibration for the current topology. Other independently
obtained profiles are reported as sensitivity checks, never selected per width.
Only T2 bare C, T2-to-T1 Cg and labeled fixture delay are fitted.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import numpy as np
from scipy.optimize import least_squares
from scipy.signal import find_peaks
from circuit_common import read_snp
from geometry.taper.minimum_inductance import OUT,half_power_features
from geometry.taper.compact_stages import ROOT
from geometry.taper.staged_fitting import circuit,LT1,CT1,LR,CR,GOUT
from geometry.taper.locked_taper2 import FIXED
from geometry.taper.paired_coupling import abcd
from geometry.taper.run_bridged_stage import digest

joint=json.loads((ROOT/'output_data/taper_exploration/frank_coupling_controls_2026-10-23/bridged_T1_existing/joint_port_fit/fit.json').read_text())
full_t1=json.loads((ROOT/'output_data/taper_exploration/frank_coupling_controls_2026-10-23/full_T1_free_junction/N4_independent/fit.json').read_text())
PROFILES={
    'selected_T1_primary':dict(ct=CT1,cr=FIXED['C_regular_fF'],ce=FIXED['C_regular_fF'],
                              gj=FIXED['Cg_junction_fF'],go=FIXED['Cg_output_fF']),
    'staged_N4':dict(ct=CT1,cr=CR,ce=CR,gj=8.40642914296533,go=GOUT),
    'bridged_joint':dict(ct=CT1,cr=joint['C0_regular_fF'],ce=CR,
                        gj=joint['Cg_junction_fF'],go=joint['Cg_output_fF']),
    'full_T1_staged_calibration':dict(ct=full_t1['C0_T1_fF'],cr=CR,ce=CR,
                                     gj=full_t1['Cg_junction_fF'],go=GOUT)}


def predictor(f,l,profile):
    p=PROFILES[profile];w=2*np.pi*f
    n=abcd(circuit(f,[LT1,LR,LR],[p['ct'],p['cr'],p['ce']],[None,p['gj'],p['go']]))
    def predict(x):
        c,g,t=x;z=1/(1j*w*g*1e-15)
        a=n[:,0,0]+z*n[:,1,0];b=n[:,0,1]+z*n[:,1,1]
        y=1j*w*c*1e-15+1/(1j*w*l*1e-9)
        den=(a+b/50)*(1+50*y)+50*n[:,1,0]+n[:,1,1]
        return 2/den*np.exp(-1j*w*t*1e-12)
    return predict


def assessment(f,actual,model,r):
    a=half_power_features(f,actual);b=half_power_features(f,model)
    errors=[];paired=len(a)==len(b) and len(a)>=2
    if paired:
        for p,q in zip(a,b):
            # A merged absolute half-power region has no individual linewidth.
            # Do not invent one or reject an otherwise matching curve for that.
            width=None if p['power_FWHM_MHz'] is None or q['power_FWHM_MHz'] is None else abs(p['power_FWHM_MHz']-q['power_FWHM_MHz'])/p['power_FWHM_MHz']
            frequency=abs(p['GHz']-q['GHz'])*1000
            tolerance=.1*p['power_FWHM_MHz'] if p['power_FWHM_MHz'] is not None else 2.
            errors.append(dict(height_relative=abs(p['magnitude']-q['magnitude'])/p['magnitude'],
                               width_relative=width,center_error_MHz=frequency,center_tolerance_MHz=tolerance))
    passed=(r['complex_RMS']<.01 and r['complex_max']<.05 and paired and not any(r['active_bounds'])
            and all(e['height_relative']<.1 and (e['width_relative'] is None or e['width_relative']<.1)
                    and e['center_error_MHz']<=e['center_tolerance_MHz'] for e in errors))
    return dict(pass_curve=bool(passed),observed_peaks=a,model_peaks=b,peak_errors=errors,
                version=2,merged_width_policy='Individual half-power width omitted; use curve errors, heights and 2 MHz center limit',
                scope='4.65-5.4 GHz no-IDC fixed-reference curve agreement, not final T2 target acceptance')


def fit(snp,l,out,profile,data_label='Sonnet simulation'):
    out.mkdir(parents=True,exist_ok=True);n=read_snp(snp)
    marker=out/'fit.json';provenance=dict(source_sha256=digest(snp),L_nH=l,profile=profile,constants=PROFILES[profile],
                                       code_sha256=digest(__file__))
    if data_label!='Sonnet simulation':provenance['data_label']=data_label
    if marker.exists():
        old=json.loads(marker.read_text())
        if old['provenance']==provenance:return old
    f=n.f;actual=n.s[:,1,0];ix=np.unique(np.linspace(0,len(f)-1,1501).astype(int))
    predict=predictor(f[ix],l,profile)
    def residual(x,smooth=False):
        m=predict(x);obs=actual[ix]
        d=(1/m-1/obs)*np.maximum(abs(obs),.03)**2 if smooth else m-obs
        return np.r_[d.real,d.imag]
    bounds=([0,5,-20],[2000,180,30]);seeds=[]
    for c in (0,100,300,800):
        for g in (40,56.7,75):
            q=least_squares(lambda x:residual(x,True),[c,g,5],bounds=bounds,x_scale=[100,57,5],max_nfev=150)
            seeds.append(q.x)
    fits=[least_squares(residual,x,bounds=bounds,x_scale=[100,57,5],max_nfev=300,
                        ftol=1e-11,xtol=1e-11,gtol=1e-11)
          for x in sorted(seeds,key=lambda x:np.linalg.norm(residual(x)))[:4]]
    best=min(fits,key=lambda q:q.fun@q.fun);full=predictor(f,l,profile)
    def full_residual(x):
        d=full(x)-actual;return np.r_[d.real,d.imag]
    best=least_squares(full_residual,best.x,bounds=bounds,x_scale=[100,57,5],max_nfev=200)
    model=full(best.x);error=abs(model-actual)
    r=dict(provenance=provenance,C_bare_fF=float(best.x[0]),Cg_T2_T1_fF=float(best.x[1]),delay_ps=float(best.x[2]),
           complex_RMS=float(np.sqrt(np.mean(error**2))),complex_max=float(error.max()),
           magnitude_RMS=float(np.sqrt(np.mean((abs(model)-abs(actual))**2))),active_bounds=best.active_mask.tolist(),
           T1_parameters_refitted=False,free_parameters=['C_bare_T2','Cg_T2_T1','fixture_delay'])
    r['assessment']=assessment(f,actual,model,r)
    marker.write_text(json.dumps(r,indent=2)+'\n');np.savez(out/'response.npz',f=f,actual=actual,model=model)
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(12,7))
    ax.plot(f/1e9,abs(actual),color='#b22222',lw=2.8,label=data_label)
    ax.plot(f/1e9,abs(model),'k--',lw=2.2,label='Fixed downstream circuit + T2 fit')
    text=(f"Measured L = {l:.4f} nH\nBare C = {r['C_bare_fF']:.3f} fF\n"
          f"T2-to-T1 Cg = {r['Cg_T2_T1_fF']:.3f} fF\nComplex S21 RMS = {r['complex_RMS']:.4f}\n"
          f"Curve check: {'PASS' if r['assessment']['pass_curve'] else 'FAIL'}\nReference: {profile}")
    ax.text(.03,.96,text,transform=ax.transAxes,va='top',fontsize=14,bbox=dict(boxstyle='round,pad=.5',fc='white',ec='gray'))
    ax.set_ylim(0,max(abs(actual).max(),abs(model).max())*1.7)
    ax.set_title('Minimum-Inductance Test: No Input IDC',fontsize=23,pad=14)
    ax.set_xlabel('Frequency (GHz)',fontsize=19);ax.set_ylabel(r'$|S_{21}|$',fontsize=23)
    ax.tick_params(labelsize=16);ax.legend(fontsize=12,loc='upper right');ax.grid(alpha=.25)
    fig.tight_layout();fig.savefig(out/'s21_fit.png',dpi=190);plt.close(fig)
    print(profile,l,r['complex_RMS'],r['complex_max'],r['assessment']['pass_curve'],flush=True)
    return r


def review_width(width):
    case=OUT/f'meander_w{width}';cal=json.loads((case/'Lcal/results.json').read_text());l=cal['L0_nH']
    source=ROOT/f'sonnet_files/taper_exploration/minimum_inductance_2026-10-23/meander_w{width}/direct/generated_ncell_3.snp'
    status=json.loads((case/'direct/run_status.json').read_text())
    if status['status']!='complete' or status['snp_sha256']!=digest(source):raise ValueError('Incomplete response')
    rows={profile:fit(source,l,case/'review'/profile,profile) for profile in PROFILES}
    from geometry.taper.minimum_inductance_paired import compare
    paired=compare(width,l,rows)
    r=dict(width_um=width,L_nH=l,primary_profile='full_T1_staged_calibration',
           primary_pass=rows['full_T1_staged_calibration']['assessment']['pass_curve'],
           all_profiles_pass=all(v['assessment']['pass_curve'] for v in rows.values()),profiles=rows,
           paired_ideal_control=paired)
    (case/'review/summary.json').write_text(json.dumps(r,indent=2)+'\n')
    return r


if __name__=='__main__':
    import sys
    review_width(int(sys.argv[1]))
