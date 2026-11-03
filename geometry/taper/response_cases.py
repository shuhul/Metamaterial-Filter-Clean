"""Prepare and solve band-limited six-port T2 response-screening cases."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json,runpy,shutil
import xml.etree.ElementTree as ET
from functools import partial
import numpy as np
from geometry.prepare_sonnet import main as prepare,tag,write_xml
from geometry.taper.locked_taper2 import prepare_case,locked_validate,require_native_locked
from geometry.taper.second_layout import from_dict
from geometry.taper.optimize_response import ROOT,OUT,A,theory,THEORY_VERSION

G=runpy.run_path(str(ROOT/'investigations/taper2_fit_2026-10-22/gap_probes.py'))
W=G['W']

def case(name,p,*,frequencies=(4.5,4.9,5.,5.5),analyze_screen=True,
         prepared_gds=None,geometry_validator=None,native_validator=None):
    if prepared_gds is None:
        gds,out=prepare_case(name,p,allow_t2_junction=True)
        validator=partial(locked_validate,allow_t2_junction=True)
    else:
        gds=Path(prepared_gds);out=ROOT/'output_data/taper_exploration'/name
        out.mkdir(parents=True,exist_ok=True);validator=geometry_validator
        if validator is None:raise ValueError('Custom geometry requires its validator')
        validator(gds,json.loads(gds.with_suffix('.json').read_text()))
    native=ROOT/'sonnet_files/taper_exploration'/name
    project=native/'generated_ncell_3.sonx'
    if not project.exists():prepare(3,gds_path=gds,output_dir=native,geometry_validator=validator)
    meta=json.loads(gds.with_suffix('.json').read_text())
    (native_validator or require_native_locked)(project,meta)
    (out/'response_intent.json').write_text(json.dumps({'goal':'Direct theoretical S21 agreement',
        'T1_geometry_locked':True,'T2_junction_is_design_variable':True,
        'screen_frequencies_GHz':list(frequencies),'bare_inductance_not_resimulated_for_this_screen':True},indent=2)+'\n')
    # Reuse the reviewed constructor in a separate namespace and destination.
    fn=G['prepare'];fn.__globals__['OUT']=out;fn.__globals__['CASE']=name
    probe=fn() # returns out/series_gap_probes/series_probe.sonx
    marker=probe.parent/'screen_completed.json'
    manifest=json.loads((probe.parent/'manifest.json').read_text())
    if W['digest'](project)!=manifest['source_sha256']:
        raise ValueError('Native source changed after probe construction; revalidate using a new case')
    if marker.exists():
        saved=json.loads(marker.read_text())
        if W['digest'](probe)!=saved['project_sha256'] or W['digest'](probe.with_suffix('.snp'))!=saved['snp_sha256']:
            raise ValueError('Completed screening artifact changed; refusing stale cache reuse')
        completed=saved['frequencies_GHz']
    else:completed=[]
    if not set(frequencies).issubset(completed):
        tree=ET.parse(probe);freq=tree.find(f"{tag('Sweeps')}/{tag('Set')}/{tag('Frequencies')}")
        for q in list(freq):freq.remove(q)
        for f in frequencies:
            ET.SubElement(freq,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(f),Stop=str(f),Step='1')
        write_xml(tree,probe);n=W['em'](probe,6)
        marker.write_text(json.dumps({'nports':6,'frequencies_GHz':(n.f/1e9).tolist(),
            'project_sha256':W['digest'](probe),'snp_sha256':W['digest'](probe.with_suffix('.snp'))},indent=2)+'\n')
    return analyze(name) if analyze_screen else W['read_network'](probe.with_suffix('.snp'),6)

def analyze(name):
    out=ROOT/'output_data/taper_exploration'/name
    n=W['read_network'](out/'series_gap_probes/series_probe.snp',6)
    p=np.diag([1,1,-1,-1,-1,-1]);z=A['s_to_z'](p@n.s@p)
    r=n.f/5e9;powers=np.arange(-1,2*len(r)-2,2)
    coef=np.linalg.solve(r[:,None]**powers,z.reshape(len(r),-1))
    f=np.linspace(4.5e9,5.5e9,2001)
    zz=(((f/5e9)[:,None]**powers)@coef).reshape(len(f),6,6)
    s=A['short_internal_z'](zz);target=theory(f)
    from scipy.signal import find_peaks,peak_widths
    from scipy.optimize import minimize_scalar
    def features(v):
        ix,_=find_peaks(abs(v),prominence=.015)
        widths=peak_widths(abs(v)**2,ix,rel_height=.5)[0]*(f[1]-f[0])/1e6
        return [{'GHz':float(f[i]/1e9),'height':float(abs(v[i])),
                 'FWHM_MHz':float(width)} for i,width in zip(ix,widths)]
    d=s[:,1,0]-target[:,1,0]
    def delay_cost(tau):return np.mean(abs(s[:,1,0]-target[:,1,0]*np.exp(-2j*np.pi*f*tau*1e-12))**2)
    delay=minimize_scalar(delay_cost,bounds=(-30,30),method='bounded')
    result={'case':name,'theory_version':THEORY_VERSION,'screen_prediction_only':True,'physical_acceptance':False,
            'S21_complex_RMSE':float(np.sqrt(np.mean(abs(d)**2))),
            'S21_magnitude_RMSE':float(np.sqrt(np.mean((abs(s[:,1,0])-abs(target[:,1,0]))**2))),
            'extra_reference_delay_ps':float(delay.x),'delay_aligned_S21_complex_RMSE':float(np.sqrt(delay.fun)),
            'peaks':features(s[:,1,0]),'target_peaks':features(target[:,1,0]),
            'max_singular_value':float(np.linalg.svd(s,compute_uv=False).max())}
    np.savez(out/'response_screen.npz',f=f,s=s,target=target,coefficients=coef,powers=powers)
    (out/'response_screen.json').write_text(json.dumps(result,indent=2)+'\n')
    plt=A['plt'];fig,ax=plt.subplots(figsize=(11,6))
    ax.plot(f/1e9,abs(s[:,1,0]),'r',lw=2,label='EM multiport prediction')
    ax.plot(f/1e9,abs(target[:,1,0]),'k--',lw=2,label='Specified theoretical filter')
    ax.set(xlabel='Frequency (GHz)',ylabel='|S21|',title=f'T2 response screen: {name}')
    ax.legend();ax.grid(alpha=.3);fig.tight_layout();fig.savefig(out/'response_screen.png',dpi=180);plt.close(fig)
    print(json.dumps(result,indent=2),flush=True)
    return result

if __name__=='__main__':
    import sys
    name=sys.argv[1]
    p=from_dict(json.loads((ROOT/'geometry/taper/output'/f'{name}.json').read_text())['second'])
    case(name,p)
