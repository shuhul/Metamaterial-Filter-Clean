"""T2 iteration with three fitted-C targets and a fixed older-circuit reference.

2um mesh, one4-thread native solve. Theoretical-S21 error is diagnostic only.
Historical projects/results are immutable; each geometry gets a new case.
"""
import os
os.environ['SONNET_MAX_JOBS']='1'
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
from dataclasses import asdict,replace
from hashlib import sha256
import json
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from geometry.taper import taper2_plain as workflow
from geometry.taper.taper2_grid4 import ROOT,fit_policy,OLDER_REFERENCE
from geometry.taper.progressive_fitting import fit,transmission
from geometry.taper.sonnet import configure,simulate
from geometry.prepare_sonnet import tag
from geometry.specifications import TAPER_2,TAPER_1,TAPER_2_CAPACITANCE_TOLERANCE_PERCENT
from circuit_common import read_snp

CAMPAIGN='three_tapers_v3_T2_cap_targets_2026-10-29'
OUT=ROOT/'output_data/taper_exploration'/CAMPAIGN
OLD_OUT=ROOT/'output_data/taper_exploration'/workflow.CAMPAIGN
NATIVE=ROOT/'sonnet_files/taper_exploration'/workflow.CAMPAIGN
TARGET=np.array([TAPER_2.C0_fF,TAPER_2.Cg_input_fF,TAPER_1.Cg_input_fF])
NAMES=['T2 C0','T2 input Cg','T2-to-T1 Cg']
POLICY=fit_policy()


def fit_case(case,*,l0=None,provisional=False):
    project=NATIVE/case/'response/generated_ncell_3.sonx';snp=project.with_suffix('.snp')
    previous=OLD_OUT/case/'response/extraction/fit.json'
    if l0 is None:
        if not previous.exists():raise ValueError('Need explicit or declared provisional L')
        l0=json.loads(previous.read_text())['L_nH'][0]
    ls=list(POLICY['fixed_L_nH']);ls[0]=l0
    tree=ET.parse(project);box=tree.find(f'{tag("Geometry")}/{tag("Box")}')
    size=box.find(tag('Size'));cells=box.find(tag('NumCells'))
    if any(float(size.get(a))/float(cells.get(a))!=2 for a in ('X','Y')):raise ValueError('Not a2um mesh')
    net=read_snp(snp)
    if len(net.f)!=7501 or not np.allclose([net.f[0],net.f[-1]],[4.65e9,5.4e9],rtol=0,atol=1):raise ValueError('Wrong band/sample grid')
    out=OUT/case;out.mkdir(parents=True,exist_ok=True)
    seeds=[TARGET.tolist()]
    if previous.exists():
        old=json.loads(previous.read_text());seeds.append([old['C_fF'][0],old['G_fF'][0],old['G_fF'][1]])
    fitted=fit(snp,out/'fit',ls,POLICY['initial_C_fF'],POLICY['initial_G_fF'],
        [('C',0),('G',0),('G',1)],title='T2: Three-Capacitance Extraction',capacitance_seeds=seeds)
    values=np.array([fitted['C_fF'][0],fitted['G_fF'][0],fitted['G_fF'][1]])
    errors=100*(values/TARGET-1)
    result=dict(case=case,project=str(project),fit=str(out/'fit/fit.json'),
        values_fF=values.tolist(),targets_fF=TARGET.tolist(),target_errors_percent=errors.tolist(),
        max_capacitance_error_percent=float(max(abs(errors))),
        capacitance_pass=bool(np.all(abs(errors)<=TAPER_2_CAPACITANCE_TOLERANCE_PERCENT)),
        fit_RMS=fitted['complex_RMS'],
        fit_quality_pass=bool(fitted['complex_RMS']<=.025 and fitted['optimizer_success']
                              and fitted['interior_free_parameters']),
        L0_nH=l0,provisional_L=provisional,fixed_older_reference=str(OLDER_REFERENCE),
        fixed_older_reference_sha256=sha256(OLDER_REFERENCE.read_bytes()).hexdigest(),
        project_sha256=sha256(project.read_bytes()).hexdigest(),snp_sha256=sha256(snp.read_bytes()).hexdigest(),
        objective='Three fitted taper capacitances within5% of nominal. Fit-to-Sonnet RMS checks model quality; theoretical-target S21 RMS is not acceptance.',
        geometry=json.loads((project.parent/'geometry.json').read_text())['trial_settings'])
    (out/'result.json').write_text(json.dumps(result,indent=2)+'\n')
    fig,ax=plt.subplots(figsize=(14,8));data=np.load(out/'fit/response.npz')
    ax.plot(net.f/1e9,abs(net.s[:,1,0]),lw=2.8,color='#2768aa',label='Sonnet EM (2um mesh)')
    ax.plot(data['f']/1e9,abs(data['model'][:,1,0]),'k--',lw=2.2,label='Three-capacitance complex fit')
    text='\n'.join(f'{n} = {v:.3f} fF (target {t:g}; {e:+.2f}%)' for n,v,t,e in zip(NAMES,values,TARGET,errors))
    text+=f'\nL T2 = {l0:.5f} nH'+(' (provisional)' if provisional else ' (fixed calibration)')
    text+=f'\nFit-to-Sonnet RMS = {fitted["complex_RMS"]:.6f}\nCapacitance tolerance: +/-5%\nOlder circuit values fixed; fixture delay fitted'
    ax.text(.025,.97,text,transform=ax.transAxes,va='top',fontsize=14,bbox=dict(boxstyle='round,pad=.5',fc='white',ec='gray'))
    ax.set(xlabel='Frequency (GHz)',ylabel=r'$|S_{21}|$',title='Second Taper Cell: Capacitance Target Check',
           xlim=(4.65,5.4),ylim=(0,max(abs(net.s[:,1,0]).max(),abs(data['model'][:,1,0]).max())*1.6))
    ax.xaxis.label.set_size(20);ax.yaxis.label.set_size(23);ax.title.set_size(23);ax.tick_params(labelsize=16)
    ax.grid(alpha=.2);ax.legend(fontsize=12,loc='upper right');fig.tight_layout()
    for ext in ('png','pdf'):fig.savefig(out/f's21_fit.{ext}',dpi=190)
    plt.close(fig);update_ledger();print('CAP_RESULT',case,values.tolist(),errors.tolist(),fitted['complex_RMS'],flush=True)
    return result


def update_ledger():
    rows=[json.loads(p.read_text()) for p in OUT.glob('*/result.json')]
    rows.sort(key=lambda r:r['max_capacitance_error_percent'])
    (OUT/'ledger.json').write_text(json.dumps(dict(objective='three new fitted Cs within5%',
        fixed_older_reference=str(OLDER_REFERENCE),cases=rows),indent=2)+'\n')


def refit_saved():
    OUT.mkdir(parents=True,exist_ok=True)
    for meta in sorted(NATIVE.glob('*/response/geometry.json')):
        settings=json.loads(meta.read_text()).get('trial_settings',{});c=settings.get('coupling',{})
        if (settings.get('meander',{}).get('w')!=148 or c.get('d_regular')!=74
                or c.get('t_regular_ground')!=2 or not (meta.parent/'generated_ncell_3.snp').exists()):continue
        case=meta.parts[-3]
        if not (OLD_OUT/case/'response/extraction/fit.json').exists():continue
        fit_case(case)


def run_case(case,settings,*,l0=POLICY['fixed_L_nH'][0]):
    project,_,meta=workflow.prepare(case,settings)
    out=OUT/case;out.mkdir(parents=True,exist_ok=True)
    (out/'proposal.json').write_text(json.dumps(dict(settings=asdict(settings),source_geometry=meta,
        provisional_L_nH=l0,fixed_older_reference=str(OLDER_REFERENCE)),indent=2)+'\n')
    status=out/'run_status.json'
    if status.exists() and json.loads(status.read_text()).get('status')=='complete':
        record=json.loads(status.read_text())
        if record['project_sha256']!=sha256(project.read_bytes()).hexdigest() or record['snp_sha256']!=sha256(project.with_suffix('.snp').read_bytes()).hexdigest():
            raise ValueError('Completed run changed')
    else:
        configure(project,4.65,5.4,7500);simulate(project,out)
    return fit_case(case,l0=l0,provisional=True)


if __name__=='__main__':refit_saved()
