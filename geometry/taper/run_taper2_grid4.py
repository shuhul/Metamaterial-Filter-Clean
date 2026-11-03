"""Run the prepared coarse T2 once, fit three Cs, and record timing evidence."""
import os
os.environ['SONNET_MAX_JOBS']='1'
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
from hashlib import sha256
import json,re,time
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from geometry.taper import taper2_grid4 as layout
from geometry.taper.sonnet import simulate
from geometry.taper.progressive_fitting import fit
from circuit_common import read_snp

OUT=layout.ROOT/'output_data/taper_exploration/three_tapers_v3_T2_grid4_2026-10-29/best_003853_grid4/response'


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    project=layout.DEST/'generated_ncell_3.sonx'
    report=json.loads((layout.DEST/'validation.json').read_text())
    if sha256(project.read_bytes()).hexdigest()!=report['files_sha256'][project.name]:
        raise ValueError('Reviewed coarse project changed')
    policy=layout.fit_policy()
    (layout.DEST/'fit_policy.json').write_text(json.dumps(policy,indent=2)+'\n')
    report['files_sha256']['fit_policy.json']=sha256((layout.DEST/'fit_policy.json').read_bytes()).hexdigest()
    (layout.DEST/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    manifest=OUT/'benchmark_start.json';status=OUT/'run_status.json'
    if not manifest.exists():
        manifest.write_text(json.dumps(dict(project=str(project),project_sha256=sha256(project.read_bytes()).hexdigest(),
            sondata_cache_existed_before=(project.parent/'sondata').exists(),
            snp_existed_before=project.with_suffix('.snp').exists(),
            mesh_um=[4,4],frequency_range_GHz=[4.65,5.4],samples=7501,thread_limit=4,
            policy=policy),indent=2)+'\n')
    if status.exists() and json.loads(status.read_text()).get('status')=='complete':
        saved=json.loads(status.read_text())
        if sha256(project.with_suffix('.snp').read_bytes()).hexdigest()!=saved['snp_sha256']:
            raise ValueError('Completed coarse response changed')
        snp=project.with_suffix('.snp')
    else:snp=simulate(project,OUT)
    start=time.perf_counter()
    fitted=fit(snp,OUT/'three_cap_fit',policy['fixed_L_nH'],policy['initial_C_fF'],policy['initial_G_fF'],
        [tuple(q) for q in policy['free_capacitances']],title='T2: Three-Capacitance Fit on 4 um Mesh')
    fit_seconds=time.perf_counter()-start
    baseline=layout.OLD_FIT.parent.parent
    old_resource=json.loads((baseline/'em.log.resources.json').read_text())
    new_resource=json.loads((OUT/'em.log.resources.json').read_text())
    timing=dict(baseline_2um_solver_seconds=old_resource['elapsed_seconds'],
        coarse_4um_solver_seconds=new_resource['elapsed_seconds'],
        solver_speedup=old_resource['elapsed_seconds']/new_resource['elapsed_seconds'],
        coarse_fit_processing_seconds=fit_seconds,
        coarse_solver_plus_fit_seconds=new_resource['elapsed_seconds']+fit_seconds,
        benchmark_start=json.loads(manifest.read_text()),
        baseline_log=str(baseline/'em.log'),coarse_log=str(OUT/'em.log'),
        baseline_subsections=[int(x.replace(',','')) for x in re.findall(r'Circuit requires ([\d,]+) subsections', (baseline/'em.log').read_text())],
        coarse_subsections=[int(x.replace(',','')) for x in re.findall(r'Circuit requires ([\d,]+) subsections', (OUT/'em.log').read_text())],
        baseline_new_EM_frequency_counts=re.findall(r'Total time for (\d+) frequencies:',(baseline/'em.log').read_text()),
        coarse_new_EM_frequency_counts=re.findall(r'Total time for (\d+) frequencies:',(OUT/'em.log').read_text()),
        comparison_limit='Historical baseline versus current run; same band, sampling and four threads. Coarse metal is snapped, including 2-to4um tip-gap and strip changes. Not a controlled identical-metal mesh-only benchmark.',
        fit=fitted)
    (OUT/'timing_comparison.json').write_text(json.dumps(timing,indent=2)+'\n')
    data=np.load(OUT/'three_cap_fit/response.npz');old=read_snp(layout.SOURCE.with_suffix('.snp'))
    fig,ax=plt.subplots(figsize=(14,8))
    ax.plot(old.f/1e9,abs(old.s[:,1,0]),color='#888',lw=1.8,label='Original 2 um simulation')
    ax.plot(data['f']/1e9,abs(data['actual'][:,1,0]),color='#2768aa',lw=2.8,label='Snapped 4 um simulation')
    ax.plot(data['f']/1e9,abs(data['model'][:,1,0]),'k--',lw=2.1,label='Three-capacitance complex fit')
    text=f'T2 C0 = {fitted["C_fF"][0]:.3f} fF\nT2 input Cg = {fitted["G_fF"][0]:.3f} fF\nT2-to-T1 Cg = {fitted["G_fF"][1]:.3f} fF\nFit RMS = {fitted["complex_RMS"]:.6f}\nOlder capacitances and all L fixed\nFixture delay fitted'
    ax.text(.025,.97,text,transform=ax.transAxes,va='top',fontsize=15,bbox=dict(boxstyle='round,pad=.5',fc='white',ec='gray'))
    ax.set(xlabel='Frequency (GHz)',ylabel=r'$|S_{21}|$',title='T2 Transmission: 4 um Simulation',
           xlim=(4.65,5.4),ylim=(0,max(abs(data['actual'][:,1,0]).max(),abs(old.s[:,1,0]).max())*1.5))
    ax.xaxis.label.set_size(20);ax.yaxis.label.set_size(23);ax.title.set_size(23);ax.tick_params(labelsize=16)
    ax.grid(alpha=.2);ax.legend(fontsize=12,loc='upper right');fig.tight_layout()
    for ext in ('png','pdf'):fig.savefig(OUT/f's21_comparison.{ext}',dpi=190)
    plt.close(fig)
    print(json.dumps(timing,indent=2),flush=True)


if __name__=='__main__':main()
