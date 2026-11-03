"""Run and archive a generated single-cell C_add sweep, then fit L0 and C0/Cg."""

from pathlib import Path
from datetime import datetime, timezone
from hashlib import sha256
import json
import re
import shutil
import subprocess

import numpy as np
from scipy.stats import linregress

import constants as cfg
from circuit_common import read_snp, write_sonx, magnitude_db, save_figure, plot_s21, single_sweep_resonance
from circuit_models import s21_single
from extract_circuit_params import (
    C_ADD_SWEEP_PF, CAPACITANCE_C_ADD_PF, MASK_HALF_WIDTH_HZ, extract_capacitances,
)
import matplotlib.pyplot as plt

FOLDER = cfg.SONNET_FILES_DIR / 'generated_ncell_1'
PROJECT = FOLDER / 'generated_inductance.sonx'
SNP = PROJECT.with_suffix('.snp')
OUTPUT = cfg.OUTPUT_DATA_DIR / 'generated_inductance'


def simulate(capacitance, label, *, project_file=None, output_dir=None,solver_env=None,cache_only=False):
    project = Path(project_file).resolve() if project_file is not None else PROJECT
    output = Path(output_dir).resolve() if output_dir is not None else OUTPUT
    folder, snp = project.parent, project.with_suffix('.snp')
    write_sonx({'C_add':str(capacitance)}, project)
    shutil.copyfile(project, output / f'{label}.sonx')
    print(f'Running {label}: C_add={capacitance*1000:g} fF', flush=True)
    with (output / f'{label}.log').open('w') as log:
        from geometry.taper.solver_guard import run as guarded_solver
        command=[str(cfg.EM_EXE),'-v']+(['-AbsNoDiscrete'] if cache_only else [])+[str(project)]
        process=guarded_solver(command,cwd=folder,stdout=log,env=solver_env)
    process.check_returncode()
    log_text = (output / f'{label}.log').read_text()
    if 'Em errors detected: 0' not in log_text or 'Em simulation completed' not in log_text:
        raise RuntimeError('Sonnet log does not confirm a completed error-free run')
    if cache_only and 'Total time for 0 frequencies:' not in log_text:
        raise RuntimeError('Cache-only inductance check did not confirm zero new EM frequencies')
    header = snp.read_text(encoding='utf-8')
    match = re.search(r'^!< C_add = (\S+)', header, re.MULTILINE)
    if not match or not np.isclose(float(match[1]), capacitance, rtol=1e-9, atol=1e-12):
        raise RuntimeError('Touchstone export does not match requested C_add')
    shutil.copyfile(snp, output / f'{label}.snp')
    net = read_snp(snp)
    # Sampled maxima reproduce the existing extraction method.
    frequency = single_sweep_resonance(net)
    print(f'{label}: resonance {frequency/1e9:.7f} GHz', flush=True)
    return net, frequency


def main(*, fit_capacitance=True, project_file=None, output_dir=None, c_add_sweep_pf=None, target_l0_nh=None,solver_env=None,cache_only=False,reuse_cache_after_first=False):
    project = Path(project_file).resolve() if project_file is not None else PROJECT
    output = Path(output_dir).resolve() if output_dir is not None else OUTPUT
    snp = project.with_suffix('.snp')
    target_l0_nh = cfg.TARGET_L0_NH if target_l0_nh is None else target_l0_nh
    if not np.isfinite(target_l0_nh) or target_l0_nh <= 0:
        raise ValueError('Inductance target must be finite and positive')
    c_add_sweep = C_ADD_SWEEP_PF if c_add_sweep_pf is None else tuple(c_add_sweep_pf)
    if len(c_add_sweep) < 3 or not np.all(np.isfinite(c_add_sweep)) or not np.all(np.diff(c_add_sweep)>0) or min(c_add_sweep)<=0:
        raise ValueError('Use at least three positive increasing C_add values in pF')
    # Preserve existing native tuner's two-argument simulate override.
    run_point = simulate if project_file is None and output_dir is None and solver_env is None and not cache_only else (
        lambda value, label: simulate(value, label, project_file=project, output_dir=output,solver_env=solver_env,cache_only=cache_only))
    output.mkdir(parents=True, exist_ok=True)
    records, nets = [], []
    for index,capacitance in enumerate(c_add_sweep):
        label = f'cadd_{capacitance*1000:g}ff'
        if reuse_cache_after_first:
            net,frequency=simulate(capacitance,label,project_file=project,output_dir=output,
                                   solver_env=solver_env,cache_only=cache_only or index>0)
        else:
            net, frequency = run_point(capacitance,label)
        records.append({'C_add_fF':capacitance*1000, 'frequency_Hz':frequency})
        records[-1]['project_sha256'] = sha256((output/f'{label}.sonx').read_bytes()).hexdigest()
        records[-1]['snp_sha256'] = sha256((output/f'{label}.snp').read_bytes()).hexdigest()
        nets.append(net)
        (output/'sweep_points.json').write_text(json.dumps(records,indent=2)+'\n')
    freqs = np.array([r['frequency_Hz'] for r in records])
    if not np.all(np.diff(freqs)<0):
        raise RuntimeError('Resonance does not decrease monotonically with C_add')
    x = np.array(c_add_sweep)*1e-12
    y = 1/(2*np.pi*freqs)**2
    fit = linregress(x,y)
    if fit.slope <= 0:
        raise RuntimeError('Nonphysical inductance slope')
    result = {'L0_nH':fit.slope*1e9, 'R_squared':fit.rvalue**2,
              'regression_standard_error_nH':fit.stderr*1e9,
              'intercept_s2':fit.intercept, 'sweep':records,
              'target_L0_nH':target_l0_nh,
              'target_difference_percent':100*(fit.slope/(target_l0_nh*1e-9)-1),
              'method':'Slope of 1/(2*pi*f_peak)^2 vs C_add; sampled maxima, fixed geometry.'}
    result['adjacent_pair_L0_nH'] = (np.diff(y)/np.diff(x)*1e9).tolist()
    fitted_freqs = 1/(2*np.pi*np.sqrt(fit.intercept+fit.slope*x))
    result['max_frequency_residual_MHz'] = float(np.max(abs(fitted_freqs-freqs))/1e6)
    result['export_frequency_spacing_MHz'] = float(np.median(np.diff(nets[0].f))/1e6)
    (output/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    fig,ax=plt.subplots(figsize=(7,5))
    ax.plot(x*1e15,y*1e21,'ko',label='Sonnet')
    ax.plot(x*1e15,(fit.intercept+fit.slope*x)*1e21,'r--',label=f'R² = {fit.rvalue**2:.7f}')
    ax.set(xlabel='C_add (fF)',ylabel='1/omega² (10^-21 s²)',title=f'L0 = {fit.slope*1e9:.5f} nH')
    ax.legend(); ax.grid(alpha=.3)
    save_figure(fig,output/'inductance_fit.png')
    fig,ax=plt.subplots(figsize=(9,5))
    for record,net in zip(records,nets):
        ax.plot(net.f/1e9,magnitude_db(net.s[:,1,0]),label=f"{record['C_add_fF']:g} fF")
    ax.set(xlabel='Frequency (GHz)',ylabel='|S21| (dB)',ylim=(-60,2),title='Generated single-cell C_add sweep')
    ax.legend(); ax.grid(alpha=.3)
    save_figure(fig,output/'sweep_s21.png')

    if not fit_capacitance:
        return result

    # Leave a tiny C_add, matching the original single-cell capacitance workflow.
    net,peak=run_point(CAPACITANCE_C_ADD_PF,'cadd_0.1ff')
    c0,cg=extract_capacitances(fit.slope,project,snp,cfg.EM_EXE,simulate=False)
    response=s21_single(2*np.pi*net.f,c0*1e-15,cg*1e-15,fit.slope)
    window = (max(net.f[0], peak-MASK_HALF_WIDTH_HZ), min(net.f[-1], peak+MASK_HALF_WIDTH_HZ))
    plot_s21(net,response,window,output,fit.slope,c0,cg)
    result.update({'C0_fF':float(c0),'Cg_fF':float(cg),'final_C_add_fF':.1,
                   'complex_s21_rms_error':float(np.sqrt(np.mean(abs(response-net.s[:,1,0])**2)))})
    (output/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__ == '__main__':
    OUTPUT.mkdir(parents=True, exist_ok=True)
    status = {'state':'running', 'started_utc':datetime.now(timezone.utc).isoformat()}
    status_path = OUTPUT/'run_status.json'
    status_path.write_text(json.dumps(status,indent=2)+'\n')
    try:
        main()
    except Exception as error:
        status.update(state='failed', error=str(error))
        raise
    else:
        status.update(state='complete', results_sha256=sha256((OUTPUT/'results.json').read_bytes()).hexdigest())
    finally:
        status['updated_utc'] = datetime.now(timezone.utc).isoformat()
        status_path.write_text(json.dumps(status,indent=2)+'\n')
