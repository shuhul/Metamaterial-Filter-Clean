"""Geometry iteration scored against the specified circuit, not a refitted circuit.

Run as a module. Every candidate gets separate immutable geometry/native files.
Only a matched fixture delay may change in the target comparison. All target
L/C values, 50-ohm terminations, and the 4.65--5.40 GHz scoring band stay fixed.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
from dataclasses import replace
from pathlib import Path
import json
import xml.etree.ElementTree as ET

import numpy as np
from scipy.optimize import minimize_scalar
from circuit_common import read_snp
from geometry.specifications import (TAPER_1, TAPER_2, UNIT_CELL, specification,
    S21_COMPLEX_RMS_CUTOFF, TARGET_S21_FREQUENCY_RANGE_GHZ)
from geometry.taper.clearance_extraction import prepare, measure
from geometry.taper.second_layout import from_dict
from geometry.taper.corrected_iteration import ROOT, review
from geometry.taper.progressive_fitting import transmission
from geometry.taper.sonnet import simulate, configure
from geometry.prepare_sonnet import tag, write_xml
from geometry.taper.run_bridged_stage import digest

CAMPAIGN = 'target_response_2026-10-25'
OUT = ROOT / 'output_data/taper_exploration' / CAMPAIGN
START, STOP = TARGET_S21_FREQUENCY_RANGE_GHZ


def parameters(width=400, wire=88, span=140, gap=2, pairs=36, right_thickness=20,
               junction_length=584):
    source = ROOT / ('geometry/taper/output/clearance_extraction_2026-10-25/'
                     'frame700_wire88_left214/idc/geometry.json')
    p = from_dict(json.loads(source.read_text())['second'])
    return replace(p, taper=replace(p.taper, w=width, l_right=junction_length, t_right=right_thickness,
                   meander=replace(p.taper.meander, t=wire)),
                   taper1=replace(p.taper1,l_left=junction_length),
                   idc_span=span, idc_pairs=pairs, d_idc_top=2,
                   equal_idc_sides=True, d_t2_t1=gap)


def target_response(f, delay=0):
    return transmission(f,
        [TAPER_2.L0_nH, TAPER_1.L0_nH, UNIT_CELL.L0_nH, UNIT_CELL.L0_nH],
        [TAPER_2.C0_fF, TAPER_1.C0_fF, UNIT_CELL.C0_fF, UNIT_CELL.C0_fF],
        [TAPER_2.Cg_input_fF, TAPER_1.Cg_input_fF,
         UNIT_CELL.Cg_input_fF, UNIT_CELL.Cg_input_fF], delay)


def score(source, folder, *, fit_result=None):
    folder=Path(folder); folder.mkdir(parents=True, exist_ok=True)
    net=read_snp(source)
    if not np.isclose(net.f[0], START*1e9) or not np.isclose(net.f[-1], STOP*1e9):
        raise ValueError('Target scoring requires the fixed complete frequency band')
    if len(net.f)<7500 or not np.allclose(net.f,np.linspace(START*1e9,STOP*1e9,len(net.f)),rtol=0,atol=1):
        raise ValueError('Target scoring requires a dense uniform grid over the fixed band')
    f=net.f; actual=net.s[:,1,0]; zero=target_response(f)
    # A global grid prevents phase-wrap minima from selecting a bad delay.
    objective=lambda tau: float(np.mean(abs(zero*np.exp(-2j*np.pi*f*tau*1e-12)-actual)**2))
    grid=np.linspace(-100,100,801); i=int(np.argmin([objective(t) for t in grid]))
    opt=minimize_scalar(objective, bounds=(max(-100,grid[i]-.5),min(100,grid[i]+.5)), method='bounded')
    prediction=target_response(f,opt.x)
    result=dict(metric='complex_S21_RMS_against_fixed_target_circuit',
        complex_RMS=float(np.sqrt(opt.fun)), magnitude_RMS=float(np.sqrt(np.mean((abs(prediction)-abs(actual))**2))),
        maximum_complex_error=float(np.max(abs(prediction-actual))),
        pass_target=bool(opt.success and np.sqrt(opt.fun)<=S21_COMPLEX_RMS_CUTOFF), cutoff=S21_COMPLEX_RMS_CUTOFF,
        delay_optimizer_success=bool(opt.success),
        fixture_delay_ps=float(opt.x), fixture_delay_bounds_ps=[-100,100],
        frequency_range_GHz=[START,STOP], frequency_samples=len(f),
        all_electrical_targets_fixed=True, specification=specification(),
        source=str(Path(source).resolve()), source_sha256=digest(source),
        scope='One matched fixture delay is optimized; no electrical parameters or frequency shifts are fitted.')
    (folder/'target_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez(folder/'target_comparison.npz',f=f,actual=actual,target=prediction)
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(14,8))
    ax.plot(f/1e9,abs(actual),lw=2.8,color='#2768aa',label='Sonnet EM')
    ax.plot(f/1e9,abs(prediction),'--',lw=2.5,color='#db6d28',label='Specified target circuit')
    ax.set(xlabel='Frequency (GHz)',ylabel=r'$|S_{21}|$',xlim=(START,STOP),ylim=(0,1.25),
           title='Taper 2 + Taper 1: Simulation vs Target Circuit')
    ax.xaxis.label.set_size(20);ax.yaxis.label.set_size(23);ax.title.set_size(23)
    ax.tick_params(labelsize=16);ax.grid(alpha=.2);ax.legend(fontsize=16,loc='upper right')
    ax.text(.025,.97,f"Target complex RMS = {result['complex_RMS']:.4f} (limit 0.05)\n"
            f"Magnitude RMS = {result['magnitude_RMS']:.4f}\n"
            f"Fixture delay = {opt.x:.2f} ps; all target L/C fixed",
            transform=ax.transAxes,va='top',fontsize=15,
            bbox=dict(boxstyle='round,pad=.5',facecolor='white',edgecolor='#aaaaaa'))
    fig.tight_layout()
    for suffix in ('png','pdf'):fig.savefig(folder/f's21_target_comparison.{suffix}',dpi=180)
    plt.close(fig)
    print('TARGET SCORE',json.dumps(result),flush=True)
    return result


def response(name,p,*,inductance=None,junction=None,recess_t1=0):
    project,out,meta=prepare(name,p,'idc',campaign=CAMPAIGN,junction=junction,recess_t1=recess_t1)
    marker=out/'run_status.json'
    if marker.exists() and json.loads(marker.read_text()).get('status')=='complete':
        saved=json.loads(marker.read_text())
        if saved['project_sha256']!=digest(project) or saved['snp_sha256']!=digest(project.with_suffix('.snp')):
            raise ValueError('Saved simulation changed')
        source=project.with_suffix('.snp')
    else:
        tree=ET.parse(project);root=tree.getroot();sweeps=root.find(tag('Sweeps'))
        for child in list(sweeps):sweeps.remove(child)
        fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
        ET.SubElement(fs,tag('Sweep'),Type='ABS_ENTRY',On='TRUE',Adaptive='TRUE',
                      Start=str(START),Stop=str(STOP),Target='7500')
        for q in root.iter(tag('File')):
            if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','TRUE')
        write_xml(tree,project);configure(project,START,STOP,7500)
        source=simulate(project,out)
    result=score(source,out)
    if inductance is not None:
        # Heuristic only: these fitted capacitances do not determine success.
        review(source,out/'fit',inductance,idc=True)
    return result


if __name__=='__main__':
    import sys
    names=[s for s in sys.argv[1:] if not s.startswith('--')]
    name=names[0] if names else 'equal_w400_wire88_span140_gap2'
    p=parameters()
    if '--prepare' in sys.argv:prepare(name,p,'idc',campaign=CAMPAIGN)
    elif '--calibrate' in sys.argv:measure(name,p,'idc_Lcal',campaign=CAMPAIGN)
    else:response(name,p)
