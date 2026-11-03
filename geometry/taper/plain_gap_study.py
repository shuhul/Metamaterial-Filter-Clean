"""Controlled, one-parameter study of the user-approved plain 4 um layout.

No ground cuts, junction IDC, pads, or changes to T1/regular metal are exposed.
The saved baseline is reused. All EM jobs use the common four-thread guard.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from dataclasses import replace, asdict
from pathlib import Path
from datetime import datetime, timezone
import json
import shutil
import xml.etree.ElementTree as ET
import gdstk
import numpy as np

from circuit_common import read_snp
from geometry.prepare_sonnet import tag, write_xml, BIN
from geometry.taper.second_layout import from_dict
from geometry.taper.clearance_extraction import prepare
from geometry.taper.target_iteration import ROOT, score, START, STOP
from geometry.taper.sonnet import simulate, configure
from geometry.taper.solver_guard import run as guarded_solver
from geometry.taper.solver_resources import solver_environment
from geometry.taper.run_bridged_stage import digest

CAMPAIGN='plain_gap_sensitivity_2026-10-25'
NATIVE=ROOT/'sonnet_files/taper_exploration'/CAMPAIGN
OUT=ROOT/'output_data/taper_exploration'/CAMPAIGN
BASE_NAME='equal_w400_wire60_span140_gap4'
BASE_NATIVE=ROOT/'sonnet_files/taper_exploration/target_response_2026-10-25'/BASE_NAME/'idc/generated_ncell_3.sonx'
BASE_GDS=ROOT/'geometry/taper/output/target_response_2026-10-25'/BASE_NAME/'idc/geometry.gds'
BASE_OUT=ROOT/'output_data/taper_exploration/target_response_2026-10-25'/BASE_NAME/'idc'


def baseline():
    return from_dict(json.loads(BASE_GDS.with_suffix('.json').read_text())['second'])


def configurations():
    p=baseline()
    return {
        'input_span156':replace(p,idc_span=156),
        'facing_thickness24':replace(p,taper=replace(p.taper,t_right=24)),
        'frame_width416':replace(p,taper=replace(p.taper,w=416)),
    }


def leaves(value,prefix=''):
    result={}
    for k,v in value.items():
        name=f'{prefix}.{k}' if prefix else k
        if isinstance(v,dict):result.update(leaves(v,name))
        else:result[name]=v
    return result


def check_parameters(name,p):
    expected={'input_span156':'idc_span','facing_thickness24':'taper.t_right',
              'frame_width416':'taper.w'}[name]
    original=leaves(asdict(baseline()));actual=leaves(asdict(p))
    changed=[k for k in original if original[k]!=actual[k]]
    if changed!=[expected]:raise ValueError(f'Expected one parameter: {changed}')
    if p.d_t2_t1!=4 or p.taper.l_right!=584 or p.taper1.l_left!=584 or not p.equal_idc_sides:
        raise ValueError('Approved coupling gap/equal lengths must remain fixed')
    return {expected:{'baseline':original[expected],'variant':actual[expected]}}


def prepare_case(name):
    p=configurations()[name];change=check_parameters(name,p)
    project,out,meta=prepare(name,p,'idc',campaign=CAMPAIGN)
    if 'junction_idc' in meta or meta.get('T1_top_ground_recess_um',0):
        raise ValueError('Rejected topology entered plain-gap study')
    gd=ROOT/'geometry/taper/output'/CAMPAIGN/name/'idc/geometry.gds'
    metal=gdstk.read_gds(str(gd)).top_level()[0].get_polygons()
    top=gdstk.rectangle((0,0),(meta['w_sim'],448))
    if gdstk.boolean([top],metal,'not',precision=.001):
        raise ValueError('Ground plane is not the complete straight rectangle')
    (out/'intent.json').write_text(json.dumps(dict(changed_parameter=change,
        baseline_project=str(BASE_NATIVE),baseline_sha256=digest(BASE_NATIVE),
        ground_edge='straight',T2_T1_gap_um=4,facing_leg_lengths_um=[584,584],
        T1_and_downstream_native_lock='passed',extra_junction_IDC=False,
        input_IDC_bus_lengths_equal=True,grid_um=[2,2],goal='Fixed-target complex S21 RMS <= 0.05'),indent=2)+'\n')
    return project,out


def run_case(name):
    project,out=prepare_case(name);marker=out/'run_status.json'
    if marker.exists() and json.loads(marker.read_text()).get('status')=='complete':
        record=json.loads(marker.read_text())
        if record['project_sha256']!=digest(project) or record['snp_sha256']!=digest(project.with_suffix('.snp')):
            raise ValueError('Completed study data changed')
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
    score(source,out)
    from geometry.taper.response_diagnostics import analyze, summarize_study
    analyze(source,out)
    summarize_study()
    return source


def run_currents():
    """Same baseline metal/ports; three direct frequencies with saved currents."""
    dest=NATIVE/'baseline_currents';out=OUT/'baseline_currents'
    dest.mkdir(parents=True,exist_ok=True);out.mkdir(parents=True,exist_ok=True)
    project=dest/'generated_ncell_3.sonx';marker=out/'run_status.json'
    frequencies=[4.8441,4.9616,5.0791]  # Baseline lower peak, midpoint, upper peak.
    if marker.exists() and json.loads(marker.read_text()).get('status')=='complete':
        r=json.loads(marker.read_text())
        if r['project_sha256']!=digest(project) or r['snp_sha256']!=digest(project.with_suffix('.snp')):
            raise ValueError('Current diagnostic data changed')
        return project
    if not project.exists():
        shutil.copyfile(BASE_NATIVE,project)
        tree=ET.parse(project);root=tree.getroot();control=root.find(tag('Control'))
        ET.SubElement(control,tag('Currents'),On='TRUE')
        sweeps=root.find(tag('Sweeps'))
        for child in list(sweeps):sweeps.remove(child)
        fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
        for f in frequencies:
            ET.SubElement(fs,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(f),Stop=str(f),Step='1')
        for q in root.iter(tag('File')):
            if q.get('Type')=='TOUCH':
                q.set('Name',project.with_suffix('.snp').name);q.set('UseInterpolatedData','FALSE')
        write_xml(tree,project)
    # XML element identity confirms absolutely no metal/port/stack edits.
    old=ET.parse(BASE_NATIVE).getroot().find(tag('Geometry'))
    new=ET.parse(project).getroot().find(tag('Geometry'))
    if ET.tostring(old)!=ET.tostring(new):raise ValueError('Current fixture changed physical geometry')
    record=dict(status='running',project=str(project),project_sha256=digest(project),
                baseline_sha256=digest(BASE_NATIVE),frequencies_GHz=frequencies,
                geometry_stack_ports_identical=True,started_utc=datetime.now(timezone.utc).isoformat())
    marker.write_text(json.dumps(record,indent=2)+'\n')
    try:
        with (out/'em.log').open('w') as log:
            guarded_solver([str(BIN/'em.exe'),'-v',str(project)],cwd=dest,stdout=log,
                           env=solver_environment()).check_returncode()
        text=(out/'em.log').read_text()
        if 'Em errors detected: 0' not in text or 'Em simulation completed' not in text:
            raise RuntimeError('Current diagnostic simulation did not complete cleanly')
        net=read_snp(project.with_suffix('.snp'));original=read_snp(BASE_NATIVE.with_suffix('.snp'))
        if len(net.f)!=3 or not np.allclose(net.f,np.array(frequencies)*1e9,rtol=0,atol=1):
            raise ValueError('Current frequencies differ from request')
        differences=[]
        for i,f in enumerate(net.f):
            ix=int(np.argmin(abs(original.f-f)))
            differences.append(float(np.max(abs(net.s[i]-original.s[ix]))))
        record.update(status='complete',snp_sha256=digest(project.with_suffix('.snp')),
            direct_vs_saved_ABS_max_S_errors=differences,
            reference_check_pass=bool(max(differences)<.01),
            ended_utc=datetime.now(timezone.utc).isoformat())
        marker.write_text(json.dumps(record,indent=2)+'\n')
        print('CURRENTS COMPLETE',json.dumps(record),flush=True)
        return project
    except Exception as exc:
        record.update(status='failed',reason=str(exc));marker.write_text(json.dumps(record,indent=2)+'\n')
        raise


if __name__=='__main__':
    import sys
    if '--currents' in sys.argv:run_currents()
    elif '--prepare' in sys.argv:
        for name in configurations():prepare_case(name)
    else:
        for name in configurations():run_case(name)
