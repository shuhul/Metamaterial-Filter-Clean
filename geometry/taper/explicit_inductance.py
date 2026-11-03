"""Inductance from explicit Sonnet capacitors and cached capacitor changes.

The original fixture is immutable. Each check has its own native project and
starts one native CAP cache, reused for all later capacitor values. No Python admittance interpolant
or numerical capacitor termination is used to obtain the resonance frequencies.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
from hashlib import sha256
from datetime import datetime,timezone
import json,re,shutil,time
import xml.etree.ElementTree as ET
import numpy as np
from geometry.prepare_sonnet import tag,write_xml
from geometry.taper.sonnet import add_capacitor
import extract_generated_inductance as legacy

ROOT=Path(__file__).resolve().parents[2]
T1_CHECK=ROOT/'output_data/taper_exploration/explicit_L_checks_2026-10-28/T1/results.json'


def digest(p):return sha256(Path(p).read_bytes()).hexdigest()


def measure(source,native,out,*,name,target_L_nH,caps_fF,metadata=None,samples=1601,
            range_GHz=(2.4,4.0)):
    source=Path(source);native=Path(native);out=Path(out)
    caps=list(map(float,caps_fF));project=native/source.name
    cal_start_GHz,cal_stop_GHz=map(float,range_GHz)
    if not 0<cal_start_GHz<cal_stop_GHz:raise ValueError('Invalid explicit calibration band')
    expected=dict(source=str(source.resolve()),source_sha256=digest(source),
        C_add_fF=caps,range_GHz=[cal_start_GHz,cal_stop_GHz],samples=samples,target_L_nH=target_L_nH,
        method='Explicit Sonnet CAP components; cache-only subsequent capacitor changes')
    marker=out/'completed.json';manifest=out/'manifest.json'
    if marker.exists():
        record=json.loads(marker.read_text())
        if record['configuration']!=expected:raise ValueError('Explicit calibration settings or source changed')
        for path,hashed in record['files_sha256'].items():
            if digest(path)!=hashed:raise ValueError('Explicit calibration source/export changed: '+path)
        return json.loads((out/'results.json').read_text())
    native.mkdir(parents=True,exist_ok=True);out.mkdir(parents=True,exist_ok=True)
    if manifest.exists() and json.loads(manifest.read_text())['configuration']!=expected:
        raise ValueError('Use a new explicit calibration directory for different settings')
    if not project.exists():
        shutil.copy2(source,project)
    tree=ET.parse(project)
    if not list(tree.iter(tag('Component'))):
        if metadata is None:raise ValueError('Metadata needed to attach the explicit capacitor')
        add_capacitor(project,metadata);tree=ET.parse(project)
    components=list(tree.iter(tag('Component')))
    if len(components)!=1:raise ValueError('Expected exactly one shunt calibration component')
    comp=components[0];comp.set('Type','CAP');comp.set('Units','PF')
    for child in list(comp):
        if child.tag==tag('Value'):comp.remove(child)
    ET.SubElement(comp,tag('Value')).text='C_add'
    variables=[v for v in tree.iter(tag('Variable')) if v.get('Name')=='C_add']
    if len(variables)!=1:raise ValueError('Expected one C_add variable')
    variables[0].set('Value',str(caps[0]/1000))
    sweeps=tree.find(tag('Sweeps'))
    for q in list(sweeps):sweeps.remove(q)
    fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
    ET.SubElement(fs,tag('Sweep'),Type='ABS_ENTRY',On='TRUE',Adaptive='TRUE',Start=str(cal_start_GHz),Stop=str(cal_stop_GHz),Target=str(samples-1))
    for q in tree.iter(tag('File')):
        if q.get('Type')=='TOUCH':q.set('Name',project.with_suffix('.snp').name);q.set('UseInterpolatedData','TRUE')
    write_xml(tree,project)
    source_tree=ET.parse(source)
    metal=lambda t:[ET.tostring(p) for p in t.iter(tag('PlanarPolygon'))]
    if metal(tree)!=metal(source_tree):raise ValueError('Explicit calibration changed metal')
    manifest.write_text(json.dumps(dict(configuration=expected,metal_identical=True,
        copied_existing_cache=(native/'sondata').exists()),indent=2)+'\n')
    start=time.perf_counter();started=datetime.now(timezone.utc).isoformat()
    result=legacy.main(fit_capacitance=False,project_file=project,output_dir=out,
        c_add_sweep_pf=[c/1000 for c in caps],target_l0_nh=target_L_nH,reuse_cache_after_first=True)
    rows=[];files={str(source.resolve()):digest(source),str(project.resolve()):digest(project),
                   str(project.with_suffix('.snp').resolve()):digest(project.with_suffix('.snp'))}
    for cap in caps:
        stem=f'cadd_{cap:g}ff';log=out/(stem+'.log');text=log.read_text()
        if 'Em warnings detected: 0' not in text:raise ValueError('Explicit calibration has native warnings')
        counts=[int(v) for v in re.findall(r'Total time for (\d+) frequencies:',text)]
        resource=json.loads((out/(stem+'.log.resources.json')).read_text())
        rows.append(dict(C_add_fF=cap,new_EM_frequency_counts=counts,solver_process_seconds=resource['elapsed_seconds']))
        for suffix in ('.sonx','.snp'):files[str((out/(stem+suffix)).resolve())]=digest(out/(stem+suffix))
    elapsed=time.perf_counter()-start
    result.update(case=name,method=expected['method'],target_error_percent=result['target_difference_percent'],
        source_project=str(source.resolve()),explicit_project=str(project.resolve()),
        capacitor_connection='Native Sonnet CAP',resonance_source='Native Sonnet exported S21 sampled maxima',
        Python_frequency_approximation=False,started_utc=started,total_extraction_seconds=elapsed,timings=rows)
    if result['R_squared']<.999:raise ValueError('Explicit capacitor regression is not sufficiently linear')
    if any(f<=cal_start_GHz*1e9+2e5 or f>=cal_stop_GHz*1e9-2e5 for f in [r['frequency_Hz'] for r in result['sweep']]):
        raise ValueError('Explicit calibration resonance is at the search boundary')
    if digest(source)!=expected['source_sha256']:raise ValueError('Original calibration fixture changed')
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    marker.write_text(json.dumps(dict(configuration=expected,files_sha256=files),indent=2)+'\n')
    print('EXPLICIT L',name,result['L0_nH'],'nH;',elapsed,'seconds',flush=True)
    return result


def check_t1():
    name='plain_cw344_w220_p532_g56_j66_c590_tg6'
    source=ROOT/'sonnet_files/taper_exploration/final_designs'/('new_T1_'+name+'_final')/'Lcal/generated_ncell_1.sonx'
    saved=ROOT/'output_data/taper_exploration/three_tapers_v3_2026-10-26'/name/'Lcal/results.json'
    caps=[r['C_add_fF'] for r in json.loads(saved.read_text())['sweep']]
    return measure(source,ROOT/'sonnet_files/taper_exploration/explicit_L_checks_2026-10-28/T1',
                   T1_CHECK.parent,name='T1',target_L_nH=5.97,caps_fF=caps,samples=16001)


def t1_inductance():
    q=check_t1()
    if q.get('Python_frequency_approximation') is not False or abs(q['target_error_percent'])>2:
        raise ValueError('T1 explicit inductance check must pass')
    return q['L0_nH']


if __name__=='__main__':check_t1()
