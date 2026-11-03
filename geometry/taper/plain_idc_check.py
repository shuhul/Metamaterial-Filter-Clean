"""Independent uniform-IDC check, with explicit low/microwave frequency data.

This smaller two-terminal fixture checks Cg and dispersion; its ground C does
not replace the full T2 C0 extraction. No inductor or T1 is drawn here.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import xml.etree.ElementTree as ET
import gdstk
import numpy as np
from geometry.taper.compact_stages import NATIVE,OUT,GDS
from geometry.taper.idc_screen import build,validate
from geometry.prepare_sonnet import main as prepare,tag,write_xml
from geometry.taper.response_cases import W
from geometry.taper.run_bridged_stage import digest


def run(name,p):
    dest=NATIVE/name/'idc_check';gd=GDS/name/'idc_check';out=OUT/name/'idc_check'
    for path in (dest,gd,out):path.mkdir(parents=True,exist_ok=True)
    top,meta=build(p);gds=gd/'geometry.gds';project=dest/'generated_ncell_1.sonx'
    if not project.exists():
        lib=gdstk.Library(unit=1e-6,precision=1e-9);lib.add(top);lib.write_gds(str(gds))
        gds.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
        prepare(1,gds_path=gds,output_dir=dest,geometry_validator=validate)
        tree=ET.parse(project);sweeps=tree.find(tag('Sweeps'))
        for q in list(sweeps):sweeps.remove(q)
        fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
        for f in (.1,.2,.3,4.8,5.2):
            ET.SubElement(fs,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(f),Stop=str(f),Step='1')
        for q in tree.findall(f"{tag('OutputFiles')}//{tag('File')}"):
            if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','FALSE')
        write_xml(tree,project)
    elif json.loads(gds.with_suffix('.json').read_text())!=meta:
        raise ValueError('Changed IDC geometry requires new case name')
    marker=out/'completed.json'
    if marker.exists():
        d=json.loads(marker.read_text())
        if digest(project)!=d['project_sha256'] or digest(project.with_suffix('.snp'))!=d['snp_sha256']:
            raise ValueError('IDC source changed')
        n=W['read_network'](project.with_suffix('.snp'),2)
    else:
        n=W['em'](project,2)
        marker.write_text(json.dumps(dict(project_sha256=digest(project),snp_sha256=digest(project.with_suffix('.snp'))),indent=2)+'\n')
    eye=np.eye(2)[None];y=(eye-n.s)@np.linalg.inv(eye+n.s)/50;w=2*np.pi*n.f
    mutual=-y[:,0,1].imag/w*1e15
    result=dict(source=str(project),IDC_pairs=p.idc_pairs,finger_width_um=p.idc_t,gap_um=p.idc_d,
        overlap_um=p.idc_span-20-2*p.idc_d,frequency_GHz=(n.f/1e9).tolist(),mutual_C_fF=mutual.tolist(),
        low_frequency_mean_C_fF=float(np.mean(mutual[n.f<1e9])),
        source_ground_C_fF=((y[:,0,0]+y[:,0,1]).imag/w*1e15).tolist(),
        node_ground_C_fF=((y[:,1,1]+y[:,1,0]).imag/w*1e15).tolist(),
        accepted_full_T2=False,scope='Independent IDC and dispersion check only; terminal surroundings differ from the final chain')
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True);return result
