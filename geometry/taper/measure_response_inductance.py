"""Fast isolated inductance calibration for an ordinary T2 response case."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import xml.etree.ElementTree as ET
from geometry.taper.response_cases import ROOT
from geometry.taper.second_layout import from_dict,write,validate
from geometry.taper.sonnet import add_capacitor
from geometry.prepare_sonnet import main as prepare,tag,write_xml
from geometry.taper.measure_compact_inductance import calibrate

def measure(name):
    meta=json.loads((ROOT/'geometry/taper/output'/f'{name}.json').read_text());p=from_dict(meta['second'])
    iso=name+'_L_cal';out=ROOT/'output_data/taper_exploration'/iso;out.mkdir(parents=True,exist_ok=True)
    project=ROOT/'sonnet_files/taper_exploration'/iso/'generated_ncell_1.sonx'
    if not project.exists():
        gds=write(p,ROOT/'geometry/taper/output'/iso,isolation=True)
        prepare(1,gds_path=gds,output_dir=project.parent,geometry_validator=validate)
        add_capacitor(project,json.loads(gds.with_suffix('.json').read_text()))
        tree=ET.parse(project);comp=next(tree.iter(tag('Component')))
        comp.set('Type','PORTS');comp.attrib.pop('Units',None)
        for q in list(comp):
            if q.tag==tag('Value'):comp.remove(q)
        sweeps=tree.find(tag('Sweeps'))
        for q in list(sweeps):sweeps.remove(q)
        block=ET.SubElement(sweeps,tag('Set'),On='TRUE');fs=ET.SubElement(block,tag('Frequencies'))
        for value in (4,5,5.5,6,7):
            ET.SubElement(fs,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(value),Stop=str(value),Step='1')
        for q in tree.findall(f"{tag('OutputFiles')}//{tag('File')}"):
            if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','FALSE')
        write_xml(tree,project)
    return calibrate(project,out,name,auto_loads=True)

if __name__=='__main__':
    import sys
    measure(sys.argv[1])
