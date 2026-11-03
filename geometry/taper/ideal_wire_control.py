"""Diagnostic only: remove the physical wire and load its terminals with L.

Open terminals must reproduce the independently simulated wire-absent layout.
An ideal load is not a manufacturable geometry and cannot be a final design.
"""
from copy import deepcopy
import json
import xml.etree.ElementTree as ET
import gdstk
import numpy as np
from circuit_common import read_snp
from geometry.prepare_sonnet import tag,points,write_xml
from geometry.taper.compact_stages import ROOT,NATIVE,GDS
from geometry.taper.response_cases import W
from geometry.taper.fit_restored_t2 import features


def load_inductor(s,f,l_nH):
    z=2j*np.pi*np.asarray(f)*l_nH*1e-9
    load=np.empty((len(f),2,2),complex)
    load[:,0,0]=load[:,1,1]=z/(z+100)
    load[:,0,1]=load[:,1,0]=100/(z+100)
    return s[:,:2,:2]+s[:,:2,2:]@load@np.linalg.solve(np.eye(2)[None]-s[:,2:,2:]@load,s[:,2:,:2])


def prepare():
    source=NATIVE/'w500_wire52_gap6_wire_absent/direct/generated_ncell_3.sonx'
    meta=json.loads((GDS/'w500_wire52_gap6_wire_absent/direct/geometry.json').read_text())
    folder=ROOT/'sonnet_files/taper_exploration/ideal_wire_diagnostic';folder.mkdir(exist_ok=True)
    out=ROOT/'output_data/taper_exploration/ideal_wire_diagnostic';out.mkdir(exist_ok=True)
    p=folder/'wire_terminal_control.sonx'
    if p.exists():return p,out,source
    tree=ET.parse(source);root=tree.getroot();x=meta['cell_centers'][1]
    metal=[gdstk.Polygon(points(q)) for q in root.iter(tag('PlanarPolygon'))]
    if gdstk.inside([(x,465),(x,467),(x,1065),(x,1067)],metal)!=(True,False,False,True):raise ValueError('Wire terminals are not exposed')
    prototype=next(ET.parse(ROOT/'sonnet_files/filter_v6.sonx').iter(tag('Component')))
    comp=deepcopy(prototype);comp.set('Type','PORTS');comp.set('Name','idealWireControl');comp.set('MacroID','idealWireControl')
    comp.set('Id',str(max(int(q.get('Id')) for q in root.iter(tag('PlanarPolygon')))+1));comp.attrib.pop('Units',None)
    for e in list(comp):
        if e.tag==tag('Value'):comp.remove(e)
    comp.find(tag('LabelPos')).text=f'({x+6},760)';comp.find(tag('SymbolBox')).text=f'({x-2},470,{x+2},1062)'
    for i,(pin,y,orientation) in enumerate(zip(comp.findall(tag('Pin')),(466,1066),('TOP','BOTTOM'))):
        pin.set('MacroID',f'wireControlPin{i}');pin.set('PortNum',str(3+i));pin.set('Orientation',orientation)
        pin.find(tag('Center')).text=f'({x},{y})'
    level=next(q for q in root.find(tag('Geometry')).findall(tag('Level')) if q.get('Number')=='0')
    level.find(tag('Objects')).append(comp)
    for q in root.iter(tag('File')):
        if q.get('Type')=='TOUCH':q.set('Name',p.with_suffix('.snp').name)
    write_xml(tree,p)
    (folder/'DIAGNOSTIC_ONLY.md').write_text('Physical wire removed. Four-port field data support an ideal 0.74 nH load in Python. This is not a physical/final T2 design.\n')
    return p,out,source


def run():
    p,out,source=prepare();n=W['em'](p,4)
    ref=read_snp(source.with_suffix('.snp'))
    if not np.allclose(n.f,ref.f,rtol=0,atol=1):raise ValueError('Control grids differ')
    opened=n.s[:,:2,:2]+n.s[:,:2,2:]@np.linalg.solve(np.eye(2)[None]-n.s[:,2:,2:],n.s[:,2:,:2])
    error=abs(opened-ref.s)
    r=dict(diagnostic_only=True,physical_design_accepted=False,
        open_S21_RMS=float(np.sqrt(np.mean(error[:,1,0]**2))),open_S21_max_error=float(error[:,1,0].max()),
        open_all_S_RMS=float(np.sqrt(np.mean(error**2))),
        control_passed=bool(np.sqrt(np.mean(error[:,1,0]**2))<.005 and error[:,1,0].max()<.05))
    (out/'open_control.json').write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2),flush=True)
    if not r['control_passed']:raise ValueError('Ideal-wire terminal control perturbs the reference; no interpretation')
    s=load_inductor(n.s,n.f,.74);r['loaded_peaks']=features(n.f,s)
    np.savez(out/'ideal_wire_response.npz',f=n.f,s=s,opened=opened)
    # Touchstone export explicitly labels this as a synthetic ideal-load result.
    import skrf
    skrf.Network(frequency=skrf.Frequency.from_f(n.f,unit='hz'),s=s,z0=50).write_touchstone(str(out/'diagnostic_ideal_wire.s2p'))
    from geometry.taper.compact_fit import fit_direct
    fit_direct(out/'diagnostic_ideal_wire.s2p',out/'fit',.74,data_label='EM layout + ideal wire (diagnostic)')
    (out/'results.json').write_text(json.dumps(r,indent=2)+'\n')
    return r


if __name__=='__main__':run()
