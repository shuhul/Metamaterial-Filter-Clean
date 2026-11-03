"""Equal-L ideal control of the new direct-T1-output fixture, with open check."""
from copy import deepcopy
import json
import xml.etree.ElementTree as ET
import numpy as np
import gdstk
import skrf
from geometry.taper.progressive_fixtures import ROOT,NATIVE,OUT,CASES,run
from geometry.taper.response_cases import W
from geometry.taper.ideal_wire_control import load_inductor
from geometry.prepare_sonnet import tag,points,write_xml
from geometry.taper.run_bridged_stage import digest


def prepare():
    source=NATIVE/'02a_T2_absent_T1/generated_ncell_2.sonx'
    folder=NATIVE/'02b_T2_ideal_T1';folder.mkdir(parents=True,exist_ok=True)
    out=OUT/'02b_T2_ideal_T1';out.mkdir(parents=True,exist_ok=True)
    project=folder/'wire_terminal_control.sonx'
    if not project.exists():
        tree=ET.parse(source);root=tree.getroot()
        ref=ET.parse(ROOT/'sonnet_files/taper_exploration/ideal_wire_diagnostic/wire_terminal_control.sonx')
        comp=deepcopy(next(ref.iter(tag('Component'))))
        comp.set('Id',str(max(int(q.get('Id')) for q in root.iter(tag('PlanarPolygon')))+1))
        metal=[gdstk.Polygon(points(q)) for q in root.iter(tag('PlanarPolygon'))]
        if gdstk.inside([(836,465),(836,467),(836,1065),(836,1067)],metal)!=(True,False,False,True):
            raise ValueError('Wrong ideal-control wire terminals')
        level=next(q for q in root.find(tag('Geometry')).findall(tag('Level')) if q.get('Number')=='0')
        level.find(tag('Objects')).append(comp)
        for q in root.iter(tag('File')):
            if q.get('Type')=='TOUCH':q.set('Name',project.with_suffix('.snp').name)
        write_xml(tree,project)
    return project,out,source


def main():
    source,_=run('02a_T2_absent_T1',*CASES['02a_T2_absent_T1'])
    project,out,_=prepare();marker=out/'operator_completed.json'
    if marker.exists():
        m=json.loads(marker.read_text())
        if m['project_sha256']!=digest(project) or m['snp_sha256']!=digest(project.with_suffix('.snp')):
            raise ValueError('Ideal operator changed')
        net=W['read_network'](project.with_suffix('.snp'),4)
    else:
        net=W['em'](project,4)
        marker.write_text(json.dumps(dict(project_sha256=digest(project),snp_sha256=digest(project.with_suffix('.snp'))),indent=2)+'\n')
    reference=W['read_network'](source,2)
    if not np.allclose(reference.f,net.f,rtol=0,atol=1):raise ValueError('Open-control grids differ')
    opened=net.s[:,:2,:2]+net.s[:,:2,2:]@np.linalg.solve(np.eye(2)[None]-net.s[:,2:,2:],net.s[:,2:,:2])
    delta=abs(opened-reference.s)
    r=dict(S21_RMS=float(np.sqrt(np.mean(delta[:,1,0]**2))),S21_max=float(delta[:,1,0].max()))
    r['passed']=r['S21_RMS']<.005 and r['S21_max']<.05
    (out/'open_control.json').write_text(json.dumps(r,indent=2)+'\n')
    if not r['passed']:raise ValueError('Open control perturbs output fixture')
    from geometry.taper.progressive_fitting import fit
    from geometry.taper.staged_fitting import LT1
    baseline=json.loads((OUT/'01_T1_regular/twoC/fit.json').read_text())
    for l in (.7327695472953403,2.5200619105580704):
        path=out/f'ideal_{l:.6f}.s2p'
        skrf.Network(frequency=skrf.Frequency.from_f(net.f,unit='hz'),s=load_inductor(net.s,net.f,l),z0=50).write_touchstone(str(path))
        fit(path,out/f'L{l:.6f}',[l,LT1],[0,baseline['C_fF'][0]],[None,baseline['G_fF'][0]],
            [('C',1),('G',1)],title='Ideal T2: Direct T1 Output-Claw Control',label='EM fixture + ideal L')


if __name__=='__main__':main()
