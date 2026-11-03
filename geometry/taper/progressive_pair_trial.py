"""Matched ideal/physical full-assembly comparison for a selected T2 trial.

The ideal control removes only the straight wire, retains all other metal,
and replaces it with the independently measured L. Open validation precedes
interpretation; it never substitutes an ideal component for the final design.
"""
from copy import deepcopy
import json
import xml.etree.ElementTree as ET
import numpy as np
import gdstk
import skrf
from geometry.taper.progressive_fixtures import ROOT,NATIVE,GDS,OUT,run as run_output
from geometry.taper.response_cases import W
from geometry.taper.ideal_wire_control import load_inductor
from geometry.taper.progressive_fitting import fit
from geometry.taper.staged_fitting import LT1,LR
from geometry.taper.run_bridged_stage import digest
from geometry.prepare_sonnet import tag,points,write_xml


def run(name,*,bare=False):
    meta=json.loads((GDS/name/'geometry.json').read_text())
    source=GDS/name/'source/geometry.gds'
    if not source.exists():source=ROOT/meta['source_gds']
    absent,_=run_output(name+'/ideal_absent',str(source),3,302,True)
    physical=NATIVE/name/'generated_ncell_3.sonx'
    x=meta['cell_centers'][1];wire=meta['source_metadata']['second']['taper']['meander']['t']
    # Check the native metal change against the intended wire cut, after import.
    geometry=lambda p:[gdstk.Polygon(points(q)) for q in ET.parse(p).iter(tag('PlanarPolygon'))]
    expected=gdstk.boolean(geometry(physical),[gdstk.rectangle((x-wire/2,466),(x+wire/2,1066))],
                           'not',precision=.001)
    if gdstk.boolean(expected,geometry(absent.with_suffix('.sonx')),'xor',precision=.001):
        raise ValueError('Ideal fixture changes metal beyond the wire')
    folder=NATIVE/name/'ideal_control';folder.mkdir(parents=True,exist_ok=True)
    out=OUT/name/'paired_ideal';out.mkdir(parents=True,exist_ok=True)
    project=folder/'ideal_wire.sonx'
    if not project.exists():
        tree=ET.parse(absent.with_suffix('.sonx'));root=tree.getroot()
        ref=ET.parse(ROOT/'sonnet_files/taper_exploration/ideal_wire_diagnostic/wire_terminal_control.sonx')
        comp=deepcopy(next(ref.iter(tag('Component'))))
        comp.set('Id',str(max(int(q.get('Id')) for q in root.iter(tag('PlanarPolygon')))+1))
        comp.find(tag('LabelPos')).text=f'({x+6},760)'
        comp.find(tag('SymbolBox')).text=f'({x-2},470,{x+2},1062)'
        for pin,y in zip(comp.findall(tag('Pin')),(466,1066)):pin.find(tag('Center')).text=f'({x},{y})'
        level=next(q for q in root.find(tag('Geometry')).findall(tag('Level')) if q.get('Number')=='0')
        level.find(tag('Objects')).append(comp)
        for q in root.iter(tag('File')):
            if q.get('Type')=='TOUCH':q.set('Name',project.with_suffix('.snp').name)
        write_xml(tree,project)
    marker=out/'operator_completed.json'
    if marker.exists():
        m=json.loads(marker.read_text())
        if m['project_sha256']!=digest(project) or m['snp_sha256']!=digest(project.with_suffix('.snp')):
            raise ValueError('Changed matched ideal operator')
        net=W['read_network'](project.with_suffix('.snp'),4)
    else:
        net=W['em'](project,4)
        marker.write_text(json.dumps(dict(project_sha256=digest(project),snp_sha256=digest(project.with_suffix('.snp'))),indent=2)+'\n')
    reference=W['read_network'](absent,2)
    if not np.allclose(net.f,reference.f,rtol=0,atol=1):raise ValueError('Ideal open-control grids differ')
    opened=net.s[:,:2,:2]+net.s[:,:2,2:]@np.linalg.solve(np.eye(2)[None]-net.s[:,2:,2:],net.s[:,2:,:2])
    e=abs(opened-reference.s)
    control=dict(S21_RMS=float(np.sqrt(np.mean(e[:,1,0]**2))),S21_max=float(e[:,1,0].max()))
    control['passed']=control['S21_RMS']<.005 and control['S21_max']<.05
    (out/'open_control.json').write_text(json.dumps(control,indent=2)+'\n')
    if not control['passed']:raise ValueError('Matched ideal open control failed')
    actual_folder=OUT/name/'junction_check' if bare else OUT/name/'twoC'
    actual_fit=json.loads((actual_folder/'fit.json').read_text());l=actual_fit['L_nH'][0]
    path=out/'ideal_loaded.s2p'
    skrf.Network(frequency=skrf.Frequency.from_f(net.f,unit='hz'),s=load_inductor(net.s,net.f,l),z0=50).write_touchstone(str(path))
    baseline=json.loads((OUT/'01_T1_regular/twoC/fit.json').read_text())
    r=fit(path,out/'twoC',[l,LT1,LR],[0 if bare else 100,baseline['C_fF'][0],190.4],
          [None if bare else 459,56.7 if bare else baseline['G_fF'][0],8.2865],
          [('G',1)] if bare else [('C',0),('G',0)],
          title='Matched Ideal T2 Control',label='EM geometry + ideal inductor')
    physical_data=np.load(actual_folder/'response.npz');ideal_data=np.load(out/'twoC/response.npz')
    if not np.allclose(physical_data['f'],ideal_data['f'],rtol=0,atol=1):raise ValueError('Paired response grids differ')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(15,6),sharex=True,sharey=True)
    for ax,data,result,label in ((axes[0],physical_data,actual_fit,'Physical straight wire'),
                                  (axes[1],ideal_data,r,'Ideal inductor, same EM fixture')):
        f=data['f'];a=data['actual'][:,1,0];m=data['model'][:,1,0]
        ax.plot(f/1e9,abs(a),color='#b22222',lw=2.7,label=label)
        ax.plot(f/1e9,abs(m),'k--',lw=2,label='Same two-C fit, fixed downstream values')
        ax.set_title(f"{label}\nComplex RMS = {result['complex_RMS']:.5f}",fontsize=16)
        ax.set_xlabel('Frequency (GHz)',fontsize=16);ax.tick_params(labelsize=13);ax.grid(alpha=.25);ax.legend(fontsize=10)
    axes[0].set_ylabel(r'$|S_{21}|$',fontsize=22)
    fig.suptitle(f'Matched Ideal / Physical T2: L = {l:.5f} nH',fontsize=22)
    fig.tight_layout();fig.savefig(out/'paired_s21.png',dpi=190);plt.close(fig)
    summary=dict(physical=actual_fit,ideal=r,open_control=control,
                 direct_physical_ideal_S21_RMS=float(np.sqrt(np.mean(abs(physical_data['actual'][:,1,0]-ideal_data['actual'][:,1,0])**2))),
                 only_physical_wire_replaced=True,diagnostic_only=True)
    (out/'comparison.json').write_text(json.dumps(summary,indent=2)+'\n')
    print('PAIRED TRIAL',name,actual_fit['complex_RMS'],r['complex_RMS'],flush=True)
    return summary


if __name__=='__main__':
    import sys
    run(sys.argv[1],bare='--bare' in sys.argv)
