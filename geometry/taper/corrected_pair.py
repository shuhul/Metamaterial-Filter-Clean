"""Historical nominal-component-L diagnostic; NOT matched effective L.

The long component span was subsequently found to change the isolation
measurement substantially. Use local_ideal_control for new comparisons.
This historical control removes only the straight wire and loads its original
terminals with a nominal component L. An ordinary wire-absent EM
calculation checks the four-port operator at its predicted open-response peaks.
The physical project remains unchanged. Controls are never final geometries.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from copy import deepcopy
from pathlib import Path
import json
import xml.etree.ElementTree as ET
import numpy as np
import gdstk
import skrf
from scipy.signal import find_peaks
from geometry.taper.corrected_iteration import ROOT,OUT,review
from geometry.taper.compact_stages import prepare
from geometry.taper.corrected_design import pad_parameters
from geometry.taper.ideal_wire_control import load_inductor
from geometry.taper.response_cases import W
from geometry.taper.sonnet import configure
from geometry.prepare_sonnet import tag,points,write_xml
from geometry.taper.run_bridged_stage import digest


def prepare_operator(name,p):
    absent,_,meta=prepare('corrected_'+name+'_wire_absent',p,'idc',omit_wire=True)
    physical=absent.parents[2]/('corrected_'+name)/'idc/generated_ncell_3.sonx'
    metal=lambda path:[gdstk.Polygon(points(q)) for q in ET.parse(path).iter(tag('PlanarPolygon'))]
    x=meta['cell_centers'][1];w=p.taper.meander.t
    cut=[gdstk.rectangle((x-w/2,466),(x+w/2,1066))]
    expected=gdstk.boolean(metal(physical),cut,'not',precision=.001)
    if gdstk.boolean(expected,metal(absent),'xor',precision=.001):
        raise ValueError('Ideal fixture changes metal beyond the straight wire')
    folder=physical.parent.parent/'matched_ideal_control';folder.mkdir(exist_ok=True)
    project=folder/'ideal_wire_operator.sonx'
    out=OUT/name/'paired_ideal';out.mkdir(parents=True,exist_ok=True)
    if not project.exists():
        tree=ET.parse(absent);r=tree.getroot();actual=metal(absent)
        if gdstk.inside([(x,465),(x,467),(x,1065),(x,1067)],actual)!=(True,False,False,True):
            raise ValueError('Ideal terminals are not on the intended bar and bottom-ground edges')
        ref=ET.parse(ROOT/'sonnet_files/taper_exploration/ideal_wire_diagnostic/wire_terminal_control.sonx')
        component=deepcopy(next(ref.iter(tag('Component'))))
        component.set('Id',str(max(int(q.get('Id')) for q in r.iter(tag('PlanarPolygon')))+1))
        component.find(tag('LabelPos')).text=f'({x+6},760)'
        component.find(tag('SymbolBox')).text=f'({x-2},470,{x+2},1062)'
        for pin,y in zip(component.findall(tag('Pin')),(466,1066)):
            pin.find(tag('Center')).text=f'({x},{y})'
        level=next(q for q in r.find(tag('Geometry')).findall(tag('Level')) if q.get('Number')=='0')
        level.find(tag('Objects')).append(component)
        write_xml(tree,project);configure(project,4.65,5.4,7500)
    (out/'geometry_check.json').write_text(json.dumps(dict(physical=str(physical),absent=str(absent),
        native_wire_cut_XOR_empty=True,only_wire_removed=True,T1_unchanged=True,
        physical_sha256=digest(physical),absent_sha256=digest(absent)),indent=2)+'\n')
    (folder/'DIAGNOSTIC_ONLY.md').write_text('Wire removed and replaced numerically with an equal-L ideal inductor. Not a physical or final design.\n')
    return project,absent,out


def run(name,p,l):
    project,absent,out=prepare_operator(name,p)
    marker=out/'operator_completed.json'
    if marker.exists():
        rec=json.loads(marker.read_text())
        if rec['project_sha256']!=digest(project) or rec['snp_sha256']!=digest(project.with_suffix('.snp')):
            raise ValueError('Changed completed ideal operator')
        net=W['read_network'](project.with_suffix('.snp'),4)
    else:
        net=W['em'](project,4)
        marker.write_text(json.dumps(dict(project_sha256=digest(project),snp_sha256=digest(project.with_suffix('.snp'))),indent=2)+'\n')
    opened=net.s[:,:2,:2]+net.s[:,:2,2:]@np.linalg.solve(np.eye(2)[None]-net.s[:,2:,2:],net.s[:,2:,:2])
    peaks,_=find_peaks(abs(opened[:,1,0]),prominence=.01)
    strongest=sorted(peaks,key=lambda i:abs(opened[i,1,0]),reverse=True)[:2]
    # The two peak checks are the sensitive tests for port-induced resonance
    # shifts. Avoid an additional low-information midband EM field per trial.
    indices=sorted(set(strongest))
    if len(indices)<2:indices=sorted(set(indices+[len(net.f)//2,len(net.f)//4]))[:2]
    frequencies=(net.f[indices]/1e9).tolist()
    check_marker=out/'open_control_completed.json'
    if not check_marker.exists():
        tree=ET.parse(absent);sweeps=tree.find(tag('Sweeps'))
        for child in list(sweeps):sweeps.remove(child)
        fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
        for f in frequencies:
            ET.SubElement(fs,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(f),Stop=str(f),Step='1')
        for q in tree.iter(tag('File')):
            if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','FALSE')
        write_xml(tree,absent);ref=W['em'](absent,2)
        check_marker.write_text(json.dumps(dict(project_sha256=digest(absent),snp_sha256=digest(absent.with_suffix('.snp'))),indent=2)+'\n')
    else:
        rec=json.loads(check_marker.read_text())
        if rec['project_sha256']!=digest(absent) or rec['snp_sha256']!=digest(absent.with_suffix('.snp')):
            raise ValueError('Changed ordinary open control')
        ref=W['read_network'](absent.with_suffix('.snp'),2)
    ix=np.array([int(np.argmin(abs(net.f-f))) for f in ref.f])
    if not np.allclose(net.f[ix],ref.f,rtol=0,atol=1):raise ValueError('Open-control frequencies differ')
    delta=abs(opened[ix]-ref.s)
    check=dict(frequencies_GHz=(ref.f/1e9).tolist(),S21_RMS=float(np.sqrt(np.mean(delta[:,1,0]**2))),
               S21_max=float(delta[:,1,0].max()),all_S_max=float(delta.max()),
               includes_predicted_open_peaks=bool(strongest) and all(
                   np.min(abs(ref.f-net.f[i]))<1 for i in strongest))
    check['passed']=bool(check['includes_predicted_open_peaks'] and check['S21_RMS']<.005 and check['S21_max']<.05)
    (out/'open_control.json').write_text(json.dumps(check,indent=2)+'\n')
    if not check['passed']:raise ValueError(f'Ideal open validation failed: {check}')
    path=out/'ideal_same_L.s2p'
    skrf.Network(frequency=skrf.Frequency.from_f(net.f,unit='hz'),s=load_inductor(net.s,net.f,l),z0=50).write_touchstone(str(path))
    result=review(path,out/'fit',l,data_label='EM geometry + nominal-L long-span component')
    physical=json.loads((OUT/name/'fit/review.json').read_text())
    if not np.isclose(physical['fitted']['L_T2_nH'],l,rtol=0,atol=1e-12):
        raise ValueError('Ideal and physical comparisons must use exactly the same measured L')
    if [tuple(q) for q in physical['free_capacitances']]!=[tuple(q) for q in result['free_capacitances']]:
        raise ValueError('Ideal and physical fits used different free parameters')
    pdata=np.load(OUT/name/'fit/response.npz');idata=np.load(out/'fit/response.npz')
    if not np.allclose(pdata['f'],idata['f'],rtol=0,atol=1):
        raise ValueError('Ideal and physical comparison frequency grids differ')
    comparison=dict(physical=physical,ideal=result,inductance_nH=l,same_fit_policy=True,
                    effective_L_equivalence_verified=False,component_span_um=600,
                    scope='Nominal component-L diagnostic; use local_ideal_control for new comparisons.',
                    open_control=check,physical_wire_is_only_geometry_change=True,accepted=False,
                    direct_physical_ideal_S21_complex_RMS=float(np.sqrt(np.mean(
                        abs(pdata['actual'][:,1,0]-idata['actual'][:,1,0])**2))))
    (out/'comparison.json').write_text(json.dumps(comparison,indent=2)+'\n')
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(16,6.5),sharex=True,sharey=True)
    for ax,folder,r,label in ((axes[0],OUT/name/'fit',physical,'Physical straight wire'),
                              (axes[1],out/'fit',result,'Nominal-L long-span control')):
        d=np.load(folder/'response.npz');f=d['f']
        ax.plot(f/1e9,abs(d['actual'][:,1,0]),color='#b22222',lw=2.8,label=label)
        ax.plot(f/1e9,abs(d['model'][:,1,0]),'k--',lw=2.2,label='Same free-capacitance fit')
        ax.set_title(f"{label}\nComplex S21 RMS = {r['complex_RMS']:.5f}",fontsize=17)
        ax.set_xlabel('Frequency (GHz)',fontsize=17);ax.tick_params(labelsize=14);ax.grid(alpha=.2);ax.legend(fontsize=11)
    axes[0].set_ylabel(r'$|S_{21}|$',fontsize=23)
    fig.suptitle(f'Nominal Component-L Diagnostic: L = {l:.5f} nH',fontsize=23)
    fig.tight_layout();fig.savefig(out/'paired_s21_corrected.png',dpi=190);plt.close(fig)
    return comparison


if __name__=='__main__':
    import sys
    name='pads24_d12_wire44';p=pad_parameters()
    if '--prepare-only' in sys.argv:prepare_operator(name,p)
    else:
        cal=ROOT/'output_data/taper_exploration/compact_centered'/('corrected_'+name)/'idc_Lcal/results.json'
        run(name,p,json.loads(cal.read_text())['L0_nH'])
