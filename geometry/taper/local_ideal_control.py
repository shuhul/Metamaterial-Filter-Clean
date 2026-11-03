"""Ideal shunt-L control at the short, original node/ground calibration gap.

The physical wire is removed. No 600 um component-port span is present.
An isolated four-port check calibrates L with the same four-capacitor method;
the full operator's open response is checked against saved ordinary EM fields.
"""
from copy import deepcopy
import json
import shutil
import xml.etree.ElementTree as ET
import numpy as np
import gdstk
import skrf
from scipy.optimize import brentq
from scipy.signal import find_peaks
from geometry.prepare_sonnet import tag,points,write_xml,BIN,run as native_command
from geometry.taper.corrected_iteration import ROOT,OUT,review
from geometry.taper.calibration_reference_control import measured_l
from geometry.taper.measure_compact_inductance import W
from geometry.taper.sonnet import configure
from geometry.taper.run_bridged_stage import digest
from geometry.taper.ideal_wire_control import load_inductor
from geometry.taper.second_layout import from_dict
from geometry.taper.compact_stages import prepare as prepare_geometry


def prepare(name):
    native=ROOT/'sonnet_files/taper_exploration/compact_centered'/('corrected_'+name)
    folder=native/'local_ideal_control';folder.mkdir(exist_ok=True)
    out=OUT/name/'local_ideal';out.mkdir(parents=True,exist_ok=True)
    physical=native/'idc/generated_ncell_3.sonx'
    absent=native.parent/('corrected_'+name+'_wire_absent')/'idc/generated_ncell_3.sonx'
    metadata=json.loads((ROOT/'geometry/taper/output/compact_centered'/('corrected_'+name)/'idc/geometry.json').read_text())
    if not absent.exists():
        prepare_geometry('corrected_'+name+'_wire_absent',from_dict(metadata['second']),'idc',omit_wire=True)
    source=native/'idc_Lcal/generated_ncell_1.sonx'
    isolated=folder/'isolated_local_operator.sonx';full=folder/'full_local_operator.sonx'
    if not isolated.exists():
        tree=ET.parse(source);root=tree.getroot()
        x=metadata['cell_centers'][1];width=metadata['second']['taper']['meander']['t']
        opening=[gdstk.rectangle((x-width/2,466),(x+width/2,1066))]
        wires=[q for q in root.iter(tag('PlanarPolygon'))
               if not gdstk.boolean([gdstk.Polygon(points(q))],opening,'xor',precision=.001)]
        if len(wires)!=1:raise ValueError('Expected exactly one separate calibration wire')
        wire=wires[0];next(q for q in root.iter() if wire in list(q)).remove(wire)
        for q in root.iter(tag('File')):
            if q.get('Type')=='TOUCH':q.set('Name',isolated.with_suffix('.snp').name)
        write_xml(tree,isolated);native_command([BIN/'soncmd.exe','-Read',isolated],folder/'isolated_read.log')
    if not full.exists():
        tree=ET.parse(absent);root=tree.getroot();prototype=ET.parse(isolated)
        comp=deepcopy(next(prototype.iter(tag('Component'))))
        comp.set('Id',str(max(int(q.get('Id')) for q in root.iter(tag('PlanarPolygon')))+100))
        comp.set('MacroID','shortGapIdealL');comp.set('Name','ShortGapIdealL')
        level=next(q for q in root.find(tag('Geometry')).findall(tag('Level')) if q.get('Number')=='0')
        level.find(tag('Objects')).append(comp)
        for q in root.iter(tag('File')):
            if q.get('Type')=='TOUCH':q.set('Name',full.with_suffix('.snp').name)
        sweeps=root.find(tag('Sweeps'))
        for q in list(sweeps):sweeps.remove(q)
        fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
        ET.SubElement(fs,tag('Sweep'),Type='ABS_ENTRY',On='TRUE',Adaptive='TRUE',Start='4.65',Stop='5.4',Target='7500')
        for q in root.iter(tag('File')):
            if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','TRUE')
        write_xml(tree,full);configure(full,4.65,5.4,7500)
        native_command([BIN/'soncmd.exe','-Read',full],folder/'full_read.log')
    # Repair the first endpoint-only diagnostic export without discarding its
    # valid field cache. It inherited three direct sweeps from the ordinary
    # control; configure() changes ranges, but does not change sweep type.
    tree=ET.parse(full);root=tree.getroot();sweeps=list(root.iter(tag('Sweep')))
    if len(sweeps)!=1 or sweeps[0].get('Type')!='ABS_ENTRY':
        archive=out/'endpoint_only_attempt';archive.mkdir(exist_ok=True)
        for source_file in (full,full.with_suffix('.snp'),full.parent/'em.log',out/'full_completed.json'):
            if source_file.exists() and not (archive/source_file.name).exists():shutil.copy2(source_file,archive/source_file.name)
        marker=out/'full_completed.json'
        if marker.exists():marker.unlink()
        blocks=root.find(tag('Sweeps'))
        for q in list(blocks):blocks.remove(q)
        fs=ET.SubElement(ET.SubElement(blocks,tag('Set'),On='TRUE'),tag('Frequencies'))
        ET.SubElement(fs,tag('Sweep'),Type='ABS_ENTRY',On='TRUE',Adaptive='TRUE',Start='4.65',Stop='5.4',Target='7500')
        for q in root.iter(tag('File')):
            if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','TRUE')
        write_xml(tree,full)
    # Imported metal must be exactly the already checked wire-absent layout.
    metal=lambda p:[gdstk.Polygon(points(q)) for q in ET.parse(p).iter(tag('PlanarPolygon'))]
    if gdstk.boolean(metal(full),metal(absent),'xor',precision=.001):raise ValueError('Local ideal control changes metal')
    x=metadata['cell_centers'][1];width=metadata['second']['taper']['meander']['t']
    cut=[gdstk.rectangle((x-width/2,466),(x+width/2,1066))]
    if gdstk.boolean(gdstk.boolean(metal(physical),cut,'not',precision=.001),metal(absent),'xor',precision=.001):
        raise ValueError('Wire-absent reference changes metal beyond the physical wire')
    components=list(ET.parse(full).iter(tag('Component')))
    if len(components)!=1:raise ValueError('Unexpected full-control components')
    for pin in components[0].findall(tag('Pin')):
        if pin.findtext(tag('Center')).split(',')[-1].rstrip(')') not in ('448','456'):
            raise ValueError('Ideal control is not on the short calibration gap')
    return isolated,full,absent,out,physical


def cached_em(project,out,kind):
    marker=out/(kind+'_completed.json')
    if marker.exists():
        saved=json.loads(marker.read_text())
        if saved!=dict(project_sha256=digest(project),snp_sha256=digest(project.with_suffix('.snp'))):
            raise ValueError('Changed completed short-gap operator')
        return W['read_network'](project.with_suffix('.snp'),4)
    net=W['em'](project,4)
    marker.write_text(json.dumps(dict(project_sha256=digest(project),snp_sha256=digest(project.with_suffix('.snp'))),indent=2)+'\n')
    return net


def ordinary_reference(name,absent,out,net,opened):
    """Reuse checked ordinary fields, or calculate only the two open peaks."""
    marker=out/'ordinary_completed.json'
    legacy=OUT/name/'paired_ideal/open_control_completed.json'
    existing=marker if marker.exists() else legacy if legacy.exists() else None
    if existing:
        saved=json.loads(existing.read_text())
        if saved!=dict(project_sha256=digest(absent),snp_sha256=digest(absent.with_suffix('.snp'))):
            raise ValueError('Ordinary wire-absent reference changed')
        return W['read_network'](absent.with_suffix('.snp'),2),True
    peaks,_=find_peaks(abs(opened[:,1,0]),prominence=.01)
    indices=sorted(sorted(peaks,key=lambda i:abs(opened[i,1,0]),reverse=True)[:2])
    if len(indices)!=2:raise ValueError('Expected two resolved open-reference peaks')
    tree=ET.parse(absent);root=tree.getroot();sweeps=root.find(tag('Sweeps'))
    for q in list(sweeps):sweeps.remove(q)
    fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
    for i in indices:
        f=str(net.f[i]/1e9)
        ET.SubElement(fs,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=f,Stop=f,Step='1')
    for q in root.iter(tag('File')):
        if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','FALSE')
    write_xml(tree,absent);ref=W['em'](absent,2)
    marker.write_text(json.dumps(dict(project_sha256=digest(absent),snp_sha256=digest(absent.with_suffix('.snp'))),indent=2)+'\n')
    return ref,False


def run(name,*,screening_only=False):
    isolated,full,absent,out,physical=prepare(name)
    cal=cached_em(isolated,out,'isolation')
    physical_cal=ROOT/'output_data/taper_exploration/compact_centered'/('corrected_'+name)/'idc_Lcal/results.json'
    l=json.loads(physical_cal.read_text())['L0_nH']
    component=brentq(lambda q:measured_l(cal.s,cal.f,q,local_gap=True)['measured_loop_L_nH']-l,.5,1.1,xtol=1e-8)
    check=measured_l(cal.s,cal.f,component,local_gap=True)
    if check['R_squared']<.9999:raise ValueError('Short-gap isolation calibration failed')
    (out/'calibration.json').write_text(json.dumps(dict(physical_measured_L_nH=l,component_L_nH=component,
        calibrated_check=check,physical_calibration=str(physical_cal),physical_calibration_sha256=digest(physical_cal),
        component_gap_um=8,long_component_span_present=False),indent=2)+'\n')
    print('SHORT GAP CALIBRATION',check,flush=True)
    net=cached_em(full,out,'full')
    opened=net.s[:,:2,:2]+net.s[:,:2,2:]@np.linalg.solve(np.eye(2)[None]-net.s[:,2:,2:],net.s[:,2:,:2])
    ref,reused=ordinary_reference(name,absent,out,net,opened)
    ix=np.array([np.argmin(abs(net.f-f)) for f in ref.f])
    if not np.allclose(net.f[ix],ref.f,rtol=0,atol=1):raise ValueError('Open validation frequency mismatch')
    error=abs(opened[ix]-ref.s)
    validation=dict(frequencies_GHz=(ref.f/1e9).tolist(),S21_RMS=float(np.sqrt(np.mean(error[:,1,0]**2))),
        S21_max=float(error[:,1,0].max()),all_S_max=float(error.max()),ordinary_fields_reused=reused)
    validation['passed']=bool(validation['S21_RMS']<.005 and validation['S21_max']<.05)
    (out/'open_control.json').write_text(json.dumps(validation,indent=2)+'\n')
    if not validation['passed'] and not (screening_only and validation['S21_RMS']<.01 and validation['S21_max']<.02):
        raise ValueError(f'Short-gap open check failed: {validation}')
    path=out/'short_gap_ideal_L.s2p'
    skrf.Network(frequency=skrf.Frequency.from_f(net.f,unit='hz'),s=load_inductor(net.s,net.f,component),z0=50).write_touchstone(str(path))
    label=('Short-gap ideal L (qualitative control)' if not validation['passed']
           else 'EM geometry + short-gap calibrated ideal inductor')
    result=review(path,out/'fit',l,data_label=label)
    comparison=dict(physical=json.loads((OUT/name/'fit/review.json').read_text()),ideal=result,
        effective_L_matched_in_isolation=True,component_L_nH=component,measured_L_nH=l,
        component_gap_um=8,calibration=check,open_control=validation,physical_project_sha256=digest(physical),accepted=False,
        control_precision_passed=validation['passed'],screening_only=not validation['passed'],
        limitation='A marginal open-control check may be used only for qualitative screening, never 2% parameter certification.' if not validation['passed'] else 'Diagnostic only; physical acceptance remains separate.')
    (out/'comparison.json').write_text(json.dumps(comparison,indent=2)+'\n')
    plot_comparison(name,comparison)
    return comparison


def plot_comparison(name,comparison):
    """Keep the physical fit, ideal-control fit and target visually distinct."""
    import matplotlib.pyplot as plt
    from geometry.taper.progressive_fitting import transmission
    from geometry.taper.corrected_iteration import REGULAR_C,REGULAR_G
    from geometry.taper.staged_fitting import LT1,LR
    base=OUT/name;fig,axes=plt.subplots(1,2,figsize=(18,7),sharex=True,sharey=True)
    panels=((axes[0],base/'fit',comparison['physical'],'Physical straight wire'),
            (axes[1],base/'local_ideal/fit',comparison['ideal'],
             'Short-gap ideal control'+(' (qualitative)' if comparison['screening_only'] else '')))
    for ax,folder,r,title in panels:
        data=np.load(folder/'response.npz');f=data['f']
        target=transmission(f,[comparison['measured_L_nH'],LT1,LR,LR],
            [100,100,REGULAR_C,REGULAR_C],[459,56.7,8.2865,REGULAR_G],0)
        ax.plot(f/1e9,abs(data['actual'][:,1,0]),color='#b22222',lw=2.6,label='EM response')
        ax.plot(f/1e9,abs(data['model'][:,1,0]),'k--',lw=2,label='Free-capacitance fit')
        ax.plot(f/1e9,abs(target),color='#1565c0',ls=':',lw=2,label='Target C values; fixed measured L')
        ax.set_title(title,fontsize=18);ax.set_xlabel('Frequency (GHz)',fontsize=18)
        ax.set_ylim(0,1.35);ax.tick_params(labelsize=14);ax.grid(alpha=.2);ax.legend(fontsize=11,loc='upper right')
        ax.text(.025,.96,f"Fit complex RMS = {r['complex_RMS']:.5f}",transform=ax.transAxes,va='top',fontsize=13)
    axes[0].set_ylabel(r'$|S_{21}|$',fontsize=24)
    caption=f"Measured L = {comparison['measured_L_nH']:.5f} nH; ideal component = {comparison['component_L_nH']:.5f} nH"
    if comparison['screening_only']:
        caption+=f"\nControl open-check max error = {comparison['open_control']['S21_max']:.4f}; strict RMS check not passed"
    fig.suptitle('Taper 2: Physical Wire and Short-Gap Ideal Control\n'+caption,fontsize=18)
    fig.tight_layout();fig.savefig(base/'local_ideal/physical_vs_local_ideal.png',dpi=160);plt.close(fig)
