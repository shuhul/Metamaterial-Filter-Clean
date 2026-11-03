"""Diagnostic: compare a wire-terminal ideal L with the isolation-measured L.

This never changes the physical design or its primary capacitance fit. The
six-port fixture retains the calibration probe and removes only the straight
wire, adding a second probe pair at the exposed wire terminals. Both ideal L
and known calibration capacitors are applied numerically to one EM operator.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
from copy import deepcopy
from pathlib import Path
import json
import xml.etree.ElementTree as ET
import numpy as np
import gdstk
from scipy.optimize import minimize_scalar, brentq
from scipy.stats import linregress
from geometry.prepare_sonnet import tag, points, write_xml, BIN, run as native_command
from geometry.taper.corrected_iteration import ROOT, OUT
from geometry.taper.measure_compact_inductance import AP, W
from geometry.taper.run_bridged_stage import digest


def prepare(name):
    source = ROOT/'sonnet_files/taper_exploration/compact_centered'/('corrected_'+name)/'idc_Lcal/generated_ncell_1.sonx'
    meta = json.loads((ROOT/'geometry/taper/output/compact_centered'/('corrected_'+name)/'idc_Lcal/geometry.json').read_text())
    folder = source.parent.parent/'calibration_reference_control'
    folder.mkdir(exist_ok=True)
    project = folder/'wire_terminal_L_calibration.sonx'
    out = OUT/name/'calibration_reference_control'
    out.mkdir(parents=True, exist_ok=True)
    if not project.exists():
        tree = ET.parse(source); root = tree.getroot()
        x = meta['cell_centers'][1]; width = meta['second']['taper']['meander']['t']
        opening = [gdstk.rectangle((x-width/2,466),(x+width/2,1066))]
        polygons = list(root.iter(tag('PlanarPolygon')))
        originals = [gdstk.Polygon(points(q)) for q in polygons]
        candidates = [q for q in polygons if not gdstk.boolean([gdstk.Polygon(points(q))], opening, 'xor', precision=.001)]
        if len(candidates) != 1:
            raise ValueError('Expected one separate straight-wire polygon')
        wire = candidates[0]
        next(q for q in root.iter() if wire in list(q)).remove(wire)
        metal = [gdstk.Polygon(points(q)) for q in root.iter(tag('PlanarPolygon'))]
        expected = gdstk.boolean(originals, opening, 'not', precision=.001)
        if gdstk.boolean(metal, expected, 'xor', precision=.001):
            raise ValueError('Calibration control changes more than the wire')
        if gdstk.inside([(x,465),(x,467),(x,1065),(x,1067)],metal) != (True,False,False,True):
            raise ValueError('Invalid exposed wire-terminal edges')
        comp = deepcopy(next(root.iter(tag('Component'))))
        comp.set('MacroID','idealWireCalibration'); comp.set('Name','IdealWireCalibration')
        comp.set('Id',str(max(int(q.get('Id')) for q in root.iter(tag('PlanarPolygon')))+100))
        comp.find(tag('LabelPos')).text = f'({x+6},760)'
        comp.find(tag('SymbolBox')).text = f'({x-2},470,{x+2},1062)'
        for i,(pin,y,orientation) in enumerate(zip(comp.findall(tag('Pin')),(466,1066),('TOP','BOTTOM'))):
            pin.set('MacroID',f'wireLCalPin{i}'); pin.set('PortNum',str(5+i))
            pin.set('Orientation',orientation); pin.find(tag('Center')).text=f'({x},{y})'
        level = next(q for q in root.find(tag('Geometry')).findall(tag('Level')) if q.get('Number')=='0')
        level.find(tag('Objects')).append(comp)
        for q in root.iter(tag('File')):
            if q.get('Type')=='TOUCH': q.set('Name',project.with_suffix('.snp').name)
        write_xml(tree,project)
        native_command([BIN/'soncmd.exe','-Read',project],folder/'read.log')
    record = dict(source=str(source),source_sha256=digest(source),project_sha256=digest(project),
                  wire_only_removed=True,calibration_probe_ports=[3,4],wire_probe_ports=[5,6],diagnostic_only=True)
    marker = out/'prepared.json'
    if marker.exists() and json.loads(marker.read_text()) != record:
        raise ValueError('Prepared calibration reference changed')
    marker.write_text(json.dumps(record,indent=2)+'\n')
    return project,out


def terminate_wire(s,f,l):
    """Terminate the last two equal-50-ohm ports with a floating series L."""
    z=2j*np.pi*np.asarray(f)*l*1e-9
    load=np.empty((len(f),2,2),complex)
    load[:,0,0]=load[:,1,1]=z/(z+100)
    load[:,0,1]=load[:,1,0]=100/(z+100)
    return s[:,:-2,:-2]+s[:,:-2,-2:]@load@np.linalg.solve(np.eye(2)[None]-s[:,-2:,-2:]@load,s[:,-2:,:-2])


def measured_l(s,f,l,*,local_gap=False):
    """Apply the same four-load resonance regression as the physical check."""
    if local_gap:
        # Open the long wire-terminal probe and place L across the original
        # short node/upper-ground calibration gap instead. This avoids a
        # 600 um component reference path in the ideal reference.
        opened=(s[:,:4,:4]+s[:,:4,4:]@np.linalg.solve(np.eye(2)[None]-s[:,4:,4:],s[:,4:,:4])
                if s.shape[1]==6 else s)
        reduced,z0=AP['reduce_pairs'](opened);y=AP['s_to_y'](reduced,z0)
        y[:,2,2]+=1/(2j*np.pi*f*l*1e-9)
    else:
        reduced,z0=AP['reduce_pairs'](terminate_wire(s,f,l));y=AP['s_to_y'](reduced,z0)
    held=np.isclose(f,5.5e9,rtol=0,atol=1);train=~held
    basis=lambda f:(np.asarray(f)[:,None]/5e9)**np.array([-1,1,3,5])
    co=np.linalg.solve(basis(f[train]),y[train].reshape(4,-1))
    def operator(f):return (basis(f)@co).reshape(len(f),3,3)
    validation=float(np.max(abs(operator(f[held])-y[held])))
    if validation>1e-5:raise ValueError('Held-out ideal calibration interpolation failed')
    desired=np.array([6.4,5.8,5.2,4.6])*1e9;yn=operator(desired)
    effective=yn[:,2,2]-(yn[:,2:3,:2]@np.linalg.solve(yn[:,:2,:2]+np.eye(2)[None]/50,yn[:,:2,2:3]))[:,0,0]
    loads=-effective.imag/(2*np.pi*desired)*1e12
    if np.any(loads<=0) or not np.all(np.diff(loads)>0):raise ValueError('Invalid ideal calibration loads')
    frequencies=[]
    for cap,center in zip(loads,desired):
        def objective(freq):
            return -abs(AP['add_capacitors'](operator(np.array([freq])),np.array([freq]),[cap*1000])[0,1,0])**2
        optimum=minimize_scalar(objective,bounds=(center-.08e9,center+.08e9),method='bounded',options={'xatol':1})
        if abs(optimum.x-center)>.075e9:raise ValueError('Calibration maximum is not bracketed')
        frequencies.append(float(optimum.x))
    fit=linregress(loads*1e-12,1/(2*np.pi*np.array(frequencies))**2)
    return dict(component_L_nH=l,measured_loop_L_nH=float(fit.slope*1e9),R_squared=float(fit.rvalue**2),
                held_out_Y_error_S=validation,load_intercept_C_fF=float(fit.intercept/fit.slope*1e15),
                placement='8 um calibration gap' if local_gap else '600 um former wire terminals')


def run(name):
    project,out=prepare(name);marker=out/'operator_completed.json'
    if marker.exists():
        saved=json.loads(marker.read_text())
        if saved != dict(project_sha256=digest(project),snp_sha256=digest(project.with_suffix('.snp'))):
            raise ValueError('Completed calibration reference changed')
        net=W['read_network'](project.with_suffix('.snp'),6)
    else:
        net=W['em'](project,6)
        marker.write_text(json.dumps(dict(project_sha256=digest(project),snp_sha256=digest(project.with_suffix('.snp'))),indent=2)+'\n')
    target=json.loads((ROOT/'output_data/taper_exploration/compact_centered'/('corrected_'+name)/'idc_Lcal/results.json').read_text())['L0_nH']
    same=measured_l(net.s,net.f,target)
    # Long component-port reference paths can give a substantial offset.
    # Choose a bracket only where positive loads and the held-out check pass;
    # a small trial L need not support four resonances in this frequency band.
    trials=[]
    for candidate in (.65,.75,.9,1.1,1.3,1.5,1.8):
        try:trials.append(measured_l(net.s,net.f,candidate))
        except ValueError:continue
    brackets=[(a['component_L_nH'],b['component_L_nH']) for a,b in zip(trials,trials[1:])
              if (a['measured_loop_L_nH']-target)*(b['measured_loop_L_nH']-target)<=0]
    if not brackets:raise ValueError(f'No validated ideal-L calibration bracket: {trials}')
    equal=brentq(lambda l:measured_l(net.s,net.f,l)['measured_loop_L_nH']-target,*brackets[0],xtol=1e-8)
    result=dict(diagnostic_only=True,physical_measured_L_nH=target,equal_component_value_control=same,
                calibrated_component_L_nH=equal,calibrated_check=measured_l(net.s,net.f,equal),validated_bracket_trials=trials,
                local_gap_check=measured_l(net.s,net.f,target,local_gap=True),
                limitation='Isolation reference check; does not establish assembled capacitances or accept a design.')
    (out/'reference_check.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result


def apply_to_full(name):
    """Historical L-only correction; excess long-span capacitance remains."""
    import skrf
    from geometry.taper.corrected_iteration import review
    from geometry.taper.ideal_wire_control import load_inductor
    base=OUT/name
    check_path=base/'calibration_reference_control/reference_check.json'
    calibration=json.loads(check_path.read_text())
    l=calibration['physical_measured_L_nH'];component_l=calibration['calibrated_component_L_nH']
    if abs(calibration['calibrated_check']['measured_loop_L_nH']-l)>1e-6:
        raise ValueError('Ideal component does not reproduce the physical isolation L')
    if calibration['calibrated_check']['R_squared']<.9999:
        raise ValueError('Poor calibrated ideal inductance regression')
    if not json.loads((base/'paired_ideal/open_control.json').read_text())['passed']:
        raise ValueError('Full ideal operator open control has not passed')
    project=ROOT/'sonnet_files/taper_exploration/compact_centered'/('corrected_'+name)/'matched_ideal_control/ideal_wire_operator.sonx'
    marker=json.loads((base/'paired_ideal/operator_completed.json').read_text())
    if marker!=dict(project_sha256=digest(project),snp_sha256=digest(project.with_suffix('.snp'))):
        raise ValueError('Completed full operator changed')
    n=W['read_network'](project.with_suffix('.snp'),4)
    out=base/'calibrated_ideal';out.mkdir(exist_ok=True)
    target=out/'effective_L_matched.s2p'
    skrf.Network(frequency=skrf.Frequency.from_f(n.f,unit='hz'),
                 s=load_inductor(n.s,n.f,component_l),z0=50).write_touchstone(str(target))
    result=review(target,out/'fit',l,data_label='EM + long-span ideal component (L calibrated only)')
    result.update(ideal_component_L_nH=component_l,physical_measured_L_nH=l,
                  ideal_component_reference=str(check_path),reference_sha256=digest(check_path),
                  effective_L_matched_in_isolation=True,additional_full_EM_fields=0,
                  long_span_capacitance_offset_not_corrected=True,
                  use_for_final_ideal_comparison=False)
    (out/'comparison.json').write_text(json.dumps(dict(physical=json.loads((base/'fit/review.json').read_text()),
        calibrated_ideal=result,nominal_component_control=json.loads((base/'paired_ideal/fit/review.json').read_text()),
        calibration=calibration,accepted=False),indent=2)+'\n')
    return result
