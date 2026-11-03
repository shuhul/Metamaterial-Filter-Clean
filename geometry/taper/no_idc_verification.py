"""Recheck the original 0.733 nH straight wire before any input IDC.

Reuse the ordinary physical response. New EM is limited to a short-gap ideal
reference and its independent open check. Never relax the regular-cell values.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
from copy import deepcopy
import json
import shutil
import xml.etree.ElementTree as ET
import gdstk
import numpy as np
import skrf
from scipy.optimize import brentq
from geometry.prepare_sonnet import tag, points, write_xml, BIN, run as native_command
from geometry.taper.corrected_iteration import ROOT, review
from geometry.taper.compact_stages import prepare as prepare_geometry
from geometry.taper.second_layout import from_dict
from geometry.taper.local_ideal_control import cached_em, ordinary_reference
from geometry.taper.calibration_reference_control import measured_l
from geometry.taper.ideal_wire_control import load_inductor
from geometry.taper.run_bridged_stage import digest
from geometry.taper.assess_corrected import best_target_delay
from geometry.taper.minimum_inductance import half_power_features as peaks, assess as assess_curve
from geometry.taper.progressive_fitting import transmission
from geometry.taper.corrected_iteration import REGULAR_C, REGULAR_G
from geometry.taper.staged_fitting import LT1, LR

NAME = 'w500_top10_wire52_gap6'
OUT = ROOT/'output_data/taper_exploration/no_idc_verification_2026-10-25'
NATIVE = ROOT/'sonnet_files/taper_exploration/no_idc_verification_2026-10-25'


def metal(path):
    return [gdstk.Polygon(points(q)) for q in ET.parse(path).iter(tag('PlanarPolygon'))]


def prepare():
    OUT.mkdir(parents=True, exist_ok=True)
    NATIVE.mkdir(parents=True, exist_ok=True)
    original = ROOT/'sonnet_files/taper_exploration/compact_centered'/NAME
    meta = json.loads((ROOT/'geometry/taper/output/compact_centered'/NAME/'direct/geometry.json').read_text())
    assert not meta['IDC_present'] and not meta['wire_removed']
    physical = original/'direct/generated_ncell_3.sonx'
    absent, _, _ = prepare_geometry('no_idc_check_2026_10_25_wire_absent',
                                    from_dict(meta['second']), 'direct', omit_wire=True)
    x = meta['cell_centers'][1]
    width = meta['second']['taper']['meander']['t']
    cut = [gdstk.rectangle((x-width/2, 466), (x+width/2, 1066))]
    if gdstk.boolean(gdstk.boolean(metal(physical), cut, 'not', precision=.001),
                     metal(absent), 'xor', precision=.001):
        raise ValueError('Absent reference changes more than the straight wire')
    isolated = NATIVE/'isolated_short_gap.sonx'
    if not isolated.exists():
        tree = ET.parse(original/'Lcal/generated_ncell_1.sonx')
        root = tree.getroot()
        wires = [q for q in root.iter(tag('PlanarPolygon')) if not
                 gdstk.boolean([gdstk.Polygon(points(q))], cut, 'xor', precision=.001)]
        if len(wires) != 1:
            raise ValueError('Expected a separate calibration wire')
        wire = wires[0]
        next(q for q in root.iter() if wire in list(q)).remove(wire)
        for q in root.iter(tag('File')):
            if q.get('Type') == 'TOUCH': q.set('Name', isolated.with_suffix('.snp').name)
        write_xml(tree, isolated)
        native_command([BIN/'soncmd.exe', '-Read', isolated], NATIVE/'isolated_read.log')
    full = NATIVE/'full_short_gap.sonx'
    if not full.exists():
        tree = ET.parse(absent); root = tree.getroot()
        comp = deepcopy(next(ET.parse(isolated).iter(tag('Component'))))
        comp.set('Id', str(max(int(q.get('Id')) for q in root.iter(tag('PlanarPolygon')))+100))
        comp.set('MacroID', 'noIDCShortGap'); comp.set('Name', 'NoIDCShortGap')
        level = next(q for q in root.find(tag('Geometry')).findall(tag('Level')) if q.get('Number') == '0')
        level.find(tag('Objects')).append(comp)
        sweeps = root.find(tag('Sweeps'))
        for q in list(sweeps): sweeps.remove(q)
        fs = ET.SubElement(ET.SubElement(sweeps, tag('Set'), On='TRUE'), tag('Frequencies'))
        ET.SubElement(fs, tag('Sweep'), Type='ABS_ENTRY', On='TRUE', Adaptive='TRUE',
                      Start='4.65', Stop='5.4', Target='7500')
        for q in root.iter(tag('File')):
            if q.get('Type') == 'TOUCH':
                q.set('Name', full.with_suffix('.snp').name); q.set('UseInterpolatedData', 'TRUE')
        write_xml(tree, full)
        native_command([BIN/'soncmd.exe', '-Read', full], NATIVE/'full_read.log')
    if gdstk.boolean(metal(full), metal(absent), 'xor', precision=.001):
        raise ValueError('Ideal control changed the metal')
    comp = next(ET.parse(full).iter(tag('Component')))
    if comp.get('Type') != 'PORTS' or [q.findtext(tag('Center')) for q in comp.findall(tag('Pin'))] != [f'({x},448)', f'({x},456)']:
        raise ValueError('Ideal terminals must occupy the original short 8 um gap')
    provenance = dict(physical_project=str(physical), physical_sha256=digest(physical),
                      isolated_sha256=digest(isolated), full_sha256=digest(full),
                      IDC_present=False, only_wire_removed=True, T1_downstream_unchanged=True)
    (OUT/'geometry_check.json').write_text(json.dumps(provenance, indent=2)+'\n')
    return physical, absent, isolated, full


def report():
    """Separate free-C curve agreement from the target-C response, using no EM."""
    import matplotlib.pyplot as plt
    cases = [('physical', 'Physical straight wire'),
             ('physical_meander_2p52', 'Physical meander'), ('ideal', 'Calibrated ideal reference')]
    cases = [(key, title) for key, title in cases if (OUT/key/'fit.json').exists()]
    fig, axes = plt.subplots(1, len(cases), figsize=(8*len(cases), 7), squeeze=False, sharey=True)
    results = {}
    for ax, (key, title) in zip(axes[0], cases):
        r = json.loads((OUT/key/'fit.json').read_text())
        if (r['free_capacitances'] != [['C', 1], ['G', 1], ['G', 2]]
            or r['C_fF'][0] != 0 or r['G_fF'][0] is not None
            or r['L_nH'][1:] != [LT1, LR, LR]
            or r['C_fF'][2:] != [REGULAR_C, REGULAR_C] or r['G_fF'][3] != REGULAR_G
            or digest(r['provenance']['source']) != r['provenance']['sha256']):
            raise ValueError('No-IDC fitting constraint or source provenance changed')
        data = np.load(OUT/key/'response.npz'); f = data['f']
        observed = data['actual'][:, 1, 0]; model = data['model'][:, 1, 0]
        base = transmission(f, r['L_nH'], [0, 100, REGULAR_C, REGULAR_C],
                            [None, 56.7, 8.2865, REGULAR_G], 0)
        delay, target = best_target_delay(f, base, observed)
        error = abs(target-observed)
        results[key] = dict(curve_pass=r['curve_pass'], complex_RMS=r['complex_RMS'],
            complex_max=r['complex_max'], L_T2_nH=r['L_nH'][0], C0_T1_fF=r['C_fF'][1],
            Cg_T2_T1_fF=r['G_fF'][1], Cg_T1_regular_fF=r['G_fF'][2],
            source_sha256=r['provenance']['sha256'], physical_peaks=peaks(f, observed),
            fitted_peaks=peaks(f, model), target_peaks=peaks(f, target),
            target_only_delay_ps=delay, target_complex_RMS=float(np.sqrt(np.mean(error**2))),
            target_complex_max=float(error.max()),
            scope='No IDC; T2 C0 fixed zero; all L and regular-cell parameters fixed. Target curve uses measured L, not 0.74 nH for the 2.52 nH experiment.')
        results[key]['fixed_parameter_and_source_checks_passed'] = True
        results[key]['peak_curve_assessment'] = assess_curve(dict(
            peaks=[], model_peaks=[], S21_complex_RMSE=r['complex_RMS'],
            S21_max_complex_error=r['complex_max'], active_bounds=r['active_bounds']), data)
        ax.plot(f/1e9, abs(observed), color='#b22222', lw=2.7, label='EM response' if key != 'ideal' else 'EM + calibrated ideal L')
        ax.plot(f/1e9, abs(model), 'k--', lw=2.1, label='Free T1 capacitance fit')
        ax.plot(f/1e9, abs(target), ':', color='#1565c0', lw=2.1, label='Target capacitances; measured L')
        box = (f"L = {r['L_nH'][0]:.5f} nH (fixed)\n"
               f"T1 C0 = {r['C_fF'][1]:.2f} fF\nT2–T1 Cg = {r['G_fF'][1]:.2f} fF\n"
               f"T1–unit Cg = {r['G_fF'][2]:.3f} fF\nComplex RMS = {r['complex_RMS']:.5f}")
        ax.text(.025, .97, box, transform=ax.transAxes, va='top', fontsize=13,
                bbox=dict(boxstyle='round,pad=.5', fc='white', ec='gray', alpha=.95))
        ax.set(title=title, xlabel='Frequency (GHz)', ylim=(0, 1.25))
        ax.title.set_size(20); ax.xaxis.label.set_size(17); ax.tick_params(labelsize=14)
        ax.legend(loc='upper right', fontsize=10); ax.grid(alpha=.22)
    axes[0, 0].set_ylabel(r'$|S_{21}|$', fontsize=23)
    fig.suptitle('No-IDC Verification: Fixed Inductances and Regular Cells', fontsize=22)
    fig.tight_layout(); fig.savefig(OUT/'no_idc_comparison.png', dpi=170); plt.close(fig)
    (OUT/'assessment.json').write_text(json.dumps(results, indent=2)+'\n')
    return results


def reuse_ordinary(absent):
    """Reuse a completed identical-physics ordinary sweep, with explicit provenance."""
    old = ROOT/'sonnet_files/taper_exploration/compact_centered/w500_wire52_gap6_wire_absent/direct/generated_ncell_3.sonx'
    record = json.loads((ROOT/'output_data/taper_exploration/compact_centered/w500_wire52_gap6_wire_absent/direct/run_status.json').read_text())
    if record['status'] != 'complete' or record['project_sha256'] != digest(old) or record['snp_sha256'] != digest(old.with_suffix('.snp')):
        raise ValueError('Historical ordinary response provenance failed')
    a, b = ET.parse(old).getroot(), ET.parse(absent).getroot()
    ignored = {tag(q) for q in ('Header', 'Sweeps', 'OutputFiles')}
    if a.attrib != b.attrib or [ET.tostring(q) for q in a if q.tag not in ignored] != [ET.tostring(q) for q in b if q.tag not in ignored]:
        raise ValueError('Ordinary reference differs in geometry, materials, ports, mesh or solver controls')
    # The existing reader recognizes native .snp files using a temporary .s2p.
    from geometry.taper.response_cases import W
    net = W['read_network'](old.with_suffix('.snp'), 2)
    if len(net.f) != 7501 or not np.allclose(net.f, np.linspace(4.65e9, 5.4e9, 7501), rtol=0, atol=1):
        raise ValueError('Ordinary reference lacks the complete matching frequency grid')
    shutil.copy2(old.with_suffix('.snp'), absent.with_suffix('.snp'))
    reuse = dict(source_project=str(old), source_project_sha256=digest(old),
                 source_snp_sha256=digest(old.with_suffix('.snp')), matching_physics=True,
                 excluded_sections=['Header', 'Sweeps', 'OutputFiles'], samples=7501,
                 new_ordinary_EM_frequencies=0)
    (OUT/'ordinary_reuse.json').write_text(json.dumps(reuse, indent=2)+'\n')
    (OUT/'ordinary_completed.json').write_text(json.dumps(dict(
        project_sha256=digest(absent), snp_sha256=digest(absent.with_suffix('.snp'))), indent=2)+'\n')


def run():
    physical, absent, isolated, full = prepare()
    reuse_ordinary(absent)
    cal_path = ROOT/'output_data/taper_exploration/compact_centered'/NAME/'Lcal/results.json'
    l = json.loads(cal_path.read_text())['L0_nH']
    source = ROOT/'output_data/taper_exploration/compact_centered'/NAME/'direct/response.snp'
    physical_fit = review(source, OUT/'physical', l, idc=False)
    meander = ROOT/'sonnet_files/taper_exploration/minimum_inductance_2026-10-23/meander_w80/direct/generated_ncell_3.snp'
    review(meander, OUT/'physical_meander_2p52', 2.5200619105580704, idc=False)
    report()
    cal = cached_em(isolated, OUT, 'isolation')
    component = brentq(lambda v: measured_l(cal.s, cal.f, v, local_gap=True)['measured_loop_L_nH']-l,
                       .5, 1.1, xtol=1e-8)
    check = measured_l(cal.s, cal.f, component, local_gap=True)
    if check['R_squared'] < .9999:
        raise ValueError('Ideal isolation calibration failed')
    (OUT/'calibration.json').write_text(json.dumps(check, indent=2)+'\n')
    print('NO IDC IDEAL CALIBRATION', check, flush=True)
    net = cached_em(full, OUT, 'full')
    opened = net.s[:, :2, :2]+net.s[:, :2, 2:]@np.linalg.solve(
        np.eye(2)[None]-net.s[:, 2:, 2:], net.s[:, 2:, :2])
    ref, reused = ordinary_reference('no_idc_check_2026_10_25', absent, OUT, net, opened)
    ix = np.array([np.argmin(abs(net.f-f)) for f in ref.f])
    if not np.allclose(net.f[ix], ref.f, rtol=0, atol=1):
        raise ValueError('Open-check frequencies differ')
    err = abs(opened[ix]-ref.s)
    validation = dict(S21_RMS=float(np.sqrt(np.mean(err[:, 1, 0]**2))),
                      S21_max=float(err[:, 1, 0].max()), all_S_max=float(err.max()),
                      frequencies_GHz=(ref.f/1e9).tolist(), ordinary_fields_reused=reused)
    validation['passed'] = validation['S21_RMS'] < .005 and validation['S21_max'] < .05
    (OUT/'open_control.json').write_text(json.dumps(validation, indent=2)+'\n')
    if not validation['passed']:
        raise ValueError(f'No-IDC ideal control failed validation: {validation}')
    path = OUT/'calibrated_ideal.s2p'
    skrf.Network(frequency=skrf.Frequency.from_f(net.f, unit='hz'),
                 s=load_inductor(net.s, net.f, component), z0=50).write_touchstone(str(path))
    ideal_fit = review(path, OUT/'ideal', l, idc=False,
                       data_label='EM + calibrated short-gap ideal L')
    result = dict(physical=physical_fit, ideal=ideal_fit, independent_L_nH=l,
                  ideal_component_L_nH=component, open_check=validation,
                  final_T2_accepted=False, IDC_present=False)
    (OUT/'comparison.json').write_text(json.dumps(result, indent=2)+'\n')
    report()
    print('NO IDC COMPLETE', dict(physical_RMS=physical_fit['complex_RMS'],
          ideal_RMS=ideal_fit['complex_RMS'], open_RMS=validation['S21_RMS'],
          open_max=validation['S21_max'], output=str(OUT)), flush=True)
    return result


if __name__ == '__main__':
    import sys
    prepare() if '--prepare' in sys.argv else run()
