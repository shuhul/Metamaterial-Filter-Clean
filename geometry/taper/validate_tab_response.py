"""Independent ordinary two-port validation of a tabbed T2 candidate.

The six-port screen is never used as the final simulation. This saves its
original source, changes only the ordinary source's sweep/output settings,
and compares a new two-port sweep with the fixed theoretical circuit.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
import json
import xml.etree.ElementTree as ET
import numpy as np
from geometry.taper.regular_cap_tab_case import ROOT, validate_native_tabbed
from geometry.taper.response_cases import W, A
from geometry.taper.tune_wide_response import impedance
from geometry.taper.sonnet import configure, simulate
from geometry.taper.compare_theory import compare
from geometry.prepare_sonnet import tag
from circuit_common import read_snp

NAME = 'response_match_span290_pairs31_rise32_right548_tab116'


def validate(name=NAME):
    out = ROOT / 'output_data/taper_exploration' / name
    project = ROOT / 'sonnet_files/taper_exploration' / name / 'generated_ncell_3.sonx'
    meta = json.loads((ROOT / 'geometry/taper/output' / (name+'.json')).read_text())
    validate_native_tabbed(project, meta)
    snapshot = out / 'screen_source_before_validation.sonx'
    manifest = json.loads((out / 'series_gap_probes/manifest.json').read_text())
    if not snapshot.exists():
        if W['digest'](project) != manifest['source_sha256']:
            raise ValueError('Screen source already changed without a saved original')
        snapshot.write_bytes(project.read_bytes())
    if W['digest'](snapshot) != manifest['source_sha256']:
        raise ValueError('Screen source snapshot is stale')
    original = ET.parse(snapshot)
    # Focus the ordinary adaptive solve on both narrow resonances and their
    # shoulders. The wider 4.5--5.5 GHz screening evidence is kept separately.
    configure(project, 4.8, 5.18, 3800)
    revised = ET.parse(project)
    for before, after in zip(list(original.getroot()), list(revised.getroot())):
        if before.tag not in (tag('Sweeps'), tag('OutputFiles')):
            if ET.tostring(before) != ET.tostring(after):
                raise ValueError('Ordinary validation modified a non-sweep/output setting')
    if len(list(revised.find(tag('Geometry')).iter(tag('Port')))) != 2:
        raise ValueError('Ordinary source must have exactly two geometry ports')
    runout = out / 'ordinary_validation'
    record = runout / 'run_status.json'
    done = json.loads(record.read_text()) if record.exists() else {}
    if done.get('status') == 'complete':
        if done['project_sha256'] != W['digest'](project) or done['snp_sha256'] != W['digest'](project.with_suffix('.snp')):
            raise ValueError('Completed ordinary validation was changed')
        snp = project.with_suffix('.snp')
    else:
        snp = simulate(project, runout)
    report = compare(snp, out)
    actual = read_snp(snp)
    predicted = A['short_internal_z'](impedance(name, actual.f))
    error = predicted - actual.s
    comparison = {
        'screen_source_sha256': W['digest'](snapshot),
        'ordinary_source_sha256': W['digest'](project),
        'only_sweep_and_output_settings_changed': True,
        'ordinary_geometry_ports': 2,
        'S21_screen_vs_ordinary_complex_RMSE': float(np.sqrt(np.mean(abs(error[:,1,0])**2))),
        'S21_screen_vs_ordinary_max_complex_error': float(abs(error[:,1,0]).max()),
        'all_S_screen_vs_ordinary_max_complex_error': float(abs(error).max()),
        'ordinary_snp_sha256': W['digest'](snp),
        'physical_acceptance': False,
        'note': 'Independent validation recorded. Acceptance requires reviewing the actual theoretical-response errors.'}
    (out / 'screen_vs_ordinary.json').write_text(json.dumps(comparison, indent=2)+'\n')
    print(json.dumps(comparison, indent=2), flush=True)
    return report


if __name__ == '__main__':
    validate()
