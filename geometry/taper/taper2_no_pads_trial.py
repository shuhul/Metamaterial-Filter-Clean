r"""Remove only T2's two ground pads from the current best case.

Run with: .venv\Scripts\python.exe -B -m geometry.taper.taper2_no_pads_trial
One four-thread solver at a time. Recheck L using explicit Sonnet capacitors,
then simulate full S21 even if removing the pads changes the extracted L.
"""
import os
os.environ['SONNET_MAX_JOBS'] = '1'
import json
from dataclasses import replace, asdict
import gdstk
from geometry.taper.taper2_idc import SETTINGS as baseline
from geometry.taper import taper2_plain as w

BASELINE = replace(baseline, coupling=replace(baseline.coupling, d_t1=50))
SETTINGS = replace(BASELINE, include_t2_pads=False)
CASE = 'idc_rail8_right30_dl256_w148_NO_T2_PADS_i42n70_j50_r74tg2'


def verify_only_pads_removed():
    """Boolean proof: no added metal; removed metal equals two original pads."""
    before, meta = w.build(BASELINE)
    after, _ = w.build(SETTINGS)
    assert not gdstk.boolean(after.get_polygons(), before.get_polygons(), 'not', precision=.001)
    removed = gdstk.boolean(before.get_polygons(), after.get_polygons(), 'not', precision=.001)
    body_ref = next(r for r in before.references if r.cell.name == 'ASYMMETRIC_T2')
    pad_refs = [r for r in body_ref.cell.references if r.cell.name == 'T2_PAD']
    assert len(pad_refs) == 2
    expected = [p.translate(*body_ref.origin) for r in pad_refs for p in r.get_polygons()]
    assert not gdstk.boolean(removed, expected, 'xor', precision=.001)
    expected_area = 2 * BASELINE.pad.w * BASELINE.pad.l
    assert abs(sum(p.area() for p in removed) - expected_area) < .001
    return dict(removed_polygons=2, removed_area_um2=expected_area,
                added_metal=False, all_other_metal_unchanged=True)


if __name__ == '__main__':
    proof = verify_only_pads_removed()
    project, out, _ = w.prepare(CASE, SETTINGS)
    (out/'proposal.json').write_text(json.dumps(dict(
        case=CASE, settings=asdict(SETTINGS), geometry_difference=proof,
        reason='User requested removing only the two T2 ground pads from the best j50 case.',
        max_concurrent_solvers=1), indent=2)+'\n')
    print(json.dumps(dict(project=str(project), geometry_proof=proof)), flush=True)
    measured = w.measure(CASE, SETTINGS)
    print(json.dumps(dict(stage='inductance', **measured), indent=2), flush=True)
    score, fitted = w.response(CASE, SETTINGS, measured['L0_nH'])
    print(json.dumps(dict(target=score, fit=fitted), indent=2), flush=True)
