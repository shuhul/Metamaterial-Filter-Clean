"""Incremental T3 design around the selected, immutable T2 assembly.

All dimensions are um. Native EM uses the existing guarded serial four-thread
runner. Calibration uses explicit Sonnet capacitors, never a sparse Y model.
"""
import os
for key in ('OPENBLAS_NUM_THREADS', 'OMP_NUM_THREADS'):
    os.environ[key] = '1'
os.environ['SONNET_MAX_JOBS'] = '1'
from dataclasses import dataclass, field, asdict, replace
from pathlib import Path
from types import SimpleNamespace
from hashlib import sha256
import json
import shutil
import xml.etree.ElementTree as ET
import numpy as np
import gdstk
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from geometry.taper import taper2_plain as previous
from geometry.taper.taper2_asymmetric_preview import Cap, Clearance
from geometry.taper.taper2_v3 import InputCoupler
from geometry.meander import MeanderParams
from geometry.pad import PadParams, make_pad
from geometry.placement import place
from geometry.prepare_sonnet import main as import_native, tag, points
from geometry.taper.audit import check_geometry
from geometry.taper.sonnet import configure, simulate
from geometry.taper.explicit_inductance import measure as explicit_measure
from geometry.taper.progressive_fitting import fit as fit_circuit
from geometry.specifications import TAPER_3, TAPER_2, TAPER_1, UNIT_CELL

ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = 'three_tapers_v3_T3_2026-10-30'
NATIVE = ROOT/'sonnet_files/taper_exploration'/CAMPAIGN
OUT = ROOT/'output_data/taper_exploration'/CAMPAIGN
REFERENCE_FILE = ROOT/'geometry/taper/selected_taper2_cap5.json'
REFERENCE = json.loads(REFERENCE_FILE.read_text())
SOURCE = Path(REFERENCE['project'])
SOURCE_GDS = SOURCE.parent/'geometry.gds'
SOURCE_META = json.loads((SOURCE.parent/'geometry.json').read_text())
REFERENCE_FIT = json.loads((SOURCE.parent.parent/'evidence/fit/fit.json').read_text())
TOLERANCE_PERCENT = 5.0


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def body_default():
    # The two cap rails remain thin. Independent clearances protect the input
    # IDC from the meander. Existing rectangular ground pads tune C0.
    return previous.Settings(
        cap=Cap(l=590, t_top=14, t_left=8, t_right=8),
        clearance=Clearance(d_left=256, d_right=124),
        meander=MeanderParams(w=84), pad=PadParams(w=24, l=586),
        d_pad=6, d_pad_left=148,
        coupling=InputCoupler(d=56, n=0, t_input=8, d_t1=66,
                              d_regular=74, t_regular_ground=2))


@dataclass(frozen=True)
class Params:
    body: previous.Settings = field(default_factory=body_default)
    d_input: int = 36
    n_input: int = 70
    t_input_rail: int = 8
    d_t2: int = 22
    n_t2: int = 70
    t_finger: int = 4
    d_finger: int = 4
    d_tip: int = 2
    d_top: int = 8  # T3 upper metal to the unchanged flat top ground
    w_side_left: int = 412  # Full-layout left ground rectangle; weak L fixture unchanged
    pad_right: PadParams | None = None  # None shares body.pad; optional independent right ground pad


def settings_dict(p):
    """Keep legacy default metadata byte-for-byte compatible in content."""
    values = asdict(p)
    if p.pad_right is None:
        values.pop("pad_right")
    if p.d_top == 8:
        values.pop("d_top")
    if p.w_side_left == 412:
        values.pop("w_side_left")
    return values


def from_dict(d):
    q = dict(d)
    q['body'] = previous.from_dict(q['body'])
    if q.get('pad_right') is not None:
        q['pad_right'] = PadParams(**q['pad_right'])
    return Params(**q)


def require_reference():
    if not REFERENCE['T2_capacitance_goal_met']:
        raise ValueError('Accepted T2 is required')
    if digest(REFERENCE_FILE) != REFERENCE_HASH:
        raise ValueError('T2 reference changed during this run')
    if digest(SOURCE) != REFERENCE['result']['project_sha256'] or digest(SOURCE.with_suffix('.snp')) != REFERENCE['result']['snp_sha256']:
        raise ValueError('Selected T2 native files changed')
    proof = json.loads((SOURCE.parent/'validation.json').read_text())
    if digest(SOURCE_GDS) != proof['source_gds_sha256']:
        raise ValueError('Selected T2 source GDS changed')


REFERENCE_HASH = digest(REFERENCE_FILE)


def isolated_body(s, *, d_top=8, pad_right=None):
    # The old scaffold assumes an8um top gap. Adjust its cap dimensions so
    # its underside/bottom/coil/pads are already at the desired coordinates,
    # then change just the central body's upper metal edge. The weak fixture
    # feeds and every production downstream cell retain their old positions.
    if not isinstance(d_top, int) or d_top < 2 or d_top % 2:
        raise ValueError('T3 top ground clearance must be positive even um')
    lift = 8-d_top
    scaffold = replace(s, cap=replace(s.cap, l=s.cap.l-lift, t_top=s.cap.t_top-lift))
    if scaffold.cap.t_top < 2 or scaffold.cap.l <= scaffold.cap.t_top:
        raise ValueError('Invalid cap dimensions for the requested top clearance')
    if pad_right is not None:
        # A small temporary pad lets the established scaffold validate both
        # independent gaps. Replace only its two pad references below.
        scaffold=replace(scaffold,pad=PadParams(w=min(s.pad.w,pad_right.w),l=min(s.pad.l,pad_right.l)))
    top, meta = previous.build(scaffold, isolation=True)
    body = next(r.cell for r in top.references if r.cell.name == 'ASYMMETRIC_T2')
    # These cells coexist with T2 cells in production. Unique names prevent a
    # new coil silently replacing a same-named dependency during GDS export.
    if d_top != 8:
        if len(body.polygons) != 1:
            raise ValueError('Expected one U-shaped cap polygon')
        old = body.polygons[0]
        xy = old.points.copy()
        mask = np.isclose(xy[:,1], -8)
        if int(mask.sum()) != 2:
            raise ValueError('Unexpected cap upper edge')
        xy[mask,1] = -d_top
        body.remove(old)
        body.add(gdstk.Polygon(xy, layer=old.layer, datatype=old.datatype))
        meta['unit']['cap'].update(asdict(s.cap))
        meta['unit']['d_top'] = d_top
        meta['parameters']['unit']['d_top'] = d_top
        meta['parameters']['unit']['cap']['t_top'] = s.cap.t_top
        meta['parameters']['d_feed_top'] = 8
    if pad_right is not None:
        width=meta['unit']['cap']['w']
        left_gap=s.d_pad if s.d_pad_left is None else s.d_pad_left
        for pad,gap,clearance in ((s.pad,left_gap,s.clearance.d_left),(pad_right,s.d_pad,s.clearance.d_right)):
            if gap<2 or gap%2 or clearance-gap-pad.w<2 or pad.l>s.meander.l-2:
                raise ValueError('Independent pad must clear its cap rail, coil and cap underside')
        refs=[r for r in body.references if r.cell.name=='T2_PAD']
        if len(refs)!=2:raise ValueError('Expected exactly two scaffold pad references')
        for ref in refs:body.remove(ref)
        bottom=-d_top-s.cap.t_top-s.meander.l
        body.add(place(make_pad(s.pad,name='PAD_LEFT'),origin=(-width/2+s.cap.t_left+left_gap,bottom),anchor='bottom_left'),
                 place(make_pad(pad_right,name='PAD_RIGHT'),origin=(width/2-s.cap.t_right-s.d_pad,bottom),anchor='bottom_right'))
        for unit in (meta['unit'],meta['parameters']['unit']):
            unit['pad']=asdict(s.pad)
            unit['pad_right']=asdict(pad_right)
        if any(not np.allclose(q.points/2,np.round(q.points/2)) for q in body.get_polygons()):
            raise ValueError('Independent T3 pads are off the2um grid')
    body.name = 'T3_BODY' 
    for c in body.dependencies(True):
        c.name = 'T3_'+c.name
    return top, body, meta


def build(p=Params(), *, stage='response'):
    require_reference()
    if not isinstance(p.w_side_left, int) or p.w_side_left < 412 or p.w_side_left % 2:
        raise ValueError('Left ground width must be even and at least412um for the400um reference planes')
    iso, body, raw = isolated_body(p.body, d_top=p.d_top, pad_right=p.pad_right)
    if stage == 'Lcal':
        raw.update(target_stage='taper_3', layout_kind='T3_weak_fixture_explicit_L',
                   t3_settings=settings_dict(p), T2_reference_sha256=REFERENCE_HASH, stage=stage)
        return iso, raw
    if stage not in ('direct', 'response'):
        raise ValueError('Unknown T3 stage')
    for v in [p.d_input, p.d_t2, p.t_finger, p.d_finger, p.d_tip, p.t_input_rail]:
        if v < 2 or v % 2:
            raise ValueError('Coupler dimensions must be positive even um')
    for n in (p.n_input, p.n_t2):
        if not isinstance(n, int) or n < 0 or n == 1:
            raise ValueError('Use zero fingers for plain rails, or at least two for an IDC')
    if p.body.cap.t_right != 8:
        raise ValueError('T3 right rail must match the selected T2 feed rail width8')
    width = raw['unit']['cap']['w']
    # Honor the existing coupler feed-span control in the full T3 assembly,
    # as the underlying weak fixture already does. None preserves legacy metal.
    input_span = p.body.coupling.w_input_stub
    if input_span is None:
        input_span = width
    if input_span < 2 or input_span % 2:
        raise ValueError('Input feed span must be positive even um')
    input_edge = 426+input_span
    if p.w_side_left > input_edge-p.t_input_rail-2:
        raise ValueError("Left ground would touch the input rail")
    left = input_edge+p.d_input
    right = left+width
    old_left = SOURCE_META['cell_lefts'][1]
    t2_left = right+p.d_t2
    shift = t2_left-old_left
    if shift % 2:
        raise ValueError('Frozen downstream group must translate by even um')
    box_width = SOURCE_META['w_sim']+shift
    source_top = gdstk.read_gds(str(SOURCE_GDS)).top_level()[0]
    top = gdstk.Cell('TOP')
    for ref in source_top.references:
        if ref.bounding_box()[0][0] < old_left-.001:
            raise ValueError('Unexpected upstream component in frozen T2 GDS')
        copied = ref.copy()
        copied.origin = (ref.origin[0]+shift, ref.origin[1])
        top.add(copied)
    # Retain exactly the selected downstream polygons. Replace only the old
    # feed electrode and its owned input coupling, now the T3-to-T2 junction.
    cut = gdstk.rectangle((old_left, -1052), (SOURCE_META['w_sim'], 448))
    downstream = gdstk.boolean(source_top.polygons, [cut], 'and', precision=.001, layer=1)
    for poly in downstream:
        poly.translate(shift, 0)
    parts = downstream+[
        gdstk.rectangle((0, 0), (box_width, 448)),
        gdstk.rectangle((0, -1052), (box_width, -618)),
        gdstk.rectangle((0, -1052), (p.w_side_left, -24)),
    ]
    top.add(gdstk.Reference(body, origin=(left+width/2, 0)))

    def idc(edge, gap, n):
        if n == 0:
            return
        if gap <= 2*p.d_tip:
            raise ValueError('IDC gap must exceed twice its tip clearance')
        upper = -p.d_top-p.body.cap.t_top-2
        lower = upper-(n-1)*(p.t_finger+p.d_finger)-p.t_finger
        if lower < -p.d_top-p.body.cap.l+2:
            raise ValueError('IDC extends below its rails')
        for i in range(n):
            y = upper-i*(p.t_finger+p.d_finger)
            x0, x1 = ((edge+p.d_tip, edge+gap) if i % 2 else (edge, edge+gap-p.d_tip))
            parts.append(gdstk.rectangle((x0, y-p.t_finger), (x1, y)))
    if stage == 'response':
        parts += [gdstk.rectangle((0, -18), (input_edge, -8)),
                  gdstk.rectangle((input_edge-p.t_input_rail, -p.d_top-p.body.cap.l), (input_edge, -p.d_top))]
        idc(input_edge, p.d_input, p.n_input)
    else:
        # No IDC or redundant upstream cap scaffold in the direct-feed test.
        parts.append(gdstk.rectangle((0, -18), (left+width/2, -8)))
    idc(right, p.d_t2, p.n_t2)
    top.add(*gdstk.boolean(parts, [], 'or', precision=.001, layer=1))
    meta = dict(layout_kind='T3_'+stage+'_frozen_T2', stage=stage, isolation=False,
        unit=raw['unit'], t3_settings=settings_dict(p), T2_reference_sha256=REFERENCE_HASH,
        downstream_shift_um=shift, downstream_left_um=t2_left, w_sim=box_width,
        cell_centers=[input_edge-input_span/2, left+width/2]+[v+shift for v in SOURCE_META['cell_centers'][1:]],
        inductance_pin_x=left+width/2+raw['unit']['meander']['x_offset_from_cap_center'],
        expected_metal_components=2 if stage=='response' else 1,
        parameters=dict(n_cells=4, l_sim=1500, l_top_ground=448, w_side=412, t_feed=10,
                        unit=dict(d_top=p.d_top, cap=dict(t_top=p.body.cap.t_top))))
    if p.d_top != 8:
        meta["parameters"]["d_feed_top"] = 8
    if p.w_side_left != 412:
        meta["parameters"]["w_side_left"] = p.w_side_left
    require_downstream(top.get_polygons(), meta)
    return top, meta


def require_downstream(polys, meta, *, native=False):
    shift = meta['downstream_shift_um']
    boundary = meta['downstream_left_um']
    y0, y1 = (0, 1500) if native else (-1052, 448)
    roi = gdstk.rectangle((boundary, y0), (meta['w_sim'], y1))
    original = ([gdstk.Polygon(points(q)) for q in ET.parse(SOURCE).iter(tag('PlanarPolygon'))]
                if native else gdstk.read_gds(str(SOURCE_GDS)).top_level()[0].get_polygons())
    expected = [p.copy().translate(shift, 0) for p in original]
    a = gdstk.boolean(polys, [roi], 'and', precision=.001)
    b = gdstk.boolean(expected, [roi], 'and', precision=.001)
    if gdstk.boolean(a, b, 'xor', precision=.001):
        raise ValueError('Selected T2 and downstream metal changed')


def validate(path, meta):
    top, expected = build(from_dict(meta['t3_settings']), stage=meta['stage'])
    if meta != expected:
        raise ValueError('T3 geometry metadata differs')
    lib = gdstk.read_gds(str(path)); tops = lib.top_level()
    if len(tops) != 1 or tops[0].name != 'TOP' or not np.isclose(lib.unit, 1e-6):
        raise ValueError('Wrong top cell or GDS units')
    polys = tops[0].get_polygons()
    if any(q.layer != 1 or q.datatype != 0 for q in polys):
        raise ValueError('Wrong metal layer')
    if gdstk.boolean(polys, top.get_polygons(), 'xor', precision=.001):
        raise ValueError('GDS differs from parameters')
    if len(gdstk.boolean(polys, [], 'or', precision=.001)) != meta['expected_metal_components']:
        raise ValueError('Unexpected T3 metal connectivity')
    if meta['stage'] == 'Lcal' and any(not np.allclose(q.points/2, np.round(q.points/2)) for q in polys):
        raise ValueError('T3 layout is not on the2um grid')
    if meta['stage'] != 'Lcal':
        require_downstream(polys, meta)
    return SimpleNamespace(n_cells=meta['parameters']['n_cells'], l_sim=1500, l_top_ground=448), polys


def prepare(case, p=Params(), *, stage='response'):
    native = NATIVE/case/stage; out = OUT/case/stage
    native.mkdir(parents=True, exist_ok=True); out.mkdir(parents=True, exist_ok=True)
    top, meta = build(p, stage=stage)
    # Old isolation builder calls its top TOP already; production does too.
    gds = native/'geometry.gds'; saved = gds.with_suffix('.json')
    if saved.exists() and json.loads(saved.read_text()) != meta:
        raise ValueError('Use a fresh case name for changed geometry')
    if not gds.exists():
        lib = gdstk.Library(unit=1e-6, precision=1e-9)
        lib.add(top, *top.dependencies(True)); lib.write_gds(str(gds))
        saved.write_text(json.dumps(meta, indent=2)+'\n')
    validate(gds, meta)
    n = meta['parameters']['n_cells']
    project = native/f'generated_ncell_{n}.sonx'
    if not project.exists():
        import_native(n, gds_path=gds, output_dir=native, geometry_validator=validate)
    check_geometry(project, gds, meta, geometry_validator=validate)
    if stage != 'Lcal':
        require_downstream([gdstk.Polygon(points(q)) for q in ET.parse(project).iter(tag('PlanarPolygon'))], meta, native=True)
    fig, ax = plt.subplots(figsize=(19, 6))
    for poly in top.get_polygons():
        ax.add_patch(Polygon(poly.points, fc='#297b54', ec='none'))
    ax.set(xlim=(360, meta['w_sim']-360), ylim=(-650, 45), aspect='equal',
           xlabel='x (um)', ylabel='y (um)', title=f'{case}: {stage} — T3 then selected T2/T1/regular')
    fig.tight_layout(); fig.savefig(out/'geometry.png', dpi=170); plt.close(fig)
    (out/'settings.json').write_text(json.dumps(settings_dict(p), indent=2)+'\n')
    return project, out, meta


def measure(case, p=Params(), *, caps_fF=(500,700,900,1100)):
    project, out, meta = prepare(case, p, stage='Lcal')
    for saved in sorted(OUT.glob('*/inductance.json')):
        if saved.parent.name == case:
            continue
        old = json.loads(saved.read_text())
        source_case = old.get('reused_from_case', saved.parent.name)
        if source_case == case:
            continue  # A sibling reuse record must not create self provenance.
        source_native = NATIVE/source_case/'Lcal'
        source_meta = json.loads((source_native/'geometry.json').read_text())
        if source_meta['unit'] != meta['unit']:
            continue
        source_gds = gdstk.read_gds(str(source_native/'geometry.gds')).top_level()[0]
        own_gds = gdstk.read_gds(str(project.parent/'geometry.gds')).top_level()[0]
        if gdstk.boolean(source_gds.get_polygons(), own_gds.get_polygons(), 'xor', precision=.001):
            continue
        proof_path = OUT/source_case/'Lcal_explicit/completed.json'
        proof = json.loads(proof_path.read_text())
        for path, hashed in proof['files_sha256'].items():
            if digest(path) != hashed:
                raise ValueError('Source explicit calibration changed')
        result = json.loads((OUT/source_case/'Lcal_explicit/results.json').read_text())
        if result['Python_frequency_approximation'] or abs(result['target_L0_nH']-TAPER_3.L0_nH)>1e-12:
            raise ValueError('Wrong source calibration method/target')
        result.update(case=case, reused_from_case=source_case,
                      identical_body_and_weak_fixture=True,
                      source_completed_sha256=digest(proof_path))
        (OUT/case/'inductance.json').write_text(json.dumps(result, indent=2)+'\n')
        print('REUSE EXPLICIT L', source_case, result['L0_nH'], flush=True)
        return result
    # Geometry-specific case identity prevents accidentally mixing bodies.
    result = explicit_measure(project, project.parent.parent/'Lcal_explicit', out.parent/'Lcal_explicit',
        name=case, target_L_nH=TAPER_3.L0_nH, caps_fF=list(caps_fF), metadata=meta)
    (OUT/case/'inductance.json').write_text(json.dumps(result, indent=2)+'\n')
    return result


def fit_response(case, p, l0, *, stage='response'):
    project, out, meta = prepare(case, p, stage=stage)
    status = out/'run_status.json'
    if status.exists() and json.loads(status.read_text()).get('status') == 'complete':
        old = json.loads(status.read_text())
        if digest(project) != old['project_sha256'] or digest(project.with_suffix('.snp')) != old['snp_sha256']:
            raise ValueError('Completed T3 response changed')
    else:
        configure(project, 4.65, 5.4, 7500); simulate(project, out)
    ls = [l0]+REFERENCE_FIT['L_nH']
    cs = [TAPER_3.C0_fF]+REFERENCE_FIT['C_fF']
    gs = [None if stage=='direct' else TAPER_3.Cg_input_fF]+REFERENCE_FIT['G_fF']
    free = [('C',0),('G',1)] if stage=='direct' else [('C',0),('G',0),('G',1)]
    seeds = ([[100,153],[500,153],[900,220],[251,500]] if stage=='direct' else
             [[100,261,153],[176,261,153],[400,261,153],[600,400,220],[900,261,153]])
    fitted = fit_circuit(project.with_suffix('.snp'), out/'fit', ls, cs, gs, free,
        capacitance_seeds=seeds,
        title='Third Taper Cell: '+('Direct Feed' if stage=='direct' else 'Transmission'))
    values = [fitted['C_fF'][0], fitted['G_fF'][1]] if stage=='direct' else [fitted['C_fF'][0],fitted['G_fF'][0],fitted['G_fF'][1]]
    targets = [TAPER_3.C0_fF,TAPER_2.Cg_input_fF] if stage=='direct' else [TAPER_3.C0_fF,TAPER_3.Cg_input_fF,TAPER_2.Cg_input_fF]
    errors = (100*(np.array(values)/targets-1)).tolist()
    result = dict(case=case, stage=stage, values_fF=values, targets_fF=targets, errors_percent=errors,
        max_cap_error_percent=max(abs(v) for v in errors), fit_RMS=fitted['complex_RMS'],
        components_pass=max(abs(v) for v in errors)<=TOLERANCE_PERCENT,
        fit_pass=fitted['curve_pass'] and fitted['interior_free_parameters'],
        measured_L_nH=l0, T2_reference_sha256=REFERENCE_HASH, project=str(project),
        all_older_parameters_fixed=True, objective='New Cs within5%; fit-to-EM RMS<=.025; L within2%. Nominal-response error diagnostic only.')
    (out/'result.json').write_text(json.dumps(result, indent=2)+'\n')
    rows = [json.loads(q.read_text()) for q in OUT.glob('*/'+stage+'/result.json')]
    rows.sort(key=lambda q:q['max_cap_error_percent'])
    (OUT/(stage+'_ledger.json')).write_text(json.dumps(rows, indent=2)+'\n')
    print('T3 RESULT', json.dumps(result), flush=True)
    return result


def fit_older_within_one_percent(case, l0, *, stage='direct'):
    """User-authorized alternative fit; older Cs stay within +/-1% of T2."""
    project = NATIVE/case/stage/'generated_ncell_4.sonx'
    out = OUT/case/stage/'fit_older_1percent'
    ls = [l0]+REFERENCE_FIT['L_nH']
    cs = [TAPER_3.C0_fF]+REFERENCE_FIT['C_fF']
    gs = [None if stage=='direct' else TAPER_3.Cg_input_fF]+REFERENCE_FIT['G_fF']
    new = [('C',0),('G',1)] if stage=='direct' else [('C',0),('G',0),('G',1)]
    older = [('C',i) for i in range(1,len(cs))]+[('G',i) for i in range(2,len(gs))]
    free = new+older
    bounds = [(0.,1000.) if k=='C' else (.1,3000.) for k,i in new]
    for k,i in older:
        value = (cs if k=='C' else gs)[i]
        bounds.append((.99*value,1.01*value))
    fitted = fit_circuit(project.with_suffix('.snp'), out, ls, cs, gs, free,
        capacitance_bounds=bounds,
        title='Third Taper: Older Capacitances Bounded to +/-1%')
    values = [fitted[k+'_fF'][i] for k,i in new]
    targets = [TAPER_3.C0_fF,TAPER_2.Cg_input_fF] if stage=='direct' else [TAPER_3.C0_fF,TAPER_3.Cg_input_fF,TAPER_2.Cg_input_fF]
    errors = (100*(np.array(values)/targets-1)).tolist()
    changes = [dict(kind=k,index=i,reference_fF=(cs if k=='C' else gs)[i],
        fitted_fF=fitted[k+'_fF'][i],
        difference_percent=100*(fitted[k+'_fF'][i]/(cs if k=='C' else gs)[i]-1))
        for k,i in older]
    older_pass = all(abs(q['difference_percent'])<=1.0000001 for q in changes)
    # An older C at its allowed limit is valid; a new C at a search limit is not.
    new_interior = not any(fitted['active_bounds'][:len(new)]+fitted['active_bounds'][-1:])
    summary = dict(case=case,stage=stage,fit_folder='fit_older_1percent',values_fF=values,
        targets_fF=targets,errors_percent=errors,max_cap_error_percent=max(abs(v) for v in errors),
        components_pass=max(abs(v) for v in errors)<=TOLERANCE_PERCENT,
        fit_RMS=fitted['complex_RMS'],fit_pass=fitted['curve_pass'] and new_interior and older_pass,
        measured_L_nH=l0,all_older_parameters_fixed=False,all_inductances_fixed=True,
        older_capacitance_changes=changes,older_capacitances_within_one_percent=older_pass,
        T2_reference_sha256=REFERENCE_HASH,project=str(project),
        policy='Separate user-authorized fit; +/-1% relative to saved T2, never accumulated across trials.')
    (out/'result.json').write_text(json.dumps(summary,indent=2)+'\n')
    return fitted


if __name__ == '__main__':
    project, out, meta = prepare('initial_w84_pad24x586_gap6_i36_j22', stage='Lcal')
    print(project)
