"""Tune asymmetric T2, with an optional explicitly authorized input-only IDC.

Edit SETTINGS and CASE here, then run as a module. All geometry is immutable
per case. Calibration and full response share exactly the same T2 body.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from dataclasses import dataclass, field, asdict, replace
from types import SimpleNamespace
from pathlib import Path
import json
import shutil
import xml.etree.ElementTree as ET
import gdstk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

from geometry.taper import taper2_v3 as common
from geometry.taper.taper2_asymmetric_preview import Cap, Clearance
from geometry.meander import MeanderParams, make_meander
from geometry.pad import PadParams, make_pad
from geometry.placement import place
from geometry.prepare_sonnet import main as import_native, tag, points
from geometry.taper.audit import check_geometry

ROOT = common.ROOT
CAMPAIGN = 'three_tapers_v3_T2_plain_2026-10-28'


@dataclass(frozen=True)
class Settings:
    cap: Cap = field(default_factory=Cap)
    clearance: Clearance = field(default_factory=Clearance)
    meander: MeanderParams = field(default_factory=lambda: MeanderParams(w=152))
    pad: PadParams = field(default_factory=lambda: PadParams(w=8, l=532))
    d_pad: int = 20
    coupling: common.InputCoupler = field(default_factory=lambda:
        common.InputCoupler(d=2, n=0, t_input=32, d_t1=48,
                            d_regular=62, t_regular_ground=6))
    allow_input_idc: bool = False
    d_pad_left: int|None = None  # Optional independent left cap-to-pad clearance.
    include_t2_pads: bool = True  # False removes only the two internal T2 ground pads.


SETTINGS = Settings()
CASE = 'left32_right6_dl128_dr122_w152_p532_i2_j48_r62'
RUN_INDUCTANCE = True
RUN_RESPONSE = True


def from_dict(q):
    return Settings(Cap(**q['cap']), Clearance(**q['clearance']),
                    MeanderParams(**q['meander']), PadParams(**q['pad']),
                    q['d_pad'], common.InputCoupler(**q['coupling']),q.get('allow_input_idc',False),
                    q.get('d_pad_left'), q.get('include_t2_pads', True))


def build(s=SETTINGS, *, isolation=False, input_probe=False):
    cap, gap, coil = s.cap, s.clearance, s.meander
    width = cap.t_left+gap.d_left+coil.w+gap.d_right+cap.t_right
    if s.coupling.n and not s.allow_input_idc:
        raise ValueError('An input IDC requires explicit authorization in this case')
    if any(v < 2 or v % 2 for v in [*asdict(cap).values(), *asdict(gap).values()]):
        raise ValueError('Cap dimensions/clearances must be positive and even')
    if width % 4:
        raise ValueError('Change an independent clearance by 4 um to retain assembly alignment')
    # The temporary symmetric scaffold is replaced below by the asymmetric
    # body. Its pad validation must not constrain the final independent gaps.
    scaffold_gap = min(s.d_pad, 2*np.floor(((width-2*cap.t_right-coil.w)/2-s.pad.w-2)/2))
    p = replace(common.T1, cap=replace(common.T1.cap, w=width,
                l=cap.l, t_top=cap.t_top, t_side=cap.t_right),
                meander=coil, pad=s.pad, d_pad=scaffold_gap)
    top, meta = common.build(p, isolation=isolation, coupler=s.coupling,input_probe=input_probe)
    old_ref = next(r for r in top.references
                   if np.allclose(r.origin, (meta['cell_centers'][1], 0)))
    top.remove(old_ref)
    left, right = -width/2, width/2
    x = left+cap.t_left+gap.d_left+coil.w/2
    north = -p.d_top
    underside = north-cap.t_top
    body = gdstk.Cell('ASYMMETRIC_T2')
    body.add(gdstk.Polygon([
        (left, north-cap.l), (left+cap.t_left, north-cap.l),
        (left+cap.t_left, underside), (right-cap.t_right, underside),
        (right-cap.t_right, north-cap.l), (right, north-cap.l),
        (right, north), (left, north)], layer=1))
    body.add(place(make_meander(coil, name='T2_COIL'), origin=(x, underside), anchor='top'))
    bottom = underside-coil.l
    pad = make_pad(s.pad, name='T2_PAD')
    left_pad_gap=s.d_pad if s.d_pad_left is None else s.d_pad_left
    if left_pad_gap<2 or left_pad_gap%2:
        raise ValueError('Left pad gap must be positive and even')
    pl, pr = left+cap.t_left+left_pad_gap, right-cap.t_right-s.d_pad
    if pl+s.pad.w >= x-coil.w/2 or pr-s.pad.w <= x+coil.w/2:
        raise ValueError('A ground pad touches the meander')
    if s.include_t2_pads:
        body.add(place(pad, origin=(pl, bottom), anchor='bottom_left'),
                 place(pad, origin=(pr, bottom), anchor='bottom_right'))
    if any(not np.allclose(q.points/2, np.round(q.points/2)) for q in body.get_polygons()):
        raise ValueError('New T2 body is off the 2 um grid')
    top.add(gdstk.Reference(body, origin=(meta['cell_centers'][1], 0)))
    meta['unit']['cap'] = dict(w=width, **asdict(cap))
    meta['unit']['d_pad'] = s.d_pad
    if not s.include_t2_pads:
        meta['unit']['include_t2_pads'] = False
    meta['unit']['meander']['x_offset_from_cap_center'] = x
    if s.d_pad_left is not None:
        meta['unit']['d_pad_left']=s.d_pad_left
    meta['layout_kind'] = ('asymmetric_plain_T2_input_probe' if input_probe else
                          ('asymmetric_plain_T2_calibration' if isolation else 'asymmetric_plain_T2_response'))
    meta['trial_settings'] = asdict(s)
    if s.include_t2_pads:
        meta['trial_settings'].pop('include_t2_pads')  # Preserve existing case metadata.
    if s.coupling.w_input_stub is None:
        meta['trial_settings']['coupling'].pop('w_input_stub')
    if not s.allow_input_idc:
        meta['trial_settings'].pop('allow_input_idc') # Preserve immutable plain cases.
    if s.d_pad_left is None:
        meta['trial_settings'].pop('d_pad_left')
    if s.coupling.n and (not isolation or input_probe):
        meta['layout_kind']=meta['layout_kind'].replace('plain','input_IDC')
        upper=-p.d_top-cap.t_top-2
        lower=upper-(s.coupling.n-1)*(s.coupling.t+s.coupling.d_finger)-s.coupling.t
        meta['IDC_clearance_um']=dict(meander=cap.t_left+gap.d_left,
            nearest_internal_ground_pad=cap.t_left+left_pad_gap,
            upper_ground=-upper,lower_ground=lower+618)
    meta['inductance_pin_x'] = meta['cell_centers'][1]+x
    meta['clearance_to_meander_um'] = dict(left=gap.d_left, right=gap.d_right)
    if not isolation:
        common.require_locked(top.get_polygons(), meta)
    return top, meta


def validate(path, meta):
    top, expected = build(from_dict(meta['trial_settings']), isolation=meta['isolation'],
                          input_probe=meta.get('input_probe',False))
    if meta != expected:
        raise ValueError('Case metadata changed')
    library = gdstk.read_gds(str(path)); tops = library.top_level()
    if len(tops) != 1 or tops[0].name != 'TOP' or not np.isclose(library.unit, 1e-6):
        raise ValueError('Wrong GDS units or top cell')
    polys = tops[0].get_polygons()
    if any(p.layer != 1 or p.datatype != 0 for p in polys):
        raise ValueError('Wrong metal layer/datatype')
    if gdstk.boolean(polys, top.get_polygons(), 'xor', precision=.001):
        raise ValueError('GDS differs from parameters')
    if len(gdstk.boolean(polys, [], 'or', precision=.001)) != meta['expected_metal_components']:
        raise ValueError('Unexpected metal connectivity')
    if not meta['isolation']:
        common.require_locked(polys, meta)
    return SimpleNamespace(n_cells=meta['parameters']['n_cells'], l_sim=1500, l_top_ground=448), polys


def prepare(name, s=SETTINGS, *, isolation=False,input_probe=False):
    stage = 'input_probe' if input_probe else ('Lcal' if isolation else 'response')
    native = ROOT/'sonnet_files/taper_exploration'/CAMPAIGN/name/stage
    out = ROOT/'output_data/taper_exploration'/CAMPAIGN/name/stage
    native.mkdir(parents=True, exist_ok=True); out.mkdir(parents=True, exist_ok=True)
    top, meta = build(s, isolation=isolation,input_probe=input_probe)
    gds = native/'geometry.gds'
    if gds.exists():
        if json.loads(gds.with_suffix('.json').read_text()) != meta:
            raise ValueError('Use a new CASE for changed geometry')
    else:
        library = gdstk.Library(unit=1e-6, precision=1e-9)
        library.add(top, *top.dependencies(True)); library.write_gds(str(gds))
        gds.with_suffix('.json').write_text(json.dumps(meta, indent=2)+'\n')
    validate(gds, meta)
    project = native/f"generated_ncell_{meta['parameters']['n_cells']}.sonx"
    if not project.exists():
        import_native(meta['parameters']['n_cells'], gds_path=gds, output_dir=native, geometry_validator=validate)
    # Calibration will add an analysis port; verify its unmodified metal before
    # adding it. For cached four-port files use the same metal XOR below.
    actual = [gdstk.Polygon(points(q)) for q in ET.parse(project).iter(tag('PlanarPolygon'))]
    expected = []
    for poly in top.get_polygons():
        xy=poly.points.copy(); xy[:,1]+=1052; xy=np.floor(xy/2+.5)*2; xy[:,1]=1500-xy[:,1]
        expected.append(gdstk.Polygon(xy))
    if gdstk.boolean(actual, expected, 'xor', precision=.001):
        raise ValueError('Native project metal differs from source')
    if not isolation:
        check_geometry(project, gds, meta, geometry_validator=validate)
        common.require_locked(actual, meta, native=True)
    fig, ax = plt.subplots(figsize=(17, 6))
    for poly in top.get_polygons():
        ax.add_patch(Polygon(poly.points, fc='#297b54', ec='none'))
    ax.set(xlim=(360, meta['w_sim']-360), ylim=(-650,45), aspect='equal',
           xlabel='x (um)', ylabel='y (um)', title=f'{name}: T2 - '+('L Calibration' if isolation and not input_probe else 'Input Coupling Probe' if input_probe else 'T1 - Regular'))
    fig.tight_layout(); fig.savefig(out/'geometry.png', dpi=170); plt.close(fig)
    (out/'specification.json').write_text(json.dumps(common.specification(), indent=2)+'\n')
    return project, out, meta


def measure(name, s=SETTINGS):
    """Explicit Sonnet capacitor sweep, cached by exact weak-fixture body."""
    from hashlib import sha256
    from geometry.taper import explicit_inductance as explicit
    project, out, meta = prepare(name, s, isolation=True)
    code_files=('geometry/taper/taper2_plain.py','geometry/taper/taper2_v3.py',
        'geometry/taper/three_tapers_v3.py','geometry/taper/three_tapers_v2.py',
        'geometry/unit_cell.py','geometry/meander.py','geometry/pad.py','geometry/prepare_sonnet.py')
    fingerprint=dict(unit=meta['unit'],template_sha256=common.W['digest'](ROOT/'sonnet_files/filter_v6.sonx'),
        T1_reference_sha256=common.W['digest'](common.REFERENCE_FILE),
        geometry_code_sha256={f:common.W['digest'](ROOT/f) for f in code_files})
    key=sha256(json.dumps(fingerprint,sort_keys=True).encode()).hexdigest()[:16]
    cache_native=ROOT/'sonnet_files/taper_exploration/explicit_L_cache'/key
    cache_out=ROOT/'output_data/taper_exploration/explicit_L_cache'/key
    source=cache_native/'source/generated_ncell_1.sonx'
    if not source.exists():
        source.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(project,source)
        for f in ('geometry.gds','geometry.json'):shutil.copy2(project.parent/f,source.parent/f)
        (source.parent/'fingerprint.json').write_text(json.dumps(fingerprint,indent=2)+'\n')
    if json.loads((source.parent/'fingerprint.json').read_text())!=json.loads(json.dumps(fingerprint)):
        raise ValueError('Explicit body-cache fingerprint changed')
    if gdstk.boolean(gdstk.read_gds(str(source.parent/'geometry.gds')).top_level()[0].get_polygons(),
        gdstk.read_gds(str(project.parent/'geometry.gds')).top_level()[0].get_polygons(),'xor',precision=.001):
        raise ValueError('Explicit cached fixture metal differs')
    manifest=cache_out/'manifest.json'
    if manifest.exists():caps=json.loads(manifest.read_text())['configuration']['C_add_fF']
    elif (out/'results.json').exists():
        # Cross-check an existing sparse calibration using exactly its old
        # added C values, but obtain every resonance from native Sonnet CAP.
        caps=[q['C_add_fF'] for q in json.loads((out/'results.json').read_text())['sweep']]
    else:caps=[250,350,500,700]
    result=explicit.measure(source,cache_native/'Lcal',cache_out,name=key,
        target_L_nH=common.TAPER_2.L0_nH,caps_fF=caps,metadata=meta)
    destination=project.parent.parent/'Lcal_explicit'
    destination.mkdir(exist_ok=True)
    for f in ('generated_ncell_1.sonx','generated_ncell_1.snp'):
        shutil.copy2(cache_native/'Lcal'/f,destination/f)
    explicit_out=out.parent/'Lcal_explicit';explicit_out.mkdir(exist_ok=True)
    for f in cache_out.iterdir():
        if f.is_file() and f.name!='completed.json':shutil.copy2(f,explicit_out/f.name)
    result.update(case=name,explicit_body_cache=key,source_result=str(cache_out/'results.json'))
    (explicit_out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    (explicit_out/'completed.json').write_text(json.dumps(dict(
        project_sha256=common.W['digest'](destination/'generated_ncell_1.sonx'),
        snp_sha256=common.W['digest'](destination/'generated_ncell_1.snp'),
        result_sha256=common.W['digest'](explicit_out/'results.json'),
        central_proof=str(cache_out/'completed.json'),central_proof_sha256=common.W['digest'](cache_out/'completed.json'),
        identical_weak_fixture_metal=True,body_fingerprint=fingerprint),indent=2)+'\n')
    return result


def probe_input(name,s,*,frequencies=(1,2,4)):
    """Three-frequency four-port test of the actual plain input electrode."""
    project,out,meta=prepare(name,s,isolation=True,input_probe=True)
    common.prepare_probe(project,meta,frequencies=frequencies)
    marker=out/'operator_completed.json'
    if marker.exists():
        q=json.loads(marker.read_text())
        if (common.W['digest'](project)!=q['project_sha256'] or
                common.W['digest'](project.with_suffix('.snp'))!=q['snp_sha256']):
            raise ValueError('Input probe cache changed')
        net=common.W['read_network'](project.with_suffix('.snp'),4)
    else:
        net=common.W['em'](project,4)
        marker.write_text(json.dumps(dict(project_sha256=common.W['digest'](project),
            snp_sha256=common.W['digest'](project.with_suffix('.snp'))),indent=2)+'\n')
    if not np.allclose(net.f,np.asarray(frequencies)*1e9,rtol=0,atol=1):
        raise ValueError('Probe frequencies changed; use a new case')
    y=common.AP['s_to_y'](*common.AP['reduce_pairs'](net.s))
    k,c,error=common.AP['extract_kc'](net.f,y)
    result=dict(scope='Isolated terminal-admittance proxy, not assembled S21 acceptance',
        input_Cg_proxy_fF=float(-c[2,0]),ground_C_proxy_fF=float(c[2].sum()),
        max_Y_residual_S=float(error),trial_settings=asdict(s),
        pointwise_input_Cg_fF=(-y[:,2,0].imag/(2*np.pi*net.f)*1e15).tolist(),
        frequencies_GHz=(net.f/1e9).tolist())
    if min(frequencies)>=4:
        result['scope']='Operating-band terminal coupling estimate; local K/C coefficients are not physical component acceptance'
        result['KC_local_fit_coefficients']=dict(input_Cg_fF=result.pop('input_Cg_proxy_fF'),
                                                ground_C_fF=result.pop('ground_C_proxy_fF'))
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result


def response(name, s, l):
    from geometry.taper.explicit_inductance import t1_inductance
    case = prepare(name, s)
    return common.response(name, None, l, prepared_case=case,fixed_T1_L_nH=t1_inductance())


def reuse_calibration(source, name, s):
    """Reuse EM only if both weak-fixture metal and complete T2 body match."""
    if source == name:
        return
    src = ROOT/'sonnet_files/taper_exploration'/CAMPAIGN/source/'Lcal'
    src_out = ROOT/'output_data/taper_exploration'/CAMPAIGN/source/'Lcal'
    old = json.loads((src/'geometry.json').read_text())
    top, meta = build(s, isolation=True)
    validate(src/'geometry.gds', old)
    if old['unit'] != meta['unit'] or gdstk.boolean(
            gdstk.read_gds(str(src/'geometry.gds')).top_level()[0].get_polygons(),
            top.get_polygons(), 'xor', precision=.001):
        raise ValueError('Calibration body or weak fixture changed; fresh EM is needed')
    marker=json.loads((src_out/'operator_completed.json').read_text())
    project=src/'generated_ncell_1.sonx'
    if (common.W['digest'](project) != marker['project_sha256'] or
            common.W['digest'](project.with_suffix('.snp')) != marker['snp_sha256']):
        raise ValueError('Source calibration changed')
    dst=ROOT/'sonnet_files/taper_exploration'/CAMPAIGN/name/'Lcal'
    dst_out=ROOT/'output_data/taper_exploration'/CAMPAIGN/name/'Lcal'
    if dst.exists() or dst_out.exists():
        raise ValueError('Destination calibration already exists')
    shutil.copytree(src,dst,ignore=shutil.ignore_patterns('sondata'))
    shutil.copytree(src_out,dst_out)
    (dst/'geometry.json').write_text(json.dumps(meta,indent=2)+'\n')
    result=json.loads((dst_out/'results.json').read_text())
    result.update(case=name,reused_from_case=source)
    (dst_out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    (dst_out/'reuse.json').write_text(json.dumps(dict(source_case=source,
        identical_body=True,identical_weak_fixture_metal=True,source_hashes_checked=True),indent=2)+'\n')
    prepare(name,s,isolation=True)


def finalize(name):
    """Archive only an actually passing native fixture and its source evidence."""
    native=ROOT/'sonnet_files/taper_exploration'/CAMPAIGN/name
    out=ROOT/'output_data/taper_exploration'/CAMPAIGN/name
    score=json.loads((out/'response/target_comparison.json').read_text())
    fitted=json.loads((out/'response/extraction/fit.json').read_text())
    measured=json.loads((out/'Lcal_explicit/results.json').read_text())
    if (not score['pass'] or score['complex_RMS'] > common.S21_COMPLEX_RMS_CUTOFF or
            not fitted['curve_pass'] or abs(measured['target_error_percent']) > 2):
        raise ValueError('Actual target response, circuit fit and measured L must pass')
    if score['specification'] != common.specification():
        raise ValueError('Superseded specification')
    meta=json.loads((native/'response/geometry.json').read_text())
    s=from_dict(meta['trial_settings'])
    if s.coupling.n and not s.allow_input_idc:
        raise ValueError('Unexpected IDC')
    check_geometry(native/'response/generated_ncell_3.sonx',native/'response/geometry.gds',meta,
                   geometry_validator=validate)
    common.require_reference()
    if measured.get('Python_frequency_approximation') is not False:
        raise ValueError('Final T2 requires explicit Sonnet capacitor calibration')
    for stage,n in (('Lcal_explicit',1),('response',3)):
        project=native/stage/f'generated_ncell_{n}.sonx'
        marker=json.loads((out/stage/('completed.json' if stage=='Lcal_explicit' else 'run_status.json')).read_text())
        if (common.W['digest'](project) != marker['project_sha256'] or
                common.W['digest'](project.with_suffix('.snp')) != marker['snp_sha256']):
            raise ValueError('EM sources changed')
        logs=list((out/stage).glob('cadd_*.log')) if stage=='Lcal_explicit' else [out/stage/'em.log']
        if not logs:raise ValueError('Missing explicit calibration logs')
        for log in logs:
            text=log.read_text()
            if not all(v in text for v in ('Em errors detected: 0','Em warnings detected: 0','Em simulation completed')):
                raise ValueError('Incomplete EM or native warnings/errors')
        if stage=='Lcal_explicit':
            if common.W['digest'](out/stage/'results.json')!=marker['result_sha256']:
                raise ValueError('Explicit inductance result changed')
            if common.W['digest'](Path(marker['central_proof']))!=marker['central_proof_sha256']:
                raise ValueError('Explicit calibration provenance changed')
    cal_meta=json.loads((native/'Lcal/geometry.json').read_text())
    if cal_meta['unit'] != meta['unit']:
        raise ValueError('Calibration body mismatch')
    check_geometry(native/'Lcal_explicit/generated_ncell_1.sonx',native/'Lcal/geometry.gds',cal_meta,
                   geometry_validator=validate)
    from geometry.taper.explicit_inductance import check_t1, T1_CHECK
    checked_t1=check_t1()
    if (not np.isclose(fitted['L_nH'][0],measured['L0_nH'],rtol=1e-12) or
            not np.isclose(fitted['L_nH'][1],checked_t1['L0_nH'],rtol=1e-12) or
            fitted['provenance']['sha256'] != score['source_sha256']):
        raise ValueError('Fitted source or fixed L differs')
    if list(ET.parse(native/'response/generated_ncell_3.sonx').iter(tag('Component'))):
        raise ValueError('Ideal circuit elements in production')
    kind='idc' if s.coupling.n else 'plain'
    destination=(ROOT/'sonnet_files/taper_exploration/final_designs'/f'new_T2_{kind}_{name}_final').resolve()
    if (not native.resolve().is_relative_to(ROOT) or not destination.is_relative_to(ROOT)
            or destination.exists()):
        raise ValueError('Unsafe or existing archive destination')
    destination.parent.mkdir(exist_ok=True)
    native.rename(destination)
    shutil.copytree(out/'Lcal_explicit',destination/'Lcal_explicit/evidence')
    for filename in ('s21_fit.png','s21_fit.pdf','geometry.png','target_comparison.json'):
        shutil.copy2(out/'response'/filename,destination/filename)
    for suffix in ('.gds','.json'):
        shutil.copy2(destination/'response'/('geometry'+suffix),destination/('geometry'+suffix))
    top,_=build(s)
    body=next(r.cell for r in top.references if r.cell.name=='ASYMMETRIC_T2')
    library=gdstk.Library(unit=1e-6,precision=1e-9)
    only=library.new_cell('TOP');only.add(gdstk.Reference(body))
    library.add(body,*body.dependencies(True));library.write_gds(str(destination/'t2_unit.gds'))
    selection=dict(case=name,campaign=CAMPAIGN,specification=common.specification(),
        T1_accepted=True,T2_accepted=True,T3_accepted=False,settings=asdict(s),geometry=meta,
        isolated_L=measured,T1_explicit_L=checked_t1,
        T1_explicit_L_result_sha256=common.W['digest'](T1_CHECK),
        fixed_target_comparison=score,heuristic_fit=fitted,
        T1_reference=str(common.REFERENCE_FILE),
        project=str(destination/'response/generated_ncell_3.sonx'),plot=str(destination/'s21_fit.png'),
        gds=str(destination/'geometry.gds'),unit_gds=str(destination/'t2_unit.gds'),
        audit=dict(no_IDC=not bool(s.coupling.n),input_only_IDC=bool(s.coupling.n),production_ideal_components=0,mesh_um=[2,2],
            T1_regular_output_shapes_unchanged=True,source_hashes_checked=True,
            native_GDS_XOR_empty=True,calibration_body_matches=True,native_errors=0,native_warnings=0))
    selection['code_sha256']={p:common.W['digest'](ROOT/p) for p in
        ('geometry/taper/taper2_plain.py','geometry/taper/taper2_idc.py',
         'geometry/taper/taper2_asymmetric_preview.py','geometry/taper/taper2_v3.py',
            'geometry/taper/progressive_fitting.py','geometry/taper/explicit_inductance.py',
            'extract_generated_inductance.py','geometry/specifications.py')}
    for path in (destination/'selection.json',ROOT/f'geometry/taper/selected_taper2_{kind}.json',
                 ROOT/'geometry/taper/three_tapers_status.json'):
        path.write_text(json.dumps(selection,indent=2)+'\n')
    (destination/'README.md').write_text(
        f'# Accepted asymmetric {kind}-capacitor T2 fixture\n\n'
        'Design: `response/generated_ncell_3.sonx`. Shareable comparison: `s21_fit.png`.\n'
        f"Fixed-target complex S21 RMS: {score['complex_RMS']:.8f} (limit0.025).\n"
        f"Independent T2 L: {measured['L0_nH']:.8f} nH ({measured['target_error_percent']:+.3f}% versus4.33).\n"
        'The same 7501 frequencies span4.65-5.40GHz. All nominal electrical targets are fixed; only one matched phase delay is adjusted.\n'
        'The separate heuristic fit varies all taper capacitances and delay while holding every L and regular extraction constant fixed.\n'
        'T1, regular and output component shapes are preserved. Assembly gap/strip parameters are explicit in selection.json.\n'
        f'Input-only IDC: {bool(s.coupling.n)}; flat top ground; mesh2x2um. Input and T2-to-T1 gaps contain no ground strip.\n'
        'Independent left/right clearances and cap-leg thicknesses are in settings. Cap width and coil offset are derived.\n'
        'Inductance uses explicit native Sonnet CAP sweeps in Lcal_explicit, with cached subsequent capacitor changes.\n'
        'Lcal holds the matching source geometry; Lcal_explicit/evidence contains all capacitor projects, exports, logs and regression results.\n'
        't2_unit.gds is the body only; geometry.gds contains the complete simulated fixture.\n'
        'This validates the T2-T1-regular fixture, not T3 or the complete filter.\n',encoding='utf-8')
    print('FINAL T2',destination,flush=True)
    return selection


def main():
    selected=ROOT/'geometry/taper/selected_taper2_plain.json'
    if selected.exists():
        saved=json.loads(selected.read_text())
        if saved['case']==CASE and saved['settings']==asdict(SETTINGS):
            print('Already verified:',saved['project']);return
    measured = measure(CASE, SETTINGS) if RUN_INDUCTANCE else json.loads(
        (ROOT/'output_data/taper_exploration'/CAMPAIGN/CASE/'Lcal_explicit/results.json').read_text())
    print('ISOLATED L', measured['L0_nH'], 'ERROR %', measured['target_error_percent'], flush=True)
    if abs(measured['target_error_percent']) > 2:
        raise ValueError('Tune isolated inductance before full response')
    if RUN_RESPONSE:
        score, fitted = response(CASE, SETTINGS, measured['L0_nH'])
        print(json.dumps(dict(target=score, heuristic_fit=fitted), indent=2), flush=True)
        if score['pass'] and fitted['curve_pass']:
            finalize(CASE)


if __name__ == '__main__':
    main()
