"""Reduced T3+T2 and T3+T2+T1 screening; final acceptance always uses full EM.

Retain the full box/planes and all metal upstream of the last cap edge. Replace
omitted cells by a feed directly connected to the last retained cap top. This
changes loading and removes mutual coupling: extracted values are provisional.
"""
import json
import time
from copy import deepcopy
from types import SimpleNamespace
import gdstk
import numpy as np
from geometry.taper import taper3 as work
from geometry.taper.progressive_fitting import fit


def build(params, branches):
    if branches not in (2, 3):
        raise ValueError('Screen exactly T3+T2 or T3+T2+T1')
    full, original = work.build(params, stage='response')
    old_index = branches-1  # Original source indices1=T2,2=T1.
    center = work.SOURCE_META['cell_centers'][old_index]
    left = work.SOURCE_META['cell_lefts'][old_index]
    edge = 2*center-left+original['downstream_shift_um']
    center += original['downstream_shift_um']
    width = original['w_sim']
    roi = gdstk.rectangle((0, -1052), (edge, 448))
    retained = gdstk.boolean(full.get_polygons(), [roi], 'and', precision=.001, layer=1)
    parts = retained+[
        gdstk.rectangle((edge, 0), (width, 448)),
        gdstk.rectangle((edge, -1052), (width, -618)),
        gdstk.rectangle((width-412, -1052), (width, -24)),
        gdstk.rectangle((center, -18), (width, -8))]
    top = gdstk.Cell('TOP')
    top.add(*gdstk.boolean(parts, [], 'or', precision=.001, layer=1))
    if gdstk.boolean(retained, gdstk.boolean(top.get_polygons(), [roi], 'and', precision=.001),
                     'xor', precision=.001):
        raise ValueError('Retained upstream geometry changed')
    meta = deepcopy(original)
    meta.update(layout_kind=f'T3_screen_{branches}_branches', stage='screen',
        screening_branches=branches, clip_boundary_um=edge,
        retained_upstream_xor_empty=True, final_acceptance_permitted=False,
        expected_metal_components=3,
        terminal_connection='50ohm box port directly wired to last retained cap top',
        limitations='Omitted-cell mutual coupling/loading and downstream ground strips absent. '
                    'Top/bottom planes remain grounded at box walls; verify against full assembly.')
    meta['parameters']['n_cells'] = branches
    return top, meta


def validate(path, meta):
    expected, rebuilt = build(work.from_dict(meta['t3_settings']), meta['screening_branches'])
    if meta != rebuilt:
        raise ValueError('Reduced geometry metadata differs')
    lib = gdstk.read_gds(str(path)); tops = lib.top_level()
    if len(tops)!=1 or tops[0].name!='TOP' or not np.isclose(lib.unit, 1e-6):
        raise ValueError('Invalid reduced GDS top/units')
    polys = tops[0].get_polygons()
    if any(p.layer!=1 or p.datatype!=0 for p in polys):
        raise ValueError('Wrong metal mapping')
    if gdstk.boolean(polys, expected.get_polygons(), 'xor', precision=.001):
        raise ValueError('Reduced GDS differs from definition')
    if len(gdstk.boolean(polys, [], 'or', precision=.001))!=meta['expected_metal_components']:
        raise ValueError('Unexpected reduced-layout connectivity')
    return SimpleNamespace(n_cells=meta['screening_branches'],l_sim=1500,l_top_ground=448),polys


def prepare(case, params, branches):
    name = 'screen_T3_T2'+('_T1' if branches==3 else '')
    native, out = work.NATIVE/case/name, work.OUT/case/name
    native.mkdir(parents=True, exist_ok=True); out.mkdir(parents=True, exist_ok=True)
    top, meta = build(params, branches)
    gds = native/'geometry.gds'; saved = gds.with_suffix('.json')
    if saved.exists() and json.loads(saved.read_text())!=meta:
        raise ValueError('Use a new case for changed reduced geometry')
    if not gds.exists():
        lib = gdstk.Library(unit=1e-6, precision=1e-9); lib.add(top); lib.write_gds(str(gds))
        saved.write_text(json.dumps(meta, indent=2)+'\n')
    project = native/f'generated_ncell_{branches}.sonx'
    if not project.exists():
        work.import_native(branches, gds_path=gds, output_dir=native, geometry_validator=validate)
    work.check_geometry(project, gds, meta, geometry_validator=validate)
    fig, ax = work.plt.subplots(figsize=(18,6))
    for poly in top.get_polygons():
        ax.add_patch(work.Polygon(poly.points,fc='#297b54',ec='none'))
    ax.set(xlim=(360,meta['w_sim']-360),ylim=(-650,45),aspect='equal',
        xlabel='x (µm)',ylabel='y (µm)',title=name+': screening only; terminal taper directly connected to output')
    fig.tight_layout(); fig.savefig(out/'geometry.png',dpi=170);work.plt.close(fig)
    return project,out


def run(case, params, l0, branches):
    started = time.perf_counter()
    project,out = prepare(case,params,branches)
    measured = json.loads((work.OUT/case/'inductance.json').read_text())
    calmeta = json.loads((work.NATIVE/case/'Lcal/geometry.json').read_text())
    ownmeta = json.loads((project.parent/'geometry.json').read_text())
    if abs(l0-measured['L0_nH'])>1e-10 or abs(measured['target_error_percent'])>2 or calmeta['unit']!=ownmeta['unit']:
        raise ValueError('Reduced geometry must use its independently calibrated body L')
    status = out/'run_status.json'
    if status.exists() and json.loads(status.read_text()).get('status')=='complete':
        previous = json.loads(status.read_text())
        if work.digest(project)!=previous['project_sha256'] or work.digest(project.with_suffix('.snp'))!=previous['snp_sha256']:
            raise ValueError('Completed reduced response changed')
    else:
        work.configure(project,4.0,6.0,10000);work.simulate(project,out)
    ls = [l0]+work.REFERENCE_FIT['L_nH'][:branches-1]
    cs = [251]+work.REFERENCE_FIT['C_fF'][:branches-1]
    gs = [261,153]+work.REFERENCE_FIT['G_fF'][1:branches-1]
    # Fit the same band used for full-assembly acceptance; retain a broader-band
    # diagnostic to expose frequency-window/model dependence.
    from circuit_common import read_snp
    net = read_snp(project.with_suffix('.snp'))
    narrow = net[(net.f>=4.65e9)&(net.f<=5.4e9)]
    narrow.write_touchstone(filename=str(out/'same_band'),form='ri')
    band_file = out/'same_band.snp';band_file.write_bytes((out/'same_band.s2p').read_bytes())
    seeds=[[100,261,153],[251,261,153],[400,120,150],[600,300,300]]
    kwargs=dict(capacitance_seeds=seeds,title=f'{branches}-Branch T3 Screening (Not Final Verification)')
    same = fit(band_file,out/'fit_same_band',ls,cs,gs,[('C',0),('G',0),('G',1)],**kwargs)
    broad = fit(project.with_suffix('.snp'),out/'fit_broad_band',ls,cs,gs,[('C',0),('G',0),('G',1)],**kwargs)
    report=dict(case=case,branches=branches,project=str(project),screening_only=True,
        final_acceptance_permitted=False,all_older_values_fixed=True,
        same_band_values_fF=[same['C_fF'][0],same['G_fF'][0],same['G_fF'][1]],
        same_band_complex_RMS=same['complex_RMS'],same_band_active_bounds=same['active_bounds'],
        broad_band_values_fF=[broad['C_fF'][0],broad['G_fF'][0],broad['G_fF'][1]],
        broad_band_complex_RMS=broad['complex_RMS'],
        preparation_EM_and_fitting_elapsed_seconds=time.perf_counter()-started,
        elapsed_includes_solver_queue_wait=True,
        native_solver_elapsed_seconds=json.loads((out/'em.log.resources.json').read_text())['elapsed_seconds'],
        limitations='Different terminal loading; omitted mutual coupling. Transfer to full geometry must be measured.')
    (out/'result.json').write_text(json.dumps(report,indent=2)+'\n')
    print('SCREENING RESULT',json.dumps(report),flush=True)
    return report
