r"""Generate the asymmetric, plain-capacitor T2 for visual review only.

Run: .venv\Scripts\python.exe -B -m geometry.taper.taper2_asymmetric_preview
This script imports GDS into a native Sonnet project. It never runs EM or fits.
The accepted T1 and regular component shapes are checked against their source.
"""
from dataclasses import dataclass, asdict, replace
from types import SimpleNamespace
import json

import gdstk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

from geometry.meander import make_meander
from geometry.pad import make_pad
from geometry.placement import place
from geometry.taper import taper2_v3 as fixture
from geometry.taper.audit import check_geometry
from geometry.prepare_sonnet import main as import_native


@dataclass(frozen=True)
class Cap:
    """Independent side thicknesses; dimensions in microns."""
    l: int = 590
    t_top: int = 14
    t_left: int = 32
    t_right: int = 6


@dataclass(frozen=True)
class Clearance:
    """Edge-to-edge distances from each cap leg to the meander envelope.

    These exclude the cap-leg thickness and include the ground pad in the gap.
    Change either by 4 um to preserve the existing assembly's grid alignment.
    """
    d_left: int = 128
    d_right: int = 122


CAP = Cap()
CLEARANCE = Clearance()
CASE = 'plain_asym_left32_right6_dl128_dr122_w152_p532_preview'
COUPLER = fixture.InputCoupler(d=2, n=0, t_input=32, d_t1=48,
                              d_regular=62, t_regular_ground=6)


def build(cap=CAP, coupler=COUPLER, clearance=CLEARANCE):
    # Use the proven assembly for feeds, planes, T1, regular and output cells.
    # Its temporary symmetric T2 is replaced below by the asymmetric component.
    meander_w = 152
    width = cap.t_left+clearance.d_left+meander_w+clearance.d_right+cap.t_right
    p = replace(fixture.T1, cap=replace(fixture.T1.cap, w=width,
                l=cap.l, t_top=cap.t_top, t_side=cap.t_right),
                meander=replace(fixture.T1.meander, w=152),
                pad=replace(fixture.T1.pad, l=532))
    if (any(v < 2 or v % 2 for v in [*asdict(cap).values(), *asdict(clearance).values()]) or
            width % 4 or cap.t_top >= cap.l):
        raise ValueError('Cap dimensions must fit and remain on the 2 um grid')
    if coupler.n:
        raise ValueError('This review geometry must have no IDC')
    top, meta = fixture.build(p, coupler=coupler)
    ref = next(r for r in top.references if r.cell.name == 'STAGE_1')
    top.remove(ref)

    # One continuous cap with separately controlled left and right legs.
    left, right = -width/2, width/2
    # The cap's width is derived, not used to constrain both gaps together.
    # This offset keeps the requested right clearance when only d_left changes.
    meander_x = left+cap.t_left+clearance.d_left+meander_w/2
    north = -p.d_top
    underside = north-cap.t_top
    bottom = north-cap.l
    cap_polygon = gdstk.Polygon([
        (left, bottom), (left+cap.t_left, bottom),
        (left+cap.t_left, underside), (right-cap.t_right, underside),
        (right-cap.t_right, bottom), (right, bottom),
        (right, north), (left, north)], layer=1, datatype=0)
    body = gdstk.Cell('ASYMMETRIC_T2')
    body.add(cap_polygon)
    body.add(place(make_meander(p.meander, name='ASYM_T2_MEANDER'),
                   origin=(meander_x, underside), anchor='top'))
    ground_bottom = underside-p.meander.l
    pad = make_pad(p.pad, name='ASYM_T2_PAD')
    pad_left = left+cap.t_left+p.d_pad
    pad_right = right-cap.t_right-p.d_pad
    if (pad_left+p.pad.w >= meander_x-p.meander.w/2 or
            pad_right-p.pad.w <= meander_x+p.meander.w/2):
        raise ValueError('A pad touches the meander')
    body.add(place(pad, origin=(pad_left, ground_bottom), anchor='bottom_left'),
             place(pad, origin=(pad_right, ground_bottom), anchor='bottom_right'))
    if any(not np.allclose(q.points/2, np.round(q.points/2)) for q in body.get_polygons()):
        raise ValueError('Asymmetric T2 body is off the 2 um grid')
    top.add(gdstk.Reference(body, origin=(meta['cell_centers'][1], 0)))
    meta['layout_kind'] = 'asymmetric_plain_T2_review_only'
    meta['unit']['cap'] = dict(w=width, **asdict(cap))
    meta['unit']['meander']['x_offset_from_cap_center'] = meander_x
    meta['review_cap'] = asdict(cap)
    meta['review_clearance'] = asdict(clearance)
    meta['simulation_run'] = False
    meta['inductance_calibration_valid_for_this_geometry'] = False
    meta['clearance_to_meander_um'] = {
        'left_cap_leg': clearance.d_left,
        'right_cap_leg': clearance.d_right}
    fixture.require_locked(top.get_polygons(), meta)
    return top, meta


def validate(path, metadata):
    top, expected = build(Cap(**metadata['review_cap']),
                          fixture.InputCoupler(**metadata['coupler']),
                          Clearance(**metadata['review_clearance']))
    if metadata != expected:
        raise ValueError('Review metadata differs from its parameters')
    library = gdstk.read_gds(str(path))
    cells = library.top_level()
    if (len(cells) != 1 or cells[0].name != 'TOP' or
            not np.isclose(library.unit, 1e-6)):
        raise ValueError('Unexpected GDS units or top cell')
    polygons = cells[0].get_polygons()
    if any(p.layer != 1 or p.datatype != 0 for p in polygons):
        raise ValueError('Wrong layer/datatype')
    # The preserved regular-cell center leads have the original 1 um source
    # offsets. The existing importer snaps these identically to the accepted
    # project; check_geometry verifies the full native 2 um grid and metal XOR.
    if gdstk.boolean(polygons, top.get_polygons(), 'xor', precision=.001):
        raise ValueError('GDS differs from the review geometry')
    if len(gdstk.boolean(polygons, [], 'or', precision=.001)) != 2:
        raise ValueError('Input capacitor shorted or metal disconnected')
    fixture.require_locked(polygons, metadata)
    return SimpleNamespace(n_cells=3, l_sim=1500, l_top_ground=448), polygons


def main():
    top, meta = build()
    root = fixture.ROOT
    native = root/'sonnet_files/taper_exploration'/fixture.CAMPAIGN/CASE
    output = root/'output_data/taper_exploration'/fixture.CAMPAIGN/CASE
    native.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=True)
    gds = native/'geometry.gds'
    if gds.exists():
        validate(gds, json.loads(gds.with_suffix('.json').read_text()))
    else:
        library = gdstk.Library(unit=1e-6, precision=1e-9)
        library.add(top, *top.dependencies(True))
        library.write_gds(str(gds))
        gds.with_suffix('.json').write_text(json.dumps(meta, indent=2)+'\n')
    project = native/'generated_ncell_3.sonx'
    if not project.exists():
        import_native(3, gds_path=gds, output_dir=native, geometry_validator=validate)
    check_geometry(project, gds, meta, geometry_validator=validate)
    for name, limits in (
            ('geometry', (360, meta['w_sim']-360)),
            ('t2_closeup', (meta['cell_lefts'][1]-90,
                            meta['cell_lefts'][2]+80))):
        fig, ax = plt.subplots(figsize=(16, 7))
        for polygon in top.get_polygons():
            ax.add_patch(Polygon(polygon.points, fc='#297b54', ec='none'))
        ax.set(xlim=limits, ylim=(-650, 45), aspect='equal',
               xlabel='x (um)', ylabel='y (um)',
               title='T2 Plain Capacitor: 32 um Left Leg / 6 um Right Leg — Preview Only')
        fig.tight_layout()
        fig.savefig(output/f'{name}.png', dpi=170)
        plt.close(fig)
    report = dict(preview_only=True, simulation_run=False,
                  native_GDS_XOR_empty=True, grid_um=[2, 2],
                  T1_regular_output_shapes_unchanged=True,
                  project=str(project), gds=str(gds),
                  preview=str(output/'geometry.png'),
                  closeup=str(output/'t2_closeup.png'), geometry=meta)
    (native/'preview.json').write_text(json.dumps(report, indent=2)+'\n')
    (native/'README.md').write_text(
        '# Asymmetric plain T2: visual review only\n\n'
        'Open `generated_ncell_3.sonx`. No simulation or parameter extraction has been run.\n'
        f"T2 cap: derived width {meta['unit']['cap']['w']}, length {CAP.l}, top {CAP.t_top}, left leg {CAP.t_left}, right leg {CAP.t_right} um.\n"
        'Input electrode: 32 um thick, 2 um gap, no IDC. T2-to-T1 plain gap: 48 um.\n'
        'Meander: 152 x 596 um, 72 runs, 4 um trace / 4 um gaps.\n'
        f'Independent cap-leg/meander clearances: left {CLEARANCE.d_left} um, right {CLEARANCE.d_right} um.\n'
        'Edit CLEARANCE.d_left or d_right; cap width and meander offset are derived automatically.\n'
        'Changing only d_left leaves the right cap leg, right pad and meander unchanged relative to one another.\n'
        'Pads: 8 x 532 um, 20 um from their respective cap legs.\n'
        'T1/regular/output component shapes remain unchanged. Centered shared strip: 6 um in a 62 um gap.\n'
        'Upper ground remains flat. Native 2 x 2 um mesh, source materials and ports are preserved.\n'
        'Old measured L/capacitances do not apply automatically to this modified body.\n'
        'Edit `geometry/taper/taper2_asymmetric_preview.py` to regenerate a fresh named preview.\n',
        encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
