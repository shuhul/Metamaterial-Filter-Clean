"""Preflight checks for the supported 2 um-grid Sonnet workflow."""

from dataclasses import asdict
import math

import gdstk
import numpy as np

if __package__:
    from .ncell import NCellParams, make_ncell
    from .unit_cell import UnitCellParams
    from .capacitor import CapParams
    from .meander import MeanderParams
    from .pad import PadParams
else:
    from ncell import NCellParams, make_ncell
    from unit_cell import UnitCellParams
    from capacitor import CapParams
    from meander import MeanderParams
    from pad import PadParams


def require(condition, message):
    if not condition:
        raise ValueError(message)


def metadata_params(metadata):
    data = dict(metadata['parameters'])
    unit = dict(data.pop('unit'))
    unit['cap'] = CapParams(**unit['cap'])
    unit['meander'] = MeanderParams(**unit['meander'])
    unit['pad'] = PadParams(**unit['pad'])
    return NCellParams(unit=UnitCellParams(**unit), **data)


def validate_grid_params(p):
    """Reject changes that can collapse narrow features when snapping vertices."""
    def visit(data, prefix=''):
        for key, value in data.items():
            if isinstance(value, dict):
                visit(value, prefix+key+'.')
            elif key not in ('n', 'n_cells'):
                require(math.isclose(value/2, round(value/2), abs_tol=1e-9),
                        f'{prefix}{key} must be a multiple of 2 um for Sonnet import')
    visit(asdict(p))
    require(p.d_cell-p.t_ground >= 4, 'Ground strips need at least 2 um clearance on each side')
    side_clearance = p.unit.cap.w/2-p.unit.cap.t_side-p.unit.d_pad-p.unit.pad.w-p.unit.meander.w/2
    require(side_clearance >= 2, 'Pads need at least 2 um clearance to the meander')
    require(p.unit.cap.t_top+p.unit.meander.l-p.unit.cap.l >= 2,
            'Cap legs need at least 2 um clearance to bottom ground')
    require(p.w_side > 400, '400 um reference planes must end inside the uniform feed region')


def validate_gds_metadata(path, metadata):
    p = metadata_params(metadata)
    validate_grid_params(p)
    require(metadata['w_sim'] == p.w_sim and metadata['l_total'] == p.l_total,
            'Derived dimensions in metadata are stale')
    require(np.allclose(metadata['cell_centers'], p.cell_centers), 'Stale cell centers')
    require(np.allclose(metadata['intended_box'], p.box), 'Stale box metadata')
    lib = gdstk.read_gds(str(path))
    require(lib.unit == 1e-6 and lib.precision == 1e-9, 'GDS must use um units and 1 nm precision')
    tops = lib.top_level()
    require(len(tops) == 1 and tops[0].name == 'TOP', 'Expected one top-level TOP cell')
    polygons = tops[0].get_polygons()
    require(all((q.layer,q.datatype)==(1,0) for q in polygons), 'Unexpected GDS layer/datatype')
    rebuilt = make_ncell(p).get_polygons()
    require(not gdstk.boolean(polygons,rebuilt,'xor',precision=.001),
            'GDS does not match its parameter metadata; regenerate both together')
    return p, polygons
