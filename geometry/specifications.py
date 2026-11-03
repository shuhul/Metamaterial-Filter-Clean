"""Central electrical targets. L in nH, capacitances in fF, geometry in um.

Cg_input belongs to the named stage's INPUT: input -> t3 -> t2 -> t1 -> regular.
These are requested specifications, never values to pin during extraction.
"""
from dataclasses import dataclass, asdict


@dataclass(frozen=True)
class CellTarget:
    L0_nH: float
    C0_fF: float
    Cg_input_fF: float


UNIT_CELL = CellTarget(4.9194, 190.4, 8.2865)
PREVIOUS_TAPER_1 = CellTarget(6.06, 100.0, 56.7)
PREVIOUS_TAPER_2 = CellTarget(0.74, 100.0, 459.0)
ALTERNATIVE_TAPER_2 = CellTarget(1.0, 10.0, 264.0)
PREVIOUS_THREE_TAPER_TARGETS = {'taper_1': CellTarget(5.91,152.,12.3),
    'taper_2': CellTarget(4.42,155.,69.7), 'taper_3': CellTarget(2.44,116.,500.)}
SPECIFICATION_VERSION = 'three_tapers_v3_2026-10-26'
TAPER_1 = CellTarget(5.97, 148.0, 15.3)
TAPER_2 = CellTarget(4.33, 149.0, 153.0)
TAPER_3 = CellTarget(2.77, 251.0, 261.0)
TARGETS = {'unit_cell': UNIT_CELL, 'taper_1': TAPER_1,
           'taper_2': TAPER_2, 'taper_3': TAPER_3}
TOLERANCE_PERCENT = 2.0
GRID_UM = 2.0
# Current T2 geometry objective: fitted new capacitances within5% of targets.
# Kept separate from historical specification() records and L tolerances.
TAPER_2_CAPACITANCE_TOLERANCE_PERCENT = 5.0
S21_COMPLEX_RMS_CUTOFF = 0.025
TARGET_S21_FREQUENCY_RANGE_GHZ = (4.65, 5.40)


def specification():
    return {'version': SPECIFICATION_VERSION,
            'units': {'L': 'nH', 'C': 'fF', 'geometry': 'um'},
            'targets': {name: asdict(target) for name, target in TARGETS.items()},
            'tolerance_percent': TOLERANCE_PERCENT, 'grid_um': GRID_UM,
            'curve_acceptance': {'metric': 'complex_S21_RMS', 'maximum': S21_COMPLEX_RMS_CUTOFF},
            'target_response_acceptance': {
                'metric': 'complex_S21_RMS_against_fixed_target_circuit',
                'maximum': S21_COMPLEX_RMS_CUTOFF,
                'frequency_range_GHz': list(TARGET_S21_FREQUENCY_RANGE_GHZ),
                'electrical_parameters': 'All specified targets fixed',
                'nuisance_parameter': 'One matched fixture delay only',
                'extracted_capacitances': 'Geometry heuristics; fitted-curve RMS is diagnostic'},
            'order': ['input', 'taper_3', 'taper_2', 'taper_1', 'unit_cell', 'output']}
