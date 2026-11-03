r"""One controlled follow-up to the completed plate276 baseline (2 um grid).

Run from the project root with .venv\Scripts\python.exe -B -m
geometry.taper.resume_second_taper. Every solver inherits the four-thread cap.
Completed inductance records are reused; each geometry has its own directory.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '4'

from dataclasses import replace
from pathlib import Path
import json
from geometry.taper.second_layout import from_dict, build
from geometry.taper.experiment import measure_unit
from geometry.taper.second_experiment import run

ROOT = Path(__file__).resolve().parents[2]
BASELINE = 't2_plate276_slot264_tooth96_l400'
CASE = 't2_left540_ground2_t1r300_d48'
T1_CALIBRATION = 't1_bar10_r300_slot264'
RUN_SIMULATION = False  # Work is paused for extraction-model validation.


def main():
    if not RUN_SIMULATION:
        print('No EM launched. See geometry/taper/taper2_status.json and TAPER_2_README.md.')
        return
    metadata = json.loads((ROOT / 'geometry/taper/output' / f'{BASELINE}.json').read_text())
    baseline = from_dict(metadata['second'])
    params = replace(baseline,
        taper=replace(baseline.taper, l_left=540),
        taper1=replace(baseline.taper1, l_right=300),
        w_coupling_teeth=112, d_t1_regular=48, t_junction_ground=2)
    build(params)  # Check clearances and grid before any simulation.
    reference1 = ROOT / 'output_data/taper_exploration' / (T1_CALIBRATION + '_isolation/results.json')
    if not reference1.exists():
        reference1, _ = measure_unit(T1_CALIBRATION, params.taper1,
            sweep_range=(3.8, 6.4), w_top_slot=params.w_t1_ground_slot,
            l_top_slot=params.l_t1_ground_slot)
    result = run(CASE, params, inductance_t1_reference=reference1)
    print('TARGET PASS:', result['passes_2_percent'], flush=True)


if __name__ == '__main__':
    main()
