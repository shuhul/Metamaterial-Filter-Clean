r"""Editable reduced-assembly T3 tuning; full verification is a separate gate.

Run from root: .venv\Scripts\python.exe -B -m geometry.taper.run_taper3_screening
Keep each changed trial under a new CASE. No selected older design is modified.
"""
import json
from dataclasses import replace
from geometry.taper.taper3 import NATIVE, OUT, from_dict, measure, fit_response
from geometry.taper.taper3_screening import run

BASE_CASE = 'w130_t6g2n66_dl342_dr340_pad80x560_mc24_i40_j28n18'
BASE = from_dict(json.loads((NATIVE/BASE_CASE/'response/geometry.json').read_text())['t3_settings'])
CASE = 'w130_t6g2n66_dl342_dr340_pad80x560_mc24_i70ni40_j28n18'
PARAMS = replace(BASE, d_input=70, n_input=40)
CAPACITORS_FF = [450, 600, 750, 900]
RUN_T3_T2 = True
RUN_T3_T2_T1 = False
RUN_FULL_VERIFICATION = False


def main():
    measured = measure(CASE, PARAMS, caps_fF=CAPACITORS_FF)
    if abs(measured['target_error_percent'])>2:
        raise ValueError('Retune L before a response solve')
    for branches, enabled in [(2,RUN_T3_T2),(3,RUN_T3_T2_T1)]:
        if enabled:
            run(CASE, PARAMS, measured['L0_nH'], branches)
    if RUN_FULL_VERIFICATION:
        fit_response(CASE, PARAMS, measured['L0_nH'])


if __name__ == '__main__':
    main()
