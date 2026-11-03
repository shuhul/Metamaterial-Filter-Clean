"""Edit grouped T3 settings/flags here. Run as a module from the project root.

Use a fresh CASE for changed parameters. Selected T2/T1 geometry is immutable.
The direct-feed stage has no input IDC. Only the final response tests input Cg.
"""
from dataclasses import replace
from geometry.pad import PadParams
from geometry.taper.taper3 import Params, OUT, ROOT, asdict, measure, fit_response, fit_older_within_one_percent
import json

CASE = 'regular_w80_n72_dl476_dr500_padL460_R24_top2_i84ni33_stub200_FIXED'
BASE = Params()
D_TOP = 2
PARAMS = replace(BASE, d_top=D_TOP, pad_right=PadParams(w=24,l=594), body=replace(BASE.body,
    cap=replace(BASE.body.cap,l=598-D_TOP,t_top=22-D_TOP),
    coupling=replace(BASE.body.coupling,w_input_stub=200),
    clearance=replace(BASE.body.clearance,d_left=476,d_right=500),d_pad=464,d_pad_left=4,
    meander=replace(BASE.body.meander,w=80,n=72,t=4,d=4,l_top=10),
    pad=replace(BASE.body.pad,w=460,l=594)),
    d_input=84,n_input=33,d_t2=22,n_t2=70,t_finger=4,d_finger=4)
FIX_T2_JUNCTION = True  #156.60039987587515fF from the accepted T2 fit
L_CAPS_FF = [450,600,750,900]
GENERATE_PREVIEW_ONLY = False
CALIBRATE_INDUCTANCE = True
SIMULATE_DIRECT_FEED = False
SIMULATE_FULL_RESPONSE = True


def main():
    progress_file=ROOT/'geometry/taper/taper3_progress.json'
    progress=json.loads(progress_file.read_text())
    progress.update(active_case=CASE,accepted=False,state='preparing_preview' if GENERATE_PREVIEW_ONLY else 'calibrating_inductance')
    progress_file.write_text(json.dumps(progress,indent=2)+'\n')
    if GENERATE_PREVIEW_ONLY:
        from geometry.taper.taper3 import prepare
        print(prepare(CASE, PARAMS, stage="response")[0])
        return
    if CALIBRATE_INDUCTANCE:
        measured = measure(CASE, PARAMS, caps_fF=L_CAPS_FF)
    else:
        measured = json.loads((OUT/CASE/'inductance.json').read_text())
    progress.update(active_L_nH=measured['L0_nH'],active_L_error_percent=measured['target_error_percent'],
                    state='L_outside_tolerance' if abs(measured['target_error_percent'])>2 else 'L_passed_preparing_response')
    progress_file.write_text(json.dumps(progress,indent=2)+'\n')
    if abs(measured['target_error_percent']) > 2:
        raise ValueError('Tune T3 L to within2% before spending a full-response solve')
    if SIMULATE_DIRECT_FEED:
        direct = fit_response(CASE, PARAMS, measured['L0_nH'], stage='direct')
        if not direct['fit_pass']:
            fit_older_within_one_percent(CASE, measured['L0_nH'], stage='direct')
            raise ValueError('Inspect direct-feed fit before adding input IDC')
    if SIMULATE_FULL_RESPONSE and FIX_T2_JUNCTION:
        from geometry.taper.taper3_fixed_junction import run
        run(CASE, PARAMS, measured['L0_nH'])
        return
    if SIMULATE_FULL_RESPONSE:
        result = fit_response(CASE, PARAMS, measured['L0_nH'])
        # These use saved S-parameters only; they do not launch more EM jobs.
        from geometry.taper.finalize_taper3 import plot_case
        from geometry.taper.taper3_fit_sensitivity import analyze
        plot_case(CASE)
        analyze(CASE)
        if not result['fit_pass']:
            fit_older_within_one_percent(CASE, measured['L0_nH'], stage='response')


if __name__ == '__main__':
    main()
