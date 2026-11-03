"""Authorized input-only IDC on the approved asymmetric T2 body.

Edit grouped settings here; the independent side clearances are retained.
One four-thread EM job at a time. Probe coupling before the full response.
"""
from dataclasses import replace,asdict
import json
from geometry.taper import taper2_plain as workflow

SETTINGS=replace(workflow.SETTINGS,allow_input_idc=True,
    cap=replace(workflow.SETTINGS.cap,t_left=8,t_right=30),
    clearance=replace(workflow.SETTINGS.clearance,d_left=256),
    meander=replace(workflow.SETTINGS.meander,w=148),
    pad=replace(workflow.SETTINGS.pad,l=200),
    d_pad_left=148,
    coupling=replace(workflow.SETTINGS.coupling,d=42,n=70,d_finger=4,t_input=8,d_t1=42,
                     d_regular=74,t_regular_ground=2))
CASE='idc_rail8_right30_dl256_pad148_w148_p200_i42n70_j42_r74tg2'
# Explicit capacitor sweeps are reused automatically for identical T2 bodies.
RUN_INPUT_PROBE=False
RUN_RESPONSE=True


def main():
    selected=workflow.ROOT/'geometry/taper/selected_taper2_idc.json'
    if selected.exists():
        saved=json.loads(selected.read_text())
        if (saved['case']==CASE and saved['settings']==asdict(SETTINGS) and
                saved['specification']==workflow.common.specification()):
            print('Already verified:',saved['project']);print('Plot:',saved['plot']);return
    project,out,meta=workflow.prepare(CASE,SETTINGS)
    print('PROJECT',project,flush=True)
    print('IDC CLEARANCES',meta['IDC_clearance_um'],flush=True)
    if RUN_INPUT_PROBE:
        probe=workflow.probe_input(CASE,SETTINGS,frequencies=(4.65,5.025,5.4))
        if abs(probe['pointwise_input_Cg_fF'][1]/workflow.common.TAPER_2.Cg_input_fF-1)>.06:
            raise ValueError('Tune the input coupling before spending a full response sweep')
    if RUN_RESPONSE:
        measured=workflow.measure(CASE,SETTINGS)
        if abs(measured['target_error_percent'])>2:
            raise ValueError('Measured inductance is outside +/-2%; skip the full response sweep')
        score,fitted=workflow.response(CASE,SETTINGS,measured['L0_nH'])
        if score['pass'] and fitted['curve_pass'] and abs(measured['target_error_percent'])<=2:
            workflow.finalize(CASE)
        print(json.dumps(score,indent=2),flush=True)


if __name__=='__main__':main()
