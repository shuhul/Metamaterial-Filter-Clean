"""Two independent small coupling controls on the actual best native body."""
from dataclasses import replace
from geometry.taper.taper2_idc import SETTINGS as baseline
from geometry.taper.parallel_trials import run

CASES=[
 ('idc_rail8_right30_dl256_pad148_w148_p200_i42n70_j50_r74tg2',
  replace(baseline,coupling=replace(baseline.coupling,d_t1=50)),
  'Only T2/T1 gap42 to50um from the actual0.03947RMS best. Establish '
  'an uncontaminated local junction derivative, with body/input/unit '
  'junction unchanged and explicitL4.323385nH reused.'),
 ('idc_rail8_right30_dl256_pad148_w148_p200_i42n70_j42_r78tg2',
  replace(baseline,coupling=replace(baseline.coupling,d_regular=78)),
  'Only centered T1/unit gap74 to78um from the actual best. Together '
  'with saved pure input, pad-length, cap-leg and ground-width controls, '
  'this completes six independent local geometry directions.')]

if __name__=='__main__':
    run(CASES,'best_body_pure_j50_and_r78')
