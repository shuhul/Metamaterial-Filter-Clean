"""Native validation of local coupling guidance and one pad-length control."""
from dataclasses import replace
from geometry.taper.taper2_open_meander_trials import input_settings as measured
from geometry.taper.parallel_trials import run

base=replace(measured,coupling=replace(measured.coupling,d=136,d_t1=38,d_regular=70))
CASES=[
 ('idc_n24fg20_right6_dl256_dr42_pad148_w152_p594_i136_j38_r70tg2',base,
  'Matched-body extracted-C model with native inverse-response correction '
  'proposes input136/junction38/regular70um, predictedRMS0.03414. '
  'Prediction is not acceptance. This tests the proposed coupling balance '
  'against the full unchanged nominal target and native rawS21.'),
 ('idc_n24fg20_right6_dl256_dr42_pad148_w152_p200_i136_j38_r70tg2',
  replace(base,pad=replace(base.pad,l=200)),
  'Only existing T2 ground-pad length594 to200um differs from its companion. '
  'Source fitted T2 C0 is158.7fF versus149fF target; shorten pads to reduce '
  'shunt capacitance. New body requires independent explicit-CAP L before S21.')]

if __name__=='__main__':
    run(CASES,'sparse_local136_j38_r70_pad594_vs200')
