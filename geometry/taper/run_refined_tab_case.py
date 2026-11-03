"""Fresh EM check of the measured-secant geometry correction."""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
from dataclasses import replace
import json
from geometry.taper.regular_cap_tab_case import run, RegularCapTab, ROOT
from geometry.taper.second_layout import from_dict

NAME = 'response_match_span290_pairs31_rise32_right548_tab116'


def main():
    source = ROOT / 'geometry/taper/output/response_match_span262_rise20_right558_tab116.json'
    p = from_dict(json.loads(source.read_text())['second'])
    p = replace(p, idc_span=290, idc_pairs=31, l_ground_rise=32,
                taper=replace(p.taper, l_right=548))
    run(NAME, p, RegularCapTab(w=116, t=4))


if __name__ == '__main__':
    main()
