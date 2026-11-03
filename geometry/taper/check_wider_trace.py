"""Check a wider straight T2 trace without repeating a full response sweep."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from dataclasses import replace
import json
from geometry.taper.regular_cap_tab_case import ROOT,write_tabbed,RegularCapTab
from geometry.taper.second_layout import from_dict
from geometry.taper.measure_response_inductance import measure

NAME='response_w16_rise220_span290_pairs31_right548_tab116'

def parameters():
    path=ROOT/'geometry/taper/output/response_match_span290_pairs31_rise32_right548_tab116.json'
    p=from_dict(json.loads(path.read_text())['second'])
    return replace(p,l_ground_rise=220,
                   taper=replace(p.taper,meander=replace(p.taper.meander,w=52,t=16)))

if __name__=='__main__':
    write_tabbed(NAME,parameters(),RegularCapTab(w=116,t=4))
    measure(NAME)
