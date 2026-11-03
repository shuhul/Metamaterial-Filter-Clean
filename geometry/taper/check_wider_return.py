"""Keep a simple straight wire and test a wider grounded return."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from dataclasses import replace
from geometry.taper.check_wider_trace import parameters as previous_parameters
from geometry.taper.regular_cap_tab_case import write_tabbed,RegularCapTab
from geometry.taper.measure_response_inductance import measure

NAME='response_w16_return48_rise220_left12_xm148_span290_tab116'

def parameters():
    p=previous_parameters()
    return replace(p,w_ground_riser=48,x_meander=-148,
                   taper=replace(p.taper,t_left=12))

if __name__=='__main__':
    write_tabbed(NAME,parameters(),RegularCapTab(w=116,t=4))
    measure(NAME)
