"""Editable T2-only iteration. Verified T1 and right-hand geometry are locked.

Run from the project root: .venv/Scripts/python.exe -B -m geometry.taper.run_locked_taper2
No command-line arguments. Start with RUN_EM=False and inspect geometry_review.png.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from geometry.meander import MeanderParams
from geometry.taper.layout import TaperParams
from geometry.taper.second_layout import SecondParams
from geometry.taper.locked_taper2 import T1,prepare_case,run

CASE_NAME='compact_idc30_span172_top148_w156_ground80_rise358_right572'
RUN_EM=False  # Latest case is rejected; inspect it before starting a new case.
GEOMETRY=SecondParams(
    taper=TaperParams(w=156,t_top=14,t_left=24,t_right=20,l_left=584,l_right=572,
                     meander=MeanderParams(w=16,l=596,l_top=10,n=3,t=4,d=6)),
    taper1=T1,idc_pairs=30,idc_span=172,idc_t=4,idc_d=4,d_idc_top=10,
    d_t2_t1=6,d_t1_regular=56,t_junction_ground=4,
    l_ground_rise=358,w_ground_riser=80,straight_inductor=True,x_meander=0,
    w_top_left_extension=148,
)


def main():
    if RUN_EM:
        print(run(CASE_NAME,GEOMETRY),flush=True)
    else:
        gds,out=prepare_case(CASE_NAME,GEOMETRY)
        print('GDS:',gds)
        print('Visual review:',out/'geometry_review.png')
        print('Verified T1/right-hand geometry XOR: empty. No EM launched.')


if __name__=='__main__':main()
