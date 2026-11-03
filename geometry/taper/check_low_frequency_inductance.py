"""Check T2 inductance around the design band using saved EM geometry/cache.

Four larger ideal C_add loads move the calibration peaks near 5 GHz. The
cache-only engine flag forbids new direct EM frequencies. Original trials
remain untouched; the entire small isolation cache is copied to a new case.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from pathlib import Path
import json
import shutil
from geometry.taper.locked_taper2 import ROOT
from geometry.taper.second_layout import validate
from geometry.taper.sonnet import solver_environment
from extract_generated_inductance import main as extract_l

CASE_NAME='straight4_w400_xm76_idc30_span180_ground160_rise338_right536'


def check(name):
    source_name=name+'_isolation';case=name+'_lf_isolation'
    source=ROOT/'sonnet_files/taper_exploration'/source_name
    dest=ROOT/'sonnet_files/taper_exploration'/case
    out=ROOT/'output_data/taper_exploration'/case
    if (out/'results.json').exists():return json.loads((out/'results.json').read_text())
    if not dest.exists():shutil.copytree(source,dest)
    out.mkdir(parents=True,exist_ok=True)
    for ext in ('.gds','.json'):
        old=ROOT/'geometry/taper/output'/f'{source_name}{ext}'
        new=ROOT/'geometry/taper/output'/f'{case}{ext}'
        if not new.exists():shutil.copyfile(old,new)
    gds=ROOT/'geometry/taper/output'/f'{case}.gds'
    validate(gds,json.loads(gds.with_suffix('.json').read_text()))
    result=extract_l(project_file=dest/'generated_ncell_1.sonx',output_dir=out,
        fit_capacitance=False,c_add_sweep_pf=(.8,1.,1.2,1.4),target_l0_nh=.74,
        solver_env=solver_environment(),cache_only=True)
    result.update(purpose='Design-band inductance consistency check from existing EM cache',
                  source_case=name,direct_EM_frequencies_added=0)
    (out/'results.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result


if __name__=='__main__':check(CASE_NAME)
