"""Controlled bare-T2 tests with exact archived parameter reconstruction."""
import json
import gdstk
from geometry.taper.progressive_fixtures import ROOT,GDS,OUT,prepare as prepare_output,run as run_output
from geometry.taper.progressive_trial import measure
from geometry.taper.compact_stages import build
from geometry.taper.progressive_fitting import fit
from geometry.taper.staged_fitting import LT1,LR


def source(name,p,stage):
    folder=GDS/name/'source'/stage;folder.mkdir(parents=True,exist_ok=True)
    path=folder/'geometry.gds';polys,meta,_=build(p,stage)
    if path.exists():
        if json.loads(path.with_suffix('.json').read_text())!=meta:raise ValueError('Changed compact control')
    else:
        lib=gdstk.Library(unit=1e-6,precision=1e-9);lib.new_cell('TOP').add(*polys);lib.write_gds(str(path))
        path.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n')
    return path


def prepare(name,p):
    return prepare_output(name,str(source(name,p,'direct')),3,302)


def run(name,p):
    prepare(name,p);cal=measure(name,p,calibration_stage='Lcal');l=cal['L0_nH']
    if abs(l/.74-1)>.15:raise ValueError(f'Bare control L={l}; adjust wire before full simulation')
    snp,out=run_output(name,str(source(name,p,'direct')),3,302)
    b=json.loads((OUT/'01_T1_regular/twoC/fit.json').read_text())
    args=([l,LT1,LR],[0,b['C_fF'][0],190.4],[None,56.7,8.2865])
    rows={}
    for key,free in (('zero_C',[('G',1)]),('bare_C_free',[('C',0),('G',1)])):
        rows[key]=fit(snp,out/key,*args,free,title='Bare T2 Control: '+key.replace('_',' '))
    (out/'control_summary.json').write_text(json.dumps(rows,indent=2)+'\n')
    return rows
