"""Active T1: ungrounded taper input gap, centered regular-junction strip.

Electrical targets are unchanged; fixed-target complex S21 RMS must be <= .025.
Run as a module; edit the parameters and flags below for a fresh named trial.
"""
from dataclasses import replace
from functools import partial
from types import SimpleNamespace
import json
import gdstk
import numpy as np
from geometry.specifications import SPECIFICATION_VERSION,specification
from geometry.taper import three_tapers as shared
from geometry.taper import three_tapers_v2 as previous

ROOT=shared.ROOT
CAMPAIGN=SPECIFICATION_VERSION
BASE=replace(previous.BASE,cap=replace(previous.BASE.cap,w=344,l=590),
             pad=replace(previous.BASE.pad,l=532))
CASE='plain_cw344_w220_p532_g56_j66_c590_tg6'
D_INPUT=56
D_JUNCTION=66
T_JUNCTION_GROUND=6
RUN_INDUCTANCE=True
RUN_RESPONSE=True
FINALIZE_ON_PASS=True


def build(p,*,isolation=False,d_input=96,d_junction=64,t_junction_ground=4):
    if d_input<8 or d_input%2 or d_junction<8 or d_junction%2:
        raise ValueError('Cap gaps must be even and at least 8 um')
    weak_gap=96+t_junction_ground%4 if isolation else d_junction
    if (t_junction_ground<4 or t_junction_ground%2 or t_junction_ground>weak_gap-4
        or (weak_gap-t_junction_ground)%4):
        raise ValueError('Centered ground edges need (gap - strip width) divisible by 4')
    temporary=2*np.floor((d_input-4)/4) if d_input%4 else None
    regular_temporary=2*np.floor((d_junction-4)/4) if d_junction%4 else None
    original,meta=previous.build(p,isolation=isolation,d_input=d_input,d_junction=d_junction,
                                 d_input_ground=temporary,d_regular_ground=regular_temporary,
                                 d_isolation_output=weak_gap if isolation else None)
    # Remove precisely the former input ground strip, never a cap or feed.
    # Remaining strips use the original builder's exact gap centers.
    edge=426+p.cap.w
    cut=gdstk.rectangle((edge,-618),(edge+d_input,0))
    # Keep component references separate: merging the entire meander into
    # ground before GDS export needlessly fractures it into many polygons.
    top=gdstk.Cell('TOP')
    top.add(*[r.copy() for r in original.references])
    cuts=[cut]
    if isolation and temporary is not None:
        # Keep the measurement's output electrode and centered strip fixed
        # while changing only its input gap, even for 2-mod-4 input gaps.
        t1_right=meta['cell_centers'][1]+p.cap.w/2
        old_x=t1_right+temporary+2
        cuts.append(gdstk.rectangle((old_x-2,-618),(old_x+2,0)))
        mid=t1_right+weak_gap/2
        top.add(gdstk.rectangle((mid-2,-618),(mid+2,0),layer=1))
    for poly in original.polygons:
        if gdstk.boolean([poly],cuts,'and',precision=.001):
            top.add(*gdstk.boolean([poly],cuts,'not',precision=.001,layer=1))
        else:top.add(poly.copy())
    if t_junction_ground!=4:
        # Widen only the existing centered T1-output strip, never the input gap.
        t1_right=meta['cell_centers'][1]+p.cap.w/2
        mid=t1_right+weak_gap/2
        strip=gdstk.rectangle((mid-t_junction_ground/2,-618),
                              (mid+t_junction_ground/2,0),layer=1)
        merged=gdstk.boolean([*top.polygons,strip],[],'or',precision=.001,layer=1)
        top.remove(*top.polygons);top.add(*merged)
        meta['regular_ground_strip_width_um']=t_junction_ground
    meta['layout_kind']='T1_no_input_ground_centered_regular_v3'
    meta['input_ground_strip']=False
    meta['regular_ground_strip']='centered'
    meta.pop('d_input_ground',None)
    meta.pop('d_regular_ground',None)
    meta.pop('d_isolation_output',None)
    if isolation and d_input!=weak_gap:meta['isolation_output_gap_um']=weak_gap
    return top,meta


def validate(path,meta):
    top,expected=build(shared.from_dict(meta['unit']),isolation=meta['isolation'],
                        d_input=meta['d_input'],d_junction=meta['d_junction'],
                        t_junction_ground=meta.get('regular_ground_strip_width_um',4))
    if expected!=meta:raise ValueError('Stale V3 geometry metadata')
    lib=gdstk.read_gds(str(path));tops=lib.top_level()
    if len(tops)!=1 or tops[0].name!='TOP' or not np.isclose(lib.unit,1e-6):raise ValueError('Invalid GDS')
    polys=tops[0].get_polygons()
    if any(q.layer!=1 or q.datatype!=0 for q in polys):raise ValueError('Invalid layer/datatype')
    if gdstk.boolean(polys,top.get_polygons(),'xor',precision=.001):raise ValueError('GDS differs from parameters')
    if len(gdstk.boolean(polys,[],'or',precision=.001))!=meta['expected_metal_components']:
        raise ValueError('Wrong feed/ground connectivity')
    input_edge=426+meta['unit']['cap']['w']
    opening=gdstk.rectangle((input_edge,-618),(input_edge+meta['d_input'],0))
    if gdstk.boolean(polys,[opening],'and',precision=.001):
        raise ValueError('The taper input gap must contain no ground strip')
    width=meta.get('regular_ground_strip_width_um',4)
    gap=96+width%4 if meta['isolation'] else meta['d_junction']
    edge=meta['cell_centers'][1]+meta['unit']['cap']['w']/2
    mid=edge+gap/2
    junction=gdstk.rectangle((edge,-618),(edge+gap,0))
    actual=gdstk.boolean(polys,[junction],'and',precision=.001)
    strip=gdstk.rectangle((mid-width/2,-618),(mid+width/2,0))
    if gdstk.boolean(actual,[strip],'xor',precision=.001):
        raise ValueError('The junction strip must be exactly centered between cap edges')
    upper=gdstk.rectangle((0,0),(meta['w_sim'],448))
    if gdstk.boolean([upper],polys,'not',precision=.001):
        raise ValueError('The upper ground plane must remain flat and uncut')
    return SimpleNamespace(n_cells=meta['parameters']['n_cells'],l_sim=1500,l_top_ground=448),polys


def prepare(name,p,*,t_junction_ground=4,**kwargs):
    return shared.prepare(name,p,campaign=CAMPAIGN,
                          builder=partial(build,t_junction_ground=t_junction_ground),validator=validate,**kwargs)


def measure(name,p,*,d_input=96,t_junction_ground=4):
    return shared.measure(name,p,d_input=d_input,
                          prepare_case=partial(prepare,t_junction_ground=t_junction_ground))


def response(name,p,l,*,d_input=96,d_junction=64,t_junction_ground=4):
    return shared.response(name,p,l,d_input=d_input,d_junction=d_junction,
                           prepare_case=partial(prepare,t_junction_ground=t_junction_ground),input_end_cell=False)


def finalize(name):
    return shared.finalize(name,campaign=CAMPAIGN,geometry_validator=validate)


def run_case(name,p,*,d_input=96,d_junction=64,t_junction_ground=4):
    measured=measure(name,p,d_input=d_input,t_junction_ground=t_junction_ground)
    if abs(measured['target_error_percent'])>2:
        raise ValueError('Isolated L must pass before the response simulation')
    return response(name,p,measured['L0_nH'],d_input=d_input,d_junction=d_junction,
                    t_junction_ground=t_junction_ground)


def main():
    selected=ROOT/'geometry/taper/selected_three_tapers.json'
    if selected.exists():
        saved=json.loads(selected.read_text())
        if saved['case']==CASE and saved['specification']==specification():
            _,expected=build(BASE,d_input=D_INPUT,d_junction=D_JUNCTION,t_junction_ground=T_JUNCTION_GROUND)
            if expected!=saved['geometry']:raise ValueError('Use a new CASE for changed geometry')
            print('Already verified:',saved['project']);print('Plot:',saved['plot']);return
    if RUN_INDUCTANCE:measured=measure(CASE,BASE,d_input=D_INPUT,t_junction_ground=T_JUNCTION_GROUND)
    else:
        path=ROOT/'geometry/taper/output'/CAMPAIGN/CASE/'Lcal/geometry.json'
        _,expected=build(BASE,isolation=True,d_input=D_INPUT,t_junction_ground=T_JUNCTION_GROUND)
        if json.loads(path.read_text())!=expected:raise ValueError('Saved calibration geometry differs')
        measured=json.loads((ROOT/'output_data/taper_exploration'/CAMPAIGN/CASE/'Lcal/results.json').read_text())
    if abs(measured['target_error_percent'])>2:raise ValueError('Isolated L is outside tolerance')
    if RUN_RESPONSE:
        score,_=response(CASE,BASE,measured['L0_nH'],d_input=D_INPUT,d_junction=D_JUNCTION,
                         t_junction_ground=T_JUNCTION_GROUND)
        if score['pass'] and FINALIZE_ON_PASS:finalize(CASE)


if __name__=='__main__':main()
