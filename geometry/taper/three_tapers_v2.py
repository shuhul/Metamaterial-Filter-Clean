"""Current T1 fixture: feed electrode -> T1 -> regular -> output end.

There is no input-end resonator. Only three electrical capacitances are fit.
"""
from dataclasses import asdict,replace
from types import SimpleNamespace
from functools import partial
import json
import gdstk
import numpy as np
from geometry.unit_cell import make_unit_cell
from geometry.specifications import SPECIFICATION_VERSION,specification
from geometry.taper import three_tapers as shared

ROOT=shared.ROOT
CAMPAIGN='three_tapers_v2_2026-10-25'  # Historical geometry, with input ground strip.
BASE=replace(shared.BASE,meander=replace(shared.BASE.meander,w=220),
             pad=replace(shared.BASE.pad,l=420))
CASE='finger6_w220_pad8_l420_g26_j80'
D_INPUT=26
D_JUNCTION=80
D_INPUT_GROUND=12  # Clearance from the input strip to the T1 cap.
D_REGULAR_GROUND=26  # Clearance from the junction strip to the regular cap.
RUN_INDUCTANCE=True
RUN_RESPONSE=True
FINALIZE_ON_PASS=True


def build(p,*,isolation=False,d_input=40,d_junction=56,d_regular_ground=None,d_input_ground=None,d_isolation_output=None):
    if p.d_top+p.cap.t_top+p.meander.l!=618:
        raise ValueError('T1 must terminate on the flat bottom ground')
    if any(v<8 or v%2 for v in (d_input,d_junction)):
        raise ValueError('Grounded gaps must be even and >=8 um')
    if (d_input%4 and d_input_ground is None) or (d_junction%4 and d_regular_ground is None):
        raise ValueError('A 2-mod-4 gap needs an explicit on-grid strip clearance')
    units=[p,p,p] if isolation else [p,p,shared.REGULAR,shared.REGULAR]
    gaps=[d_input,d_input if d_isolation_output is None else d_isolation_output] if isolation else [d_input,d_junction,56]
    lefts=[426]
    for u,gap in zip(units,gaps):lefts.append(lefts[-1]+u.cap.w+gap)
    centers=[x+u.cap.w/2 for x,u in zip(lefts,units)]
    width=lefts[-1]+units[-1].cap.w+426
    top=gdstk.Cell('TOP');parts=[]
    # The left scaffold has ONE facing electrode connected to the feed.
    # Its meander, pads, other leg, and former cap top are all absent.
    edge=lefts[0]+p.cap.w
    parts.extend([gdstk.rectangle((edge-p.cap.t_side,-8-p.cap.l),(edge,-8)),
                  gdstk.rectangle((0,-18),(edge,-8))])
    top.add(gdstk.Reference(make_unit_cell(p,name='T1'),origin=(centers[1],0)))
    if isolation:
        # A second weak electrode terminates the isolated L fixture.
        edge=lefts[2]
        parts.extend([gdstk.rectangle((edge,-8-p.cap.l),(edge+p.cap.t_side,-8)),
                      gdstk.rectangle((edge,-18),(width,-8))])
    else:
        for i in (2,3):
            top.add(gdstk.Reference(make_unit_cell(shared.REGULAR,name=f'REGULAR_{i}'),origin=(centers[i],0)))
        parts.append(gdstk.rectangle((centers[-1],-18),(width,-8)))
    parts.extend([gdstk.rectangle((0,0),(width,448)),
                  gdstk.rectangle((0,-1052),(width,-618)),
                  gdstk.rectangle((0,-1052),(412,-24)),
                  gdstk.rectangle((width-412,-1052),(width,-24))])
    for i,gap in enumerate(gaps):
        x=lefts[i]+units[i].cap.w+gap/2
        if d_input_ground is not None and (i==0 or (isolation and i==1)):
            if d_input_ground<2 or d_input_ground%2 or gap-d_input_ground-4<2:
                raise ValueError('Input ground must clear both electrodes on grid')
            x=lefts[1]-d_input_ground-2 if i==0 else lefts[1]+p.cap.w+d_input_ground+2
        if not isolation and i==1 and d_regular_ground is not None:
            if d_regular_ground<2 or d_regular_ground%2 or gap-d_regular_ground-4<2:
                raise ValueError('Ground strip must clear both cap legs on the 2 um grid')
            x=lefts[i+1]-d_regular_ground-2
        parts.append(gdstk.rectangle((x-2,-618),(x+2,0)))
    top.add(*gdstk.boolean(parts,[],'or',precision=.001,layer=1))
    meta={'layout_kind':'T1_single_input_electrode_v2','unit':asdict(p),'isolation':isolation,
          'regular':asdict(shared.REGULAR),'d_input':d_input,'d_junction':d_junction,
          'w_sim':width,'cell_centers':centers,'input_end_resonator':False,
          'expected_metal_components':3 if isolation else 2,
          'parameters':{'n_cells':1 if isolation else 2,'l_sim':1500,'l_top_ground':448,
                        'w_side':412,'t_feed':10,'unit':{'d_top':8,'cap':{'t_top':14}}}}
    if d_regular_ground is not None:meta['d_regular_ground']=d_regular_ground
    if d_input_ground is not None:meta['d_input_ground']=d_input_ground
    if d_isolation_output is not None:meta['d_isolation_output']=d_isolation_output
    return top,meta


def validate(path,meta):
    top,rebuilt=build(shared.from_dict(meta['unit']),isolation=meta['isolation'],
                      d_input=meta['d_input'],d_junction=meta['d_junction'],
                      d_regular_ground=meta.get('d_regular_ground'),d_input_ground=meta.get('d_input_ground'),
                      d_isolation_output=meta.get('d_isolation_output'))
    if rebuilt!=meta:raise ValueError('Geometry metadata differs')
    lib=gdstk.read_gds(str(path));tops=lib.top_level()
    if len(tops)!=1 or tops[0].name!='TOP' or not np.isclose(lib.unit,1e-6):raise ValueError('Invalid GDS')
    polys=tops[0].get_polygons()
    if any(q.layer!=1 or q.datatype!=0 for q in polys):raise ValueError('Invalid layer')
    if gdstk.boolean(polys,top.get_polygons(),'xor',precision=.001):raise ValueError('GDS differs from parameters')
    if len(gdstk.boolean(polys,[],'or',precision=.001))!=(3 if meta['isolation'] else 2):
        raise ValueError('Input electrode or ground has wrong connectivity')
    return SimpleNamespace(n_cells=meta['parameters']['n_cells'],l_sim=1500,l_top_ground=448),polys


def prepare(name,p,*,d_regular_ground=None,d_input_ground=None,**kwargs):
    return shared.prepare(name,p,campaign=CAMPAIGN,
                          builder=partial(build,d_regular_ground=d_regular_ground,d_input_ground=d_input_ground),validator=validate,**kwargs)


def measure(name,p,*,d_input=40,d_input_ground=None):
    return shared.measure(name,p,d_input=d_input,prepare_case=partial(prepare,d_input_ground=d_input_ground))


def response(name,p,l,*,d_input=40,d_junction=56,d_regular_ground=None,d_input_ground=None):
    return shared.response(name,p,l,d_input=d_input,d_junction=d_junction,
                           prepare_case=partial(prepare,d_regular_ground=d_regular_ground,d_input_ground=d_input_ground),input_end_cell=False)


def finalize(name):
    return shared.finalize(name,campaign=CAMPAIGN,geometry_validator=validate)


def main():
    selected=ROOT/'geometry/taper/selected_three_tapers.json'
    if selected.exists():
        saved=json.loads(selected.read_text())
        if saved['case']==CASE:
            _,meta=build(BASE,d_input=D_INPUT,d_junction=D_JUNCTION,
                         d_regular_ground=D_REGULAR_GROUND,d_input_ground=D_INPUT_GROUND)
            if saved['geometry']!=meta or saved['specification']!=specification():
                raise ValueError('Accepted case is immutable; choose a new CASE for changed geometry or targets')
            print('Already verified:',saved['project']);print('Plot:',saved['plot'])
            return
    if RUN_INDUCTANCE:
        measured=measure(CASE,BASE,d_input=D_INPUT,d_input_ground=D_INPUT_GROUND)
    else:
        base=ROOT/'geometry/taper/output'/CAMPAIGN/CASE/'Lcal/geometry.json'
        _,expected=build(BASE,isolation=True,d_input=D_INPUT,d_input_ground=D_INPUT_GROUND)
        if json.loads(base.read_text())!=expected:raise ValueError('Saved L has different geometry')
        measured=json.loads((ROOT/'output_data/taper_exploration'/CAMPAIGN/CASE/'Lcal/results.json').read_text())
    if abs(measured['target_error_percent'])>2:
        raise RuntimeError('Adjust isolated L before running the response')
    if RUN_RESPONSE:
        score,_=response(CASE,BASE,measured['L0_nH'],d_input=D_INPUT,d_junction=D_JUNCTION,
                         d_regular_ground=D_REGULAR_GROUND,d_input_ground=D_INPUT_GROUND)
        if score['pass'] and FINALIZE_ON_PASS:finalize(CASE)


if __name__=='__main__':main()
