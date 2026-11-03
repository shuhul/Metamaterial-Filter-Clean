"""Verify that new response cases retain the selected T1 EM settings."""
import json
import xml.etree.ElementTree as ET
from geometry.prepare_sonnet import tag
from geometry.taper.locked_taper2 import VERIFIED
from geometry.taper.optimize_response import ROOT,OUT

def canon(q):
    return [q.tag,sorted(q.attrib.items()),(q.text or '').strip(),[canon(x) for x in q]]

def signature(path):
    r=ET.parse(path).getroot();g=r.find(tag('Geometry'));box=g.find(tag('Box'))
    sz=box.find(tag('Size'));cells=box.find(tag('NumCells'))
    return {'units':canon(r.find(tag('Units'))),'control':canon(r.find(tag('Control'))),
        'conductors':[canon(q) for q in g.findall(tag('Conductor'))],
        'dielectrics':[canon(q) for q in g.findall(tag('Dielectric'))],
        'stack':[[q.attrib,[canon(x) for x in q if x.tag!=tag('Objects')]] for q in g.findall(tag('Level'))],
        'planar_material_and_mesh':[[q.attrib,[canon(x) for x in q if x.tag!=tag('Objects')]] for q in r.iter(tag('PlanarTechLayer'))],
        'covers':[canon(q) for q in box.findall(tag('BoxCover'))],
        'reference_planes':[canon(q) for q in box.findall(tag('RefPlane'))],
        'height_um':sz.get('Y'),'stack_height_um':sz.get('Z'),
        'grid_um':[float(sz.get('X'))/float(cells.get('X')),float(sz.get('Y'))/float(cells.get('Y'))],
        'external_port_impedances':sorted([[q.get('Number'),q.find(tag('Type')).attrib,
            q.find(tag('Impedance')).attrib] for q in r.iter(tag('Port')) if q.get('Number') in ('1','2')])}

def main():
    reference=VERIFIED['artifacts']['sonnet_project'];old=signature(reference)
    cases=['response_w400_xm160_span206_rise0_right578','response_vertical16_l480_w328_rise358_v2','response_best_regular_tab100_t4']
    rows=[]
    for name in cases:
        p=ROOT/'sonnet_files/taper_exploration'/name/'generated_ncell_3.sonx';new=signature(p)
        checks={k:old[k]==new[k] for k in old}
        rows.append({'case':name,'checks':checks,'all_equal':all(checks.values())})
    result={'reference':reference,'checks':rows,
        'intended_difference':'Box width follows T2 geometry; both cells remain 2 um, height 1500 um. Metadata/import dates and polygon IDs differ.'}
    (OUT/'EM_settings_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    if not all(r['all_equal'] for r in rows):raise ValueError(result)
    print(json.dumps(result,indent=2))

if __name__=='__main__':main()
