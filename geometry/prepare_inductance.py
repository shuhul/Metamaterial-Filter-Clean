"""Add the baseline's ideal parallel capacitor to the generated single cell."""

from pathlib import Path
from copy import deepcopy
import json
import math
import xml.etree.ElementTree as ET
import gdstk
import numpy as np

if __package__:
    from .prepare_sonnet import BASE, BIN, tag, write_xml, run, require, points
    from .validation import validate_gds_metadata
else:
    from prepare_sonnet import BASE, BIN, tag, write_xml, run, require, points
    from validation import validate_gds_metadata


def main(*, gds_path=None, project_path=None, output_path=None, sweep=(4, 6, 20000)):
    gds = Path(gds_path) if gds_path is not None else BASE/'geometry'/'output'/'ncell_1.gds'
    project = Path(project_path) if project_path is not None else BASE/'sonnet_files'/'generated_ncell_1'/'generated_ncell_1.sonx'
    folder = project.parent
    metadata = json.loads(gds.with_suffix('.json').read_text())
    params, source_polys = validate_gds_metadata(gds, metadata)
    require(params.n_cells == 1, 'Ideal-capacitor extraction requires one interior cell')
    x = 2 * math.floor(metadata['cell_centers'][1] / 2 + .5)
    ground_y = metadata['parameters']['l_top_ground']
    cap_y = ground_y + metadata['parameters']['unit']['d_top']
    tree = ET.parse(project)
    root = tree.getroot()
    original = ET.parse(BASE / "sonnet_files" / "filter_v6.sonx")
    component = deepcopy(next(c for c in original.findall(f".//{tag('Component')}")
                              if c.findtext(tag('Value')) == 'C_add'))
    geometry = root.find(tag('Geometry'))
    expected = []
    for p in source_polys:
        xy = p.points.copy()
        xy[:,1] += params.l_sim-params.l_top_ground
        xy = np.floor(xy/2+.5)*2
        xy[:,1] = params.l_sim-xy[:,1]
        expected.append(gdstk.Polygon(xy))
    actual = [gdstk.Polygon(points(p)) for p in geometry.findall(f'.//{tag("PlanarPolygon")}')]
    require(not gdstk.boolean(expected,actual,'xor',precision=.001),
            'Single-cell Sonnet geometry is stale; rerun prepare_sonnet first')
    require(gdstk.inside([(x,ground_y-1),(x,ground_y+1),(x,cap_y-1),(x,cap_y+1)],actual)
            == (True,False,False,True), 'Component pins do not bridge the intended metal edges')
    require((not geometry.findall(f'.//{tag("Component")}')), 'Validation failed: not geometry.findall(f\'.//{tag("Component")}\')')
    component.set('MacroID', 'componentGeneratedCadd')
    ids = [int(p.get('Id')) for p in geometry.findall(f'.//{tag("PlanarPolygon")}')]
    component.set('Id', str(max(ids)+1))
    # Match the grid-snapped central lead and bridge the ground-to-cap gap.
    component.find(tag('LabelPos')).text = f'({x+6},{(ground_y+cap_y)/2})'
    component.find(tag('SymbolBox')).text = f'({x-2},{ground_y+2},{x+2},{cap_y-2})'
    for i, pin in enumerate(component.findall(tag('Pin'))):
        pin.set('MacroID', f'componentPinGenerated{i+1}')
        pin.find(tag('Center')).text = f'({x},{ground_y if i == 0 else cap_y})'
    level = next(l for l in geometry.findall(tag('Level')) if l.get('Number') == '0')
    level.find(tag('Objects')).append(component)
    variable = next(v for v in geometry.findall(tag('Variable')) if v.get('Name') == 'C_add')
    variable.set('Units', 'CAPACITANCE')
    variable.set('Value', '0.0001')
    variable.set('Description', 'Ideal shunt capacitor for inductance extraction (pF)')
    for sweep_set in root.findall(f'{tag("Sweeps")}/{tag("Set")}'):
        if sweep_set.get('On') == 'TRUE':
            frequency, = sweep_set.findall(f'{tag("Frequencies")}/{tag("Sweep")}')
            frequency.set('Start',str(sweep[0]))
            frequency.set('Stop',str(sweep[1]))
            frequency.set('Target',str(sweep[2]))
    target = Path(output_path) if output_path is not None else folder/'generated_inductance.sonx'
    for export in root.findall(f'{tag("OutputFiles")}//{tag("File")}'):
        if export.get('Type') == 'TOUCH':
            export.set('Name',target.stem+'.snp')
    write_xml(tree,target)
    run([BIN/'soncmd.exe','-Read',target], folder/'inductance_read.log')
    print(f'Prepared: {target}')


if __name__ == '__main__':
    main()
