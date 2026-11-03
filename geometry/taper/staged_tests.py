r"""Prepare Frank's staged extraction fixtures. Never launch an EM solve.

Run from the project root:
    .venv\Scripts\python.exe -B -m geometry.taper.staged_tests

All lengths are microns. Inspection candidates are not accepted designs.
Only Sonnet's GDS importer and native project-read check are invoked.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
from copy import deepcopy
from dataclasses import asdict
from hashlib import sha256
from pathlib import Path
from types import SimpleNamespace
import json
import xml.etree.ElementTree as ET
import gdstk
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from geometry.prepare_sonnet import BIN, SOURCE, tag, points, write_xml, run
from geometry.taper.locked_taper2 import T1
from geometry.meander import make_meander
from geometry.placement import place
from geometry.taper.layout import make_taper
from geometry.specifications import specification

ROOT = Path(__file__).resolve().parents[2]
NAME = 'frank_staged_tests'
NATIVE = ROOT / 'sonnet_files/taper_exploration' / NAME
GDS_ROOT = ROOT / 'geometry/taper/output' / NAME

# EDITABLE INSPECTION SETTINGS. No width is claimed to give 0.74 nH yet.
WIRE_WIDTHS = (16, 32, 64)
INTEGRATION_WIRE = 32
BOX_W, BOX_H = 2592, 1500
TOP_GROUND, BOTTOM_GROUND = 0, -618
TOP_GROUND_HEIGHT = 448
T2_X, T1_X = 916, 1286
T2 = dict(w=400, t_top=14, t_left=24, t_right=20,
          l_left=584, l_right=548)
IDC = dict(pairs=31, t=4, d=4, span=290, d_top=10, trunk_t=20)


def rect(x0, y0, x1, y1):
    return gdstk.rectangle((x0, y0), (x1, y1), layer=1, datatype=0)


def grounds():
    # Identical, flat surroundings in every diagnostic. No raised pedestal,
    # internal ground island, regular cell, or regular-cap correction tabs.
    return [rect(0, 0, BOX_W, 448), rect(0, -1052, BOX_W, -618),
            rect(0, -618, 412, 0), rect(BOX_W-412, -618, BOX_W, 0)]


def cap(x, p):
    return [rect(x-p.w/2, -8-p.t_top, x+p.w/2, -8),
            rect(x-p.w/2, -8-p.l_left, x-p.w/2+p.t_left, -8),
            rect(x+p.w/2-p.t_right, -8-p.l_right, x+p.w/2, -8)]


def t1_inductor():
    m = T1.meander
    coil = place(make_meander(m, name='T1_COIL'), origin=(T1_X, -8-T1.t_top), anchor='top')
    # Local terminal landing occupies metal already inside the eventual cap.
    landing = rect(T1_X-m.t/2, -8-T1.t_top, T1_X+m.t/2, -8)
    return coil.get_polygons() + [landing]


def wire(width):
    if width <= 0 or width % 4:
        raise ValueError('Centered straight wire width must be a positive multiple of 4 um')
    # One uniform rectangle. Its top landing becomes part of the cap when
    # the cap is added. There is no widening or pedestal at its bottom.
    return [rect(T2_X-width/2, BOTTOM_GROUND, T2_X+width/2, -8)]


def comb():
    left = T2_X-T2['w']/2
    trunk = left-IDC['span']
    first_y = -8-T2['t_top']-IDC['d_top']
    pitch = IDC['t']+IDC['d']
    node, drive = [], []
    for row in range(2*IDC['pairs']):
        y = first_y-row*pitch
        if row % 2:
            node.append(rect(trunk+IDC['trunk_t']+IDC['d'], y-IDC['t'], left+T2['t_left'], y))
        else:
            drive.append(rect(trunk, y-IDC['t'], left-IDC['d'], y))
    bottom = first_y-(2*IDC['pairs']-1)*pitch-IDC['t']
    drive.append(rect(trunk, bottom, trunk+IDC['trunk_t'], -8))
    return node, drive, trunk+IDC['trunk_t']/2


def cases():
    p = SimpleNamespace(**T2)
    c1, c2 = cap(T1_X, T1), cap(T2_X, p)
    node_comb, drive_comb, drive_x = comb()
    result = []

    def add(name, title, metal, terminals, measurement, fixed, note):
        result.append(dict(name=name, title=title, metal=grounds()+metal,
                           terminals=terminals, measurement=measurement,
                           fixed=fixed, note=note))

    add('01_t1_inductor_only', 'T1 Inductor Only', t1_inductor(), [('T1', T1_X)],
        'L_t1', [], 'Cap removed except its narrow terminal landing; no T2, IDC or next-cell electrode.')
    add('02_t1_alone', 'T1 Alone: Cap + Original Meander', c1+t1_inductor(), [('T1', T1_X)],
        'C0_t1 after L_t1 is measured', ['L_t1 from 01'],
        'Exactly the selected T1 unit metal, without any neighbouring resonator or coupling electrode.')
    for i, width in enumerate(WIRE_WIDTHS, start=3):
        add(f'{i:02d}_t2_wire_w{width}', f'T2 Straight Inductor: {width} um', wire(width), [('T2', T2_X)],
            'L_t2 versus trace width', [], 'Uniform straight wire to bottom ground. Width is the only variable in this family.')
    add('06_t2_cap_only', 'T2 Cap Electrodes Only', c2, [('T2', T2_X)],
        'Electrode-to-ground capacitance diagnostic', [],
        'No wire or IDC. This is electrode capacitance, not automatically the effective C0 of the full LC.')
    add(f'07_t2_alone_w{INTEGRATION_WIRE}', f'T2 Alone: Cap + Straight {INTEGRATION_WIRE} um Wire', c2+wire(INTEGRATION_WIRE), [('T2', T2_X)],
        'C0_t2 after L_t2 is measured', ['L_t2 from matching bare-wire test'],
        f'No T1, IDC, pedestal or raised ground. The {INTEGRATION_WIRE} um width is a starting guess, not a calibrated result.')
    add('08_input_idc_only', 'Input IDC: No Inductors or T1', c2+node_comb+drive_comb,
        [('input', drive_x), ('T2', T2_X)], 'Cg_t2 from mutual terminal admittance', [],
        'Measure off-diagonal capacitive coupling directly; retain diagonal self-capacitances as diagnostics.')
    add('09_t2_t1_coupling_only', 'T2-to-T1 Coupling: Caps Only', c2+c1,
        [('T2', T2_X), ('T1', T1_X)], 'Cg_t1 from mutual terminal admittance', [],
        'Six-micron gap, no IDC and neither inductor. Only the selected T1 cap is retained in this diagnostic.')
    add('10_t2_t1_with_idc', 'Integration Preview: IDC + T2 + Unchanged T1',
        c2+wire(INTEGRATION_WIRE)+node_comb+drive_comb+c1+t1_inductor(),
        [('input', drive_x), ('T2', T2_X), ('T1', T1_X)],
        'Consistency check after separate extraction; no simultaneous free L/C fit',
        ['Previously measured L and capacitance terms'],
        'Not a complete filter. No regular cell or output feed; these are internal terminal ports for staged diagnostics.')
    return result


def digest(path):
    return sha256(Path(path).read_bytes()).hexdigest()


def signature(q):
    return (q.tag, dict(q.attrib), (q.text or '').strip(), [signature(x) for x in q])


def native_read_only(command, log):
    """Allow import/read commands only; em.exe is explicitly unreachable."""
    name = Path(command[0]).name.lower()
    if name != 'gds.exe' and not (name == 'soncmd.exe' and str(command[1]).lower() == '-read'):
        raise ValueError('Preparation permits only gds.exe import and soncmd.exe -Read; no simulations')
    return run(command, log)


def prepare(item):
    dest, gd = NATIVE/item['name'], GDS_ROOT/item['name']
    dest.mkdir(parents=True, exist_ok=True)
    gd.mkdir(parents=True, exist_ok=True)
    polygons = item['metal']
    if not all(np.allclose(p.points/2, np.round(p.points/2)) for p in polygons):
        raise ValueError('Staged source geometry must already lie on the 2 um grid')
    gds = gd/'geometry.gds'
    lib = gdstk.Library(unit=1e-6, precision=1e-9)
    lib.new_cell('TOP').add(*polygons)
    lib.write_gds(str(gds))
    # Translate into the positive GDS box convention expected by gds.exe.
    lib = gdstk.Library(unit=1e-6, precision=1e-9)
    lib.new_cell('TOP').add(*[q.copy().translate(0, 1052) for q in polygons])
    imported_gds = dest/'geometry_grid2.gds'
    lib.write_gds(str(imported_gds))
    tree = ET.parse(SOURCE)
    metal = tree.find(f'.//{tag("PlanarTechLayer")}')
    metal.insert(0, ET.Element(tag('ImportExportMapping'), GDSStreamDataTypes='1:0'))
    mapped = dest/'mapped_template.sonx'
    write_xml(tree, mapped)
    layer = dest/'superconductor.lay'
    layer.write_text('Metal-1 metal_boundary 1 0 0 "SUPERCONDUCTOR" on 0\n')
    project = dest/'staged_test.sonx'
    native_read_only([BIN/'gds.exe', '-v', '-i', mapped, '-g', '-l', layer, '-p', 'TOP',
         '-u', 'um', '-s', f'{BOX_W},{BOX_H}', '-t', 'min', '-C', '2,2',
         '-o', project, imported_gds], dest/'import.log')
    tree = ET.parse(project)
    root, reference = tree.getroot(), ET.parse(SOURCE)
    geometry = root.find(tag('Geometry'))
    actual = [gdstk.Polygon(points(q)) for q in geometry.iter(tag('PlanarPolygon'))]
    expected = [gdstk.Polygon([(x, 448-y) for x,y in q.points]) for q in polygons]
    if gdstk.boolean(actual, expected, 'xor', precision=.001):
        raise ValueError('Imported staged geometry differs from source')
    if not all(np.allclose(q.points/2, np.round(q.points/2)) for q in actual):
        raise ValueError('Imported geometry is off grid')
    if list(geometry.iter(tag('Port'))) or list(geometry.iter(tag('Component'))):
        raise ValueError('Importer unexpectedly retained template ports/components')
    box = geometry.find(tag('Box'))
    for q in list(box.findall(tag('RefPlane'))):
        box.remove(q)  # No box-wall feed ports or 400 um feed de-embedding.
    objects = next(q for q in geometry.findall(tag('Level')) if q.get('Number')=='0').find(tag('Objects'))
    # Use the same native component-pin structure as the established C_add
    # fixture, converted to exposed ports with no ideal component value.
    cap_source = ET.parse(ROOT/'sonnet_files/filter_v6.sonx')
    prototype = next(q for q in cap_source.iter(tag('Component')) if q.findtext(tag('Value'))=='C_add')
    maxid = max(int(q.get('Id')) for q in geometry.iter(tag('PlanarPolygon')))
    terminal_records = []
    for i, (label, x) in enumerate(item['terminals']):
        if gdstk.inside([(x,447), (x,449), (x,455), (x,457)], actual) != (True,False,False,True):
            raise ValueError(f'Invalid top-gap terminal for {label}')
        component = deepcopy(prototype)
        component.set('Type', 'PORTS'); component.set('Name', label+'_terminal')
        component.set('Id', str(maxid+i+1)); component.set('MacroID', f'terminal{i}')
        component.attrib.pop('Units', None)
        for q in list(component.findall(tag('Value'))): component.remove(q)
        component.find(tag('LabelPos')).text = f'({x+8},452)'
        component.find(tag('SymbolBox')).text = f'({x-2},450,{x+2},454)'
        for j, pin in enumerate(component.findall(tag('Pin'))):
            pin.set('MacroID', f'terminal{i}Pin{j}'); pin.set('PortNum', str(2*i+j+1))
            pin.find(tag('Center')).text = f'({x},{448 if j==0 else 456})'
        objects.append(component)
        terminal_records.append(dict(node=label, ground_port=2*i+1, signal_port=2*i+2,
                                     native_ground_xy=[x,448], native_signal_xy=[x,456]))
    # Preparation settings only. These broad diagnostic frequencies have NOT
    # been simulated, and no fitting/extraction is performed by this script.
    sweeps = root.find(tag('Sweeps'))
    for q in list(sweeps): sweeps.remove(q)
    fs = ET.SubElement(ET.SubElement(sweeps, tag('Set'), On='TRUE'), tag('Frequencies'))
    for f in (.1, .5, 1, 2, 4, 5, 6):
        ET.SubElement(fs, tag('Sweep'), Type='SWEEP', On='TRUE', Adaptive='FALSE',
                      Start=str(f), Stop=str(f), Step='1')
    output = root.find(tag('OutputFiles')); output.set('Folder', '.')
    for q in output.iter(tag('File')):
        if q.get('Type')=='TOUCH':
            q.set('Name', 'staged_test.snp'); q.set('UseInterpolatedData', 'FALSE')
    for goal in root.findall(f'{tag("Optimize")}/{tag("Goal")}'): goal.set('On','FALSE')
    for name in ('Dielectric','Conductor','DielectricMaterialModel','BoxCover'):
        if [signature(q) for q in geometry.iter(tag(name))] != [signature(q) for q in reference.find(tag('Geometry')).iter(tag(name))]:
            raise ValueError(f'Inherited {name} changed')
    write_xml(tree, project)
    native_read_only([BIN/'soncmd.exe', '-Read', project], dest/'sonnet_read.log')
    metadata = {k:v for k,v in item.items() if k not in ('metal','terminals')}
    metadata.update(status='inspection_only_not_simulated', project=str(project),
                    geometry_gds=str(gds), box_um=[BOX_W,BOX_H], mesh_um=[2,2],
                    T1=asdict(T1), T2=T2, IDC=IDC, wire_widths_um=list(WIRE_WIDTHS),
                    integration_wire_width_um=INTEGRATION_WIRE,
                    terminal_pairs=terminal_records, specification=specification(),
                    geometry_xor_empty=True, material_stack_preserved=True,
                    source_template_sha256=digest(SOURCE), project_sha256=digest(project),
                    geometry_sha256=digest(gds), solver_run=False,
                    raised_ground_box=False, external_box_ports=False)
    (dest/'inspection.json').write_text(json.dumps(metadata,indent=2)+'\n')
    preview(item, dest/'geometry_review.png')
    return metadata


def preview(item, path, ax=None):
    standalone = ax is None
    if standalone: fig, ax = plt.subplots(figsize=(12,7))
    for q in item['metal']:
        ax.add_patch(Polygon(q.points, facecolor='#388563', edgecolor='none'))
    active = item['metal'][4:]
    bounds = np.array([q.bounding_box() for q in active])
    left, right = bounds[:,0,0].min()-65, bounds[:,1,0].max()+65
    for i, (label,x) in enumerate(item['terminals']):
        ax.plot([x,x], [0,-8], 'o--', color='#b51f32', ms=3, lw=1)
        ax.annotate(f'{label}: P{2*i+1} ground / P{2*i+2} signal',
                    (x,0), (x,24), ha='center', fontsize=8,
                    arrowprops=dict(arrowstyle='-', color='#b51f32'))
    ax.set(xlim=(max(0,left),min(BOX_W,right)), ylim=(-650,60),
           xlabel='x (um)', ylabel='y (um)', title=item['title'])
    ax.set_aspect('equal'); ax.tick_params(labelsize=8)
    if standalone:
        fig.tight_layout(); fig.savefig(path,dpi=180); plt.close(fig)


def main():
    if NATIVE.exists() and any(NATIVE.rglob('*.snp')):
        raise ValueError('This inspection folder now contains simulation data; choose a new NAME before regenerating')
    before = digest(SOURCE)
    records = [prepare(item) for item in cases()]
    if digest(SOURCE) != before: raise ValueError('Source template changed')
    (NATIVE/'manifest.json').write_text(json.dumps(dict(status='inspection_only_not_simulated',
        solver_run=False, cases=records),indent=2)+'\n')
    fig, axes = plt.subplots(2,5,figsize=(24,10))
    for item, ax in zip(cases(), axes.flat): preview(item,None,ax)
    fig.suptitle('Staged Taper Tests — Geometry Review Only; No EM Simulations',fontsize=20)
    fig.tight_layout(); fig.savefig(NATIVE/'all_tests.png',dpi=160); plt.close(fig)
    print(f'Prepared {len(records)} inspection fixtures. No EM simulation launched.\n{NATIVE}',flush=True)


if __name__=='__main__': main()
