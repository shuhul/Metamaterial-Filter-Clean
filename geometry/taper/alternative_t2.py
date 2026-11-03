r"""New T2 campaign: isolated ordinary meander first, then capacitances.

Run from the project root with .venv\Scripts\python.exe -B -m
geometry.taper.alternative_t2. Edit WIDTH_UM for a new width trial.
This script never runs an assembled-filter simulation.
"""
import os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['OMP_NUM_THREADS'] = '1'
from dataclasses import asdict
from pathlib import Path
import json
import xml.etree.ElementTree as ET
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import gdstk
from geometry.meander import MeanderParams
from geometry.specifications import specification, ALTERNATIVE_TAPER_2
from geometry.taper.layout import TaperParams
from geometry.taper.second_layout import SecondParams, write, validate
from geometry.taper.locked_taper2 import T1
from geometry.prepare_sonnet import main as import_native, tag, write_xml
from geometry.taper.sonnet import add_capacitor
from geometry.taper.measure_compact_inductance import calibrate

ROOT = Path(__file__).resolve().parents[2]
CAMPAIGN = 'alternative_t2_1nH_10fF_2026-10-25'
WIDTH_UM = 32
RUN_SIMULATION = False  # Historical campaign, superseded by three_tapers.py.


def parameters(width):
    # Same square-turn 4 um trace style as the regular meander. Width is
    # the only variable between trials; no widened ground tail or pedestal.
    m = MeanderParams(w=width, l=600, l_top=10, n=58, t=4, d=6)
    return SecondParams(
        taper=TaperParams(w=128, t_top=10, t_left=4, t_right=4,
                          l_left=40, l_right=40, meander=m),
        taper1=T1, straight_inductor=False,
        idc_pairs=1, idc_t=4, idc_d=12, idc_span=80,
        d_idc_top=2, d_t2_t1=40, l_cal_output=40)


def run(width=WIDTH_UM, simulate=RUN_SIMULATION):
    p = parameters(width)
    case = f'meander_w{width}_n58_t4_d6'
    generated = ROOT/'geometry/taper/output'/CAMPAIGN/case
    native = ROOT/'sonnet_files/taper_exploration'/CAMPAIGN/case
    out = ROOT/'output_data/taper_exploration'/CAMPAIGN/case
    for folder in (generated, native, out):
        folder.mkdir(parents=True, exist_ok=True)
    path = generated/'geometry.gds'
    if path.exists():
        old = json.loads(path.with_suffix('.json').read_text())
        from geometry.taper.second_layout import from_dict
        if from_dict(old['second']) != p:
            raise ValueError('Existing case differs; use a new case name')
    else:
        write(p, path, isolation=True)
    meta = json.loads(path.with_suffix('.json').read_text())
    validate(path, meta)
    intent = {'specification': specification(), 'parameters': asdict(p),
              'stage': 'isolated ordinary meander L calibration only',
              'T1_present_in_EM': False, 'final_design': False,
              'note': 'Short weak calibration electrodes; not the final IDC or C0 geometry.'}
    (out/'intent.json').write_text(json.dumps(intent, indent=2)+'\n')
    fig, ax = plt.subplots(figsize=(10, 8))
    for q in gdstk.read_gds(str(path)).top_level()[0].get_polygons():
        ax.add_patch(Polygon(q.points, facecolor='#287b51', edgecolor='none'))
    x = meta['cell_centers'][1]
    ax.set(xlim=(x-160, x+250), ylim=(-640, 30), aspect='equal',
           xlabel='x (um)', ylabel='y (um)',
           title=f'Isolated T2 Meander: {width} um Width, 4 um Trace')
    fig.tight_layout(); fig.savefig(out/'geometry.png', dpi=200); plt.close(fig)
    project = native/'generated_ncell_1.sonx'
    if not project.exists():
        import_native(1, gds_path=path, output_dir=native, geometry_validator=validate)
        add_capacitor(project, meta)
        tree = ET.parse(project)
        comp = next(tree.iter(tag('Component')))
        comp.set('Type', 'PORTS'); comp.attrib.pop('Units', None)
        for q in list(comp):
            if q.tag == tag('Value'): comp.remove(q)
        sweeps = tree.find(tag('Sweeps'))
        for q in list(sweeps): sweeps.remove(q)
        fs = ET.SubElement(ET.SubElement(sweeps, tag('Set'), On='TRUE'), tag('Frequencies'))
        for f in (4, 5, 5.5, 6, 7):
            ET.SubElement(fs, tag('Sweep'), Type='SWEEP', On='TRUE',
                          Adaptive='FALSE', Start=str(f), Stop=str(f), Step='1')
        for q in tree.findall(f"{tag('OutputFiles')}//{tag('File')}"):
            if q.get('Type') == 'TOUCH': q.set('UseInterpolatedData', 'FALSE')
        write_xml(tree, project)
    print('Isolated meander project:', project, flush=True)
    if simulate:
        return calibrate(project, out, case, auto_loads=True, target_L_nH=ALTERNATIVE_TAPER_2.L0_nH)


if __name__ == '__main__':
    run()
