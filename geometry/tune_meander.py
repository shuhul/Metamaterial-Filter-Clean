"""Tune horizontal meander run lengths on isolated native filter_v6 copies.

in_width is our meander.w. Vertical height, 72 runs, trace/gap, ports,
materials and original 2 um mesh are preserved. Trial reference planes
are 400 um to stay inside uniform feeds (original 600 um generates warnings). Each trial archives a full
C_add sweep, extracted L0, capacitance fit, solver logs and plots.
"""
from pathlib import Path
import argparse
from datetime import datetime, timezone
from hashlib import sha256
import json
import subprocess
import sys
import xml.etree.ElementTree as ET
import numpy as np

BASE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(BASE))
import constants as cfg
import extract_generated_inductance as extraction
from geometry.prepare_sonnet import BIN, tag, points, write_xml
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon

SOURCE = cfg.SONNET_FILES_DIR / 'filter_v6.sonx'
PROJECTS = cfg.SONNET_FILES_DIR / 'meander_tuning_ref400'
RESULTS = cfg.OUTPUT_DATA_DIR / 'meander_tuning_ref400'


def digest(path):
    return sha256(path.read_bytes()).hexdigest()


def signature(node):
    return (node.tag, dict(node.attrib), (node.text or '').strip(), [signature(c) for c in node])


def prepare(width):
    folder = PROJECTS / f'w{width:g}'
    folder.mkdir(parents=True, exist_ok=True)
    project = folder / f'filter_v6_w{width:g}.sonx'
    if project.exists():
        raise FileExistsError(f'Refusing to overwrite existing trial: {project}')
    parameter_file = folder / 'parameters.txt'
    parameter_file.write_text(f'in_width={width:g}\n', encoding='ascii')
    command = [str(BIN/'soncmd.exe'), '-ParamFile', str(parameter_file), '-ReadWrite', str(SOURCE), str(project)]
    process = subprocess.run(command, capture_output=True, text=True)
    (folder/'native_copy.log').write_text(process.stdout+process.stderr)
    process.check_returncode()
    original, tree = ET.parse(SOURCE), ET.parse(project)
    root = tree.getroot()
    old = {p.get('Id'): points(p) for p in original.findall('.//'+tag('PlanarPolygon'))}
    new = {p.get('Id'): points(p) for p in root.findall('.//'+tag('PlanarPolygon'))}
    var = next(v for v in original.findall('.//'+tag('Variable')) if v.get('Name')=='in_width')
    linked = {r.get('Id') for r in var.findall('.//'+tag('Reference'))}
    if old.keys() != new.keys():
        raise ValueError('Native copy changed polygon IDs')
    changed = []
    for ident, before in old.items():
        after = new[ident]
        if before.shape != after.shape:
            raise ValueError('Unexpected vertex-count change')
        if not np.allclose(before, after, atol=1e-6, rtol=0):
            if ident not in linked or not np.allclose(before[:,1],after[:,1],atol=1e-6,rtol=0):
                raise ValueError('Native parameter changed unrelated geometry')
            changed.append(ident)
    if width != 190 and not changed:
        raise ValueError('Parameter did not move geometry')
    for ident in linked:
        if not np.isclose(np.ptp(new[ident][:,0]), width, atol=.001):
            raise ValueError(f'Meander {ident} does not have requested span')
    for name in ('Box','Dielectric','Conductor','DielectricMaterialModel','Port','Component'):
        if [signature(n) for n in original.findall('.//'+tag(name))] != [signature(n) for n in root.findall('.//'+tag(name))]:
            raise ValueError(f'Native copy changed {name}')
    # Same EM frequency range; denser interpolated output for tracking narrow peaks.
    for sweep_set in root.findall(f'{tag("Sweeps")}/{tag("Set")}'):
        if sweep_set.get('On')=='TRUE':
            sweep, = sweep_set.findall(f'{tag("Frequencies")}/{tag("Sweep")}')
            sweep.set('Start','4.35'); sweep.set('Stop','4.95'); sweep.set('Target','500')
    for plane in root.findall(f'.//{tag("Box")}/{tag("RefPlane")}'):
        plane.set('RefLength','400')
    root.find(tag('OutputFiles')).set('Folder','.')
    for export in root.findall(f'{tag("OutputFiles")}//{tag("File")}'):
        if export.get('Type')=='TOUCH':
            export.set('Name',project.stem+'.snp')
    for goal in root.findall(f'{tag("Optimize")}/{tag("Goal")}'):
        goal.set('On','FALSE')
    write_xml(tree,project)
    fig, ax = plt.subplots(figsize=(12,8))
    for xy in new.values():
        ax.add_patch(Polygon(xy,facecolor='seagreen'))
    ax.set(xlim=(-20,1920),ylim=(1520,-20),aspect='equal',xlabel='x (um)',ylabel='Sonnet y (um)',
           title=f'Original filter_v6 copy: horizontal meander span {width:g} um')
    fig.tight_layout(); fig.savefig(folder/'geometry.png',dpi=180); plt.close(fig)
    (folder/'manifest.json').write_text(json.dumps({'source':str(SOURCE),'source_sha256':digest(SOURCE),
        'parameter':'in_width','meander_w_um':width,'changed_polygon_ids':changed,
        'mesh_um':[2,2],'reference_planes_um':[400,400],
        'note':'Original filter_v6 copy; not the generated single-cell geometry.'},indent=2)+'\n')
    return project


def trial(width, source_hash):
    output = RESULTS/f'w{width:g}'
    result_path = output/'results.json'
    if result_path.exists():
        status = json.loads((output/'run_status.json').read_text())
        if status['state']!='complete' or status['source_sha256']!=source_hash or status['results_sha256']!=digest(result_path):
            raise ValueError('Existing trial is incomplete, changed or from a different source')
        return json.loads(result_path.read_text())
    project = PROJECTS/f'w{width:g}'/f'filter_v6_w{width:g}.sonx'
    if project.exists():
        manifest=json.loads((project.parent/'manifest.json').read_text())
        if manifest['source_sha256']!=source_hash or manifest['meander_w_um']!=width:
            raise ValueError('Cannot resume a trial with different inputs')
    else:
        project = prepare(width)
    output.mkdir(parents=True,exist_ok=True)
    status = {'state':'running','started_utc':datetime.now(timezone.utc).isoformat(),
              'source_sha256':source_hash,'meander_w_um':width,'project':str(project)}
    status_path = output/'run_status.json'
    status_path.write_text(json.dumps(status,indent=2)+'\n')
    extraction.FOLDER, extraction.PROJECT = project.parent, project
    extraction.SNP, extraction.OUTPUT = project.with_suffix('.snp'), output
    original_simulate=extraction.simulate
    def cached_or_simulate(capacitance,label):
        snapshot=output/f'{label}.sonx'
        snp=output/f'{label}.snp'
        log=output/f'{label}.log'
        if snapshot.exists() and snp.exists() and log.exists():
            log_text=log.read_text()
            header=snp.read_text(encoding='utf-8')
            import re
            match=re.search(r'^!< C_add = (\S+)',header,re.MULTILINE)
            if 'Em simulation completed' in log_text and 'Em errors detected: 0' in log_text and match and np.isclose(float(match[1]),capacitance,rtol=1e-9,atol=1e-12):
                extraction.write_sonx({'C_add':str(capacitance)},project)
                import shutil
                shutil.copyfile(snp,extraction.SNP)
                net=extraction.read_snp(snp)
                print(f'Reusing archived {label}',flush=True)
                return net,extraction.single_sweep_resonance(net)
        return original_simulate(capacitance,label)
    extraction.simulate=cached_or_simulate
    try:
        extraction.main(fit_capacitance=False)
        result = json.loads(result_path.read_text())
        result.update(meander_w_um=width,project=str(project),source_sha256=source_hash)
        result_path.write_text(json.dumps(result,indent=2)+'\n')
    except Exception as error:
        status.update(state='failed',error=str(error)); raise
    else:
        status.update(state='complete',results_sha256=digest(result_path))
    finally:
        extraction.simulate=original_simulate
        status['updated_utc']=datetime.now(timezone.utc).isoformat()
        status_path.write_text(json.dumps(status,indent=2)+'\n')
    if digest(SOURCE)!=source_hash:
        raise ValueError('Source project changed during tuning')
    return result


def summarize(records, source_hash):
    best = min(records,key=lambda r:abs(r['L0_nH']-cfg.TARGET_L0_NH))
    data = {'target_L0_nH':cfg.TARGET_L0_NH,'source_sha256':source_hash,'best':best,'trials':records,
            'note':'Horizontal run-length tuning on original filter_v6 copies; original 2 um mesh retained.'}
    (RESULTS/'summary.json').write_text(json.dumps(data,indent=2)+'\n')
    fig,ax=plt.subplots(figsize=(7,5))
    ordered=sorted(records,key=lambda r:r['meander_w_um'])
    ax.plot([r['meander_w_um'] for r in ordered],[r['L0_nH'] for r in ordered],'o-')
    ax.axhline(cfg.TARGET_L0_NH,color='red',linestyle='--',label='Target 4.9194 nH')
    ax.set(xlabel='Horizontal meander span (um)',ylabel='Extracted L0 (nH)',title='filter_v6 meander tuning')
    ax.grid(alpha=.3); ax.legend(); fig.tight_layout(); fig.savefig(RESULTS/'length_vs_inductance.png',dpi=200); plt.close(fig)
    print(f'BEST: w={best["meander_w_um"]:g} um; L0={best["L0_nH"]:.8f} nH; error={best["target_difference_percent"]:+.4f}%',flush=True)
    return best


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-trials',type=int,default=4)
    parser.add_argument('--tolerance-percent',type=float,default=.1)
    parser.add_argument('--refine-best',action='store_true')
    args=parser.parse_args()
    if args.max_trials<2 or args.tolerance_percent<=0:
        parser.error('Use at least two trials and a positive tolerance')
    RESULTS.mkdir(parents=True,exist_ok=True)
    source_hash=digest(SOURCE)
    records=[]; width=190
    for index in range(args.max_trials):
        print(f'TRIAL {index+1}: horizontal span {width:g} um',flush=True)
        records.append(trial(width,source_hash))
        best=summarize(records,source_hash)
        if len(records)>=2 and abs(best['target_difference_percent'])<=args.tolerance_percent:
            break
        if len(records)==1:
            width=186 if records[0]['L0_nH']>cfg.TARGET_L0_NH else 194
        else:
            a,b=sorted(records,key=lambda r:abs(r['L0_nH']-cfg.TARGET_L0_NH))[:2]
            slope=(a['L0_nH']-b['L0_nH'])/(a['meander_w_um']-b['meander_w_um'])
            if slope<=0:
                raise ValueError('Unexpected sensitivity: investigate before continuing')
            predicted=a['meander_w_um']+(cfg.TARGET_L0_NH-a['L0_nH'])/slope
            tried={r['meander_w_um'] for r in records}
            choices=[w for w in range(170,211,2) if w not in tried]
            width=min(choices,key=lambda w:abs(w-predicted))
    if args.refine_best:
        refine(best,source_hash)
    print(f'Summary: {RESULTS/"summary.json"}',flush=True)


def refine(best, source_hash):
    """Check the selected candidate at 0.4 MHz spacing in the narrow sweep."""
    output=RESULTS/'best_refined'
    output.mkdir(parents=True,exist_ok=True)
    project=Path(best['project'])
    tree=ET.parse(project)
    for sweep_set in tree.findall(f'{tag("Sweeps")}/{tag("Set")}'):
        if sweep_set.get('On')=='TRUE':
            for sweep in sweep_set.findall(f'{tag("Frequencies")}/{tag("Sweep")}'):
                sweep.set('Start','4.35'); sweep.set('Stop','4.95'); sweep.set('Target','1500')
    write_xml(tree,project)
    extraction.FOLDER,extraction.PROJECT=project.parent,project
    extraction.SNP,extraction.OUTPUT=project.with_suffix('.snp'),output
    status={'state':'running','meander_w_um':best['meander_w_um'],'source_sha256':source_hash,
            'started_utc':datetime.now(timezone.utc).isoformat()}
    status_path=output/'run_status.json'
    status_path.write_text(json.dumps(status,indent=2)+'\n')
    try:
        extraction.main(fit_capacitance=False)
        result_path=output/'results.json'
        result=json.loads(result_path.read_text())
        result.update(meander_w_um=best['meander_w_um'],project=str(project),source_sha256=source_hash)
        result_path.write_text(json.dumps(result,indent=2)+'\n')
        extraction.write_sonx({'C_add':'0.0001'},project)
        restored=ET.parse(project)
        for sweep_set in restored.findall(f'{tag("Sweeps")}/{tag("Set")}'):
            if sweep_set.get('On')=='TRUE':
                for sweep in sweep_set.findall(f'{tag("Frequencies")}/{tag("Sweep")}'):
                    sweep.set('Start','4'); sweep.set('Stop','6'); sweep.set('Target','500')
        write_xml(restored,project)
        archive=project.parent/'last_sweep_cadd_50ff.snp'
        if archive.exists():
            raise FileExistsError(archive)
        extraction.SNP.rename(archive)
        result.update(project_C_add_fF=.1,last_sweep_export=str(archive),
            note='Project restored to 0.1 fF and original 4-6 GHz/500 interval sweep; no final 0.1 fF export.',
            final_project_sha256=digest(project))
        result_path.write_text(json.dumps(result,indent=2)+'\n')
        summary_path=RESULTS/'summary.json'
        overview=json.loads(summary_path.read_text())
        overview['refined_best']=result
        summary_path.write_text(json.dumps(overview,indent=2)+'\n')
    except Exception as error:
        status.update(state='failed',error=str(error));raise
    else:
        status.update(state='complete',results_sha256=digest(result_path))
    finally:
        status['updated_utc']=datetime.now(timezone.utc).isoformat()
        status_path.write_text(json.dumps(status,indent=2)+'\n')
    if digest(SOURCE)!=source_hash:
        raise ValueError('Source changed during refinement')
    print(f'REFINED: w={result["meander_w_um"]:g} um; L0={result["L0_nH"]:.8f} nH; error={result["target_difference_percent"]:+.4f}%',flush=True)


if __name__=='__main__':
    main()
