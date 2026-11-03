"""Sequential simple-layout T2 iteration, with an explicit corrected fit policy."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
from dataclasses import replace
from pathlib import Path
import json
import shutil
import numpy as np
import skrf
from circuit_common import read_snp
from geometry.taper.corrected_iteration import ROOT, OUT, review
from geometry.taper.second_layout import from_dict
from geometry.taper.compact_stages import prepare
from geometry.taper.compact_run import measure_L, response
from geometry.taper.sonnet import configure, simulate
from geometry.taper.run_bridged_stage import digest

NATIVE=ROOT/'sonnet_files/taper_exploration/corrected_iteration_2026-10-23'
BASE_NAME='plain_w500_wire52_idc26_r584_t24_gap6'


def baseline_broad(*,prepare_only=False):
    name='baseline_uniform_idc_broad'
    src=ROOT/'sonnet_files/taper_exploration/compact_centered'/BASE_NAME/'idc'
    dest=NATIVE/name
    out=OUT/name;out.mkdir(parents=True,exist_ok=True)
    project=dest/'generated_ncell_3.sonx'
    if not dest.exists():
        shutil.copytree(src,dest)
        configure(project,4,9.5,5500)
        (out/'intent.json').write_text(json.dumps(dict(source=str(src/'generated_ncell_3.sonx'),
            source_sha256=digest(src/'generated_ncell_3.sonx'), geometry_unchanged=True,
            reason='Cover the T2 response as well as the two lower resonances; narrow-band capacitor fits hit bounds.'),indent=2)+'\n')
    if prepare_only:return project
    marker=out/'run_status.json'
    if not marker.exists() or json.loads(marker.read_text()).get('status')!='complete':
        configure(project,4,9.5,5500)
        simulate(project,out)
    # Keep the already computed dense narrow peaks rather than requesting a
    # dense output grid over the entire broad band at every ABS iteration.
    dense=ROOT/'output_data/taper_exploration/compact_centered'/BASE_NAME/'idc/response.snp'
    broad=read_snp(project.with_suffix('.snp'));narrow=read_snp(dense)
    mask=(broad.f<narrow.f[0]) | (broad.f>narrow.f[-1])
    f=np.r_[broad.f[mask],narrow.f];s=np.concatenate([broad.s[mask],narrow.s]);order=np.argsort(f)
    combined=out/'combined_response.s2p'
    skrf.Network(frequency=skrf.Frequency.from_f(f[order],unit='hz'),s=s[order],z0=50).write_touchstone(str(combined))
    (out/'combined_sources.json').write_text(json.dumps(dict(broad=str(project.with_suffix('.snp')),
        broad_sha256=digest(project.with_suffix('.snp')),dense=str(dense),dense_sha256=digest(dense),
        samples=len(f),dense_interval_GHz=(narrow.f[[0,-1]]/1e9).tolist(),
        geometry_unchanged=True),indent=2)+'\n')
    return review(combined,out/'fit',.7341365873386546)


def pad_parameters(wire=44,pad=24,gap=12,span=174,pairs=26,junction=6):
    meta=ROOT/'geometry/taper/output/compact_centered'/BASE_NAME/'idc/geometry.json'
    p=from_dict(json.loads(meta.read_text())['second'])
    return replace(p, taper=replace(p.taper,meander=replace(p.taper.meander,t=wire,w=4*(3*wire//4+1))),
                   w_ground_pad=pad,l_ground_pad=586,d_ground_pad=gap,
                   idc_pairs=pairs,idc_span=span,d_t2_t1=junction)


def pads_calibrate(name,p,*,prepare_only=False):
    # The only additional metal is two bottom-ground pads alongside the wire.
    project,out,meta=prepare('corrected_'+name,p,'idc')
    OUT.mkdir(parents=True,exist_ok=True)
    (OUT/(name+'_intent.json')).write_text(json.dumps(dict(native=str(project),parameters=meta['second'],
        reason='Test whether two regular-cell-style side pads can raise T2 shunt C; straight wire stays full height.',
        T1_geometry_unchanged=True),indent=2)+'\n')
    if prepare_only:return project
    return measure_L('corrected_'+name,p,stage='idc_Lcal')


def combine_responses(broad_path,dense_path,out):
    """Merge identical-geometry sweeps; use dense data throughout its interval."""
    broad=read_snp(broad_path);dense=read_snp(dense_path)
    mask=(broad.f<dense.f[0]) | (broad.f>dense.f[-1])
    f=np.r_[broad.f[mask],dense.f];s=np.concatenate([broad.s[mask],dense.s]);ix=np.argsort(f)
    out.mkdir(parents=True,exist_ok=True);path=out/'combined_response.s2p'
    skrf.Network(frequency=skrf.Frequency.from_f(f[ix],unit='hz'),s=s[ix],z0=50).write_touchstone(str(path))
    mapping={round(v):i for i,v in enumerate(dense.f)}
    pairs=[(i,mapping[round(v)]) for i,v in enumerate(broad.f) if round(v) in mapping]
    if len(pairs)>10:
        a,b=np.array(pairs).T;delta=abs(broad.s[a,1,0]-dense.s[b,1,0])
        check=dict(shared_samples=len(a),overlap_S21_RMS=float(np.sqrt(np.mean(delta**2))),
                   overlap_S21_max=float(delta.max()))
        if check['overlap_S21_RMS']>.01 or check['overlap_S21_max']>.05:
            raise ValueError(f'Broad/dense exports disagree: {check}')
    else:check=dict(shared_samples=len(pairs))
    record=dict(broad=str(broad_path),broad_sha256=digest(broad_path),dense=str(dense_path),
                dense_sha256=digest(dense_path),samples=len(f),overlap=check)
    (out/'combined_sources.json').write_text(json.dumps(record,indent=2)+'\n')
    return path


def dense_cached_export(project,out):
    """Export the narrow peaks from the completed broad cache, with no new fields."""
    dest=project.parent/'dense_cached_export';target=dest/project.name
    result=out/'dense_cached_export';result.mkdir(parents=True,exist_ok=True)
    if not dest.exists():
        # Copy only files/cache from the primary native folder, avoiding recursion.
        dest.mkdir()
        for child in project.parent.iterdir():
            if child==dest:continue
            if child.name=='sondata':shutil.copytree(child,dest/child.name)
            elif child.is_file():shutil.copy2(child,dest/child.name)
        configure(target,4.65,5.4,7500)
    marker=result/'run_status.json'
    if not marker.exists() or json.loads(marker.read_text()).get('status')!='complete':
        simulate(target,result,cache_only=True)
    return target.with_suffix('.snp')


def pads_response(name,p,*,broad=False):
    lpath=ROOT/'output_data/taper_exploration/compact_centered'/('corrected_'+name)/'idc_Lcal/results.json'
    cal=json.loads(lpath.read_text());l=cal['L0_nH']
    if abs(l/.74-1)>.02 or cal['R_squared']<.9999:
        raise ValueError(f'Adjust the straight wire before the full response: independent L={l}')
    snp,native_out=response('corrected_'+name,p,stage='idc',
                            start=4 if broad else 4.65,stop=9.5 if broad else 5.4,
                            target=5500 if broad else 7500)
    folder=OUT/name
    if broad:
        dense=dense_cached_export(snp.with_suffix('.sonx'),native_out)
        snp=combine_responses(snp,dense,folder)
    result=review(snp,folder/'fit',l)
    from geometry.taper.assess_corrected import assess
    assess(folder/'fit')
    return result


if __name__=='__main__':
    import sys
    if sys.argv[1]=='baseline-broad':baseline_broad(prepare_only='--prepare-only' in sys.argv)
    elif sys.argv[1]=='pads-calibrate':
        pads_calibrate('pads24_d12_wire44',pad_parameters(),prepare_only='--prepare-only' in sys.argv)
    elif sys.argv[1]=='pads-response':pads_response('pads24_d12_wire44',pad_parameters())
    else:raise ValueError('Unknown action')
