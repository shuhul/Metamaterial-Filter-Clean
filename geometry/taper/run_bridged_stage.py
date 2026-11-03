"""Run one reviewed feed-driven fixture, with provenance and guarded EM."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
from hashlib import sha256
from pathlib import Path
from geometry.taper.bridged_fixtures import ROOT,NATIVE,CASES
from geometry.taper.sonnet import configure,simulate

OUT=ROOT/'output_data/taper_exploration/bridged_iteration_2026-10-22'


def digest(path):return sha256(Path(path).read_bytes()).hexdigest()


def run(index,start=4.75,stop=5.3,target=5500):
    name=CASES[index-1]
    source=NATIVE/name/f'generated_ncell_{2 if index<3 else 3}.sonx'
    out=OUT/name;out.mkdir(parents=True,exist_ok=True)
    if not (out/'prepared_source.sonx').exists():
        manifest=json.loads((NATIVE/'manifest.json').read_text())
        expected=next(q for q in manifest['cases'] if q['case']==name)
        if digest(source)!=expected['project_sha256']:raise ValueError('Reviewed project changed before simulation')
        (out/'prepared_source.sonx').write_bytes(source.read_bytes())
    configure(source,start,stop,target)
    marker=out/'run_status.json'
    saved=json.loads(marker.read_text()) if marker.exists() else {}
    if saved.get('status')=='complete':
        if digest(source)!=saved['project_sha256'] or digest(source.with_suffix('.snp'))!=saved['snp_sha256']:
            raise ValueError('Saved stage source/export changed')
        print(f'Reusing verified completed stage: {source}',flush=True)
        return source.with_suffix('.snp'),out
    print(f'Running reviewed stage {index}: {source}',flush=True)
    snp=simulate(source,out)
    print(f'Completed: {snp}',flush=True)
    return snp,out


if __name__=='__main__':
    import sys
    run(int(sys.argv[1]))
