"""Preserve narrow-band fixtures; reuse copied field caches in a wider-band clone."""
import json
from shutil import copytree
from geometry.taper.progressive_fixtures import NATIVE,OUT
from geometry.taper.sonnet import configure,simulate
from geometry.taper.run_bridged_stage import digest


def run(source_name,name,start=4,stop=9.5):
    source=NATIVE/source_name;dest=NATIVE/name;out=OUT/name
    if not dest.exists():copytree(source,dest)
    out.mkdir(parents=True,exist_ok=True)
    project=dest/'generated_ncell_2.sonx'
    manifest=out/'source.json'
    expected=dict(source=str(source),source_project_sha256=digest(source/project.name),start_GHz=start,stop_GHz=stop)
    if manifest.exists() and json.loads(manifest.read_text())!=expected:raise ValueError('Changed broadband source')
    manifest.write_text(json.dumps(expected,indent=2)+'\n')
    configure(project,start,stop,round((stop-start)*1000))
    marker=out/'run_status.json'
    if marker.exists():
        r=json.loads(marker.read_text())
        if r['status']=='complete':
            if r['project_sha256']!=digest(project) or r['snp_sha256']!=digest(project.with_suffix('.snp')):
                raise ValueError('Changed completed broadband case')
            return project.with_suffix('.snp'),out
    return simulate(project,out),out


if __name__=='__main__':
    import sys
    run(sys.argv[1],sys.argv[2])
