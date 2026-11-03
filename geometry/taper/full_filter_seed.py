"""Add one deliberately selected full-geometry EM frequency to a preserved cache."""
from pathlib import Path
from datetime import datetime,timezone
import json,shutil,re,time,hashlib
import xml.etree.ElementTree as ET
from geometry.taper.full_filter import NATIVE,OUT,PROJECT_NAME
from geometry.taper.full_filter_cache_bands import SNAPSHOT,canonical
from geometry.prepare_sonnet import tag,write_xml
from geometry.taper.sonnet import BIN,guarded_solver,solver_environment
from circuit_common import read_snp
import numpy as np

FREQUENCY_GHZ=4.8375

def run(frequency=FREQUENCY_GHZ):
    root=OUT/'solver_efficiency_review'/('seed_'+str(frequency).replace('.','p'))
    if root.exists():raise FileExistsError(root)
    root.mkdir(parents=True)
    case=root/'native';shutil.copytree(SNAPSHOT,case)
    project=case/PROJECT_NAME
    source=ET.parse(NATIVE/PROJECT_NAME).getroot();tree=ET.parse(project)
    for b in tree.findall(f'{tag("Sweeps")}/{tag("Set")}'):
        if b.get('On')=='TRUE':
            for sw in b.findall(f'{tag("Frequencies")}/{tag("Sweep")}'):
                if sw.get('On')=='TRUE':
                    sw.set('Type','SWEEP');sw.set('Adaptive','FALSE');sw.set('Start',str(frequency));sw.set('Stop',str(frequency));sw.set('Step','0.0001');sw.attrib.pop('Target',None)
    write_xml(tree,project)
    check=ET.parse(project).getroot()
    assert canonical(check.find(tag('Geometry')))==canonical(source.find(tag('Geometry')))
    assert canonical(check.find(tag('Control')))==canonical(source.find(tag('Control')))
    base=[str(BIN/'em.exe'),'-v','-AbsCacheMultiSweep','-SubFreqHz5400000000']
    with (root/'preflight.log').open('w') as log:r=guarded_solver(base+['-N',str(project)],cwd=case,stdout=log,env=solver_environment())
    r.check_returncode();pre=(root/'preflight.log').read_text(encoding='utf-8-sig')
    match=re.search(r'requires\s+([\d,]+)\s+subsections',pre,re.I)
    assert match and int(match.group(1).replace(',',''))==40959,pre
    record=dict(status='running',phase='one_explicit_frequency',started_utc=datetime.now(timezone.utc).isoformat(),requested_frequency_GHz=frequency,source_cached_frequencies=15,project=str(project),project_sha256=hashlib.sha256(project.read_bytes()).hexdigest(),geometry_and_control_unchanged=True,grid_um=2,subsectioning_frequency_GHz=5.4,subsections=40959,threads=4,original_checkpoint_snapshot=str(SNAPSHOT))
    status=root/'run_status.json';status.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
    plan_path=Path('geometry/taper/full_filter_work_plan.json');plan=json.loads(plan_path.read_text(encoding='utf-8-sig'))
    plan.update(full_state='running_explicit_seed',full_phase='one_targeted_lower_band_field',full_active_output=str(root),full_active_project=str(project),full_current_frequency_GHz=frequency,full_completed_field_frequencies_before_seed=15,full_cached_checkpoint_preserved=True)
    plan_path.write_text(json.dumps(plan,indent=2)+'\n',encoding='utf-8')
    print('FULL_FILTER_SEED_STARTED',json.dumps(record),flush=True)
    started=time.monotonic()
    try:
        with (root/'em.log').open('w') as log:r=guarded_solver(base+[str(project)],cwd=case,stdout=log,env=solver_environment())
        r.check_returncode()
        text=(root/'em.log').read_text(encoding='utf-8-sig')
        assert 'Em simulation completed' in text and 'Em errors detected: 0' in text,text[-2000:]
        assert re.search(r'Total time for 1 frequenc(?:y|ies):', text),text[-2000:]
        net=read_snp(project.with_suffix('.snp'))
        found=np.flatnonzero(np.isclose(net.f,frequency*1e9,rtol=0,atol=1));assert len(found)==1
        k=int(found[0]);record.update(status='complete',ended_utc=datetime.now(timezone.utc).isoformat(),elapsed_seconds=time.monotonic()-started,seed_S21_magnitude=float(abs(net.s[k,1,0])),seed_S21_dB=float(20*np.log10(abs(net.s[k,1,0]))),exported_samples=len(net.f),cache_directory=str(case/'sondata'/project.stem))
        status.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8');print('FULL_FILTER_SEED_COMPLETE',json.dumps(record),flush=True)
    except BaseException as e:
        record.update(status='failed',reason=str(e));status.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8');raise
    return record

if __name__=='__main__':run()
