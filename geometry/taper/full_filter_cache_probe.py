"""Zero-new-frequency ABS tests on isolated copies of completed cache data."""
from pathlib import Path
from datetime import datetime,timezone
import shutil,json,re,time,hashlib
import xml.etree.ElementTree as ET
from geometry.taper.full_filter import NATIVE,OUT,PROJECT_NAME
from geometry.taper.sonnet import BIN
from geometry.prepare_sonnet import tag,write_xml
from geometry.taper.solver_guard import _run
from geometry.taper.solver_resources import solver_environment,solver_lease

ROOT=OUT/'solver_efficiency_review'
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main(cache_source=None):
    source=NATIVE/PROJECT_NAME
    signature=digest(source)
    stamp=datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    root=ROOT/('cache_tests_'+stamp);root.mkdir(parents=True)
    snapshot=root/'snapshot';snapshot.mkdir()
    shutil.copy2(source,snapshot/PROJECT_NAME)
    src_cache=Path(cache_source) if cache_source else NATIVE/'sondata'/source.stem
    dst_cache=snapshot/'sondata'/source.stem;dst_cache.mkdir(parents=True)
    cache_before=digest(src_cache/'abs_cache_d.sid')
    for f in src_cache.rglob('*'):
        if f.is_file() and f.name!='lock.lck':
            dst=dst_cache/f.relative_to(src_cache);dst.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(f,dst)
    if cache_before!=digest(src_cache/'abs_cache_d.sid'):raise RuntimeError('Source cache changed during snapshot; retry after current frequency')
    threads=root/'one_thread.ini';threads.write_text('[Sonnet]\nThreadLimit=1\n',encoding='utf-8')
    report=dict(source_project=str(source),source_project_sha256=signature,cache_sha256=cache_before,started_utc=datetime.now(timezone.utc).isoformat(),policy='Only -AbsNoDiscrete on copied projects. No geometry, materials, reference planes or native mesh changes. One thread; exclusive resource lease. Original run must be stopped between completed frequencies first.',tests=[])
    for target in (7500,2000):
        folder=root/f'target_{target}';shutil.copytree(snapshot,folder)
        project=folder/PROJECT_NAME
        if target!=7500:
            tree=ET.parse(project)
            for b in tree.findall(f'{tag("Sweeps")}/{tag("Set")}'):
                if b.get('On')=='TRUE':
                    for sw in b.findall(f'{tag("Frequencies")}/{tag("Sweep")}'):
                        if sw.get('On')=='TRUE':sw.set('Target',str(target))
            write_xml(tree,project)
        env=solver_environment();env.update(SON_THREAD_CONTROL_INI=str(threads),OMP_NUM_THREADS='1',MKL_NUM_THREADS='1',SONNET_MAX_JOBS='1')
        command=[str(BIN/'em.exe'),'-v','-AbsCacheMultiSweep','-AbsNoDiscrete','-SubFreqHz5400000000',str(project)]
        started=time.monotonic()
        with solver_lease(max_jobs=1,wait_seconds=5):
            with (folder/'probe.log').open('w') as log:result=_run(command,cwd=folder,stdout=log,env=env)
        text=(folder/'probe.log').read_text(encoding='utf-8-sig')
        if re.search(r'Matrix fill time:|Matrix solve time:',text):raise RuntimeError('Unexpected field calculation in cache-only probe')
        entry=dict(target=target,elapsed_seconds=time.monotonic()-started,returncode=result.returncode,log=str(folder/'probe.log'),snp_exists=project.with_suffix('.snp').exists(),zero_new_frequencies_confirmed=bool(re.search(r'Total time for 0 frequencies:',text)),tail=text[-6000:])
        if entry['snp_exists']:
            from circuit_common import read_snp
            net=read_snp(project.with_suffix('.snp'));entry.update(export_samples=len(net.f),export_range_GHz=(net.f[[0,-1]]/1e9).tolist())
        report['tests'].append(entry)
        (root/'results.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
        print(json.dumps(entry),flush=True)
    report['source_project_unchanged']=digest(source)==signature
    (root/'results.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print('CACHE_PROBE_REPORT',str(root/'results.json'),flush=True)
    return report

if __name__=='__main__':main()
