"""Check narrower ABS bands using15 cached fields; prohibit new EM points."""
from pathlib import Path
import json,shutil,time,re,hashlib
import xml.etree.ElementTree as ET
from geometry.taper.full_filter_cache_probe import ROOT,PROJECT_NAME
from geometry.taper.sonnet import BIN
from geometry.prepare_sonnet import tag,write_xml
from geometry.taper.solver_guard import _run
from geometry.taper.solver_resources import solver_environment,solver_lease
from circuit_common import read_snp
import numpy as np

SNAPSHOT=ROOT/'cache_tests_20261101T005507Z/snapshot'
BANDS=[(4.65,5.025),(5.025,5.1203),(5.1203,5.2008),(5.2008,5.4)]
def canonical(e):return (e.tag,sorted(e.attrib.items()),(e.text or '').strip(),[canonical(c) for c in e])
def main(bands=BANDS,label='four_bands'):
    root=ROOT/label;root.mkdir(exist_ok=True)
    source=ET.parse(SNAPSHOT/PROJECT_NAME).getroot()
    geometry=canonical(source.find(tag('Geometry')));control=canonical(source.find(tag('Control')))
    threads=root/'one_thread.ini';threads.write_text('[Sonnet]\nThreadLimit=1\n',encoding='utf-8')
    results=[]
    for lo,hi in bands:
        target=max(50,round((hi-lo)/.0001));folder=root/f'{lo:.4f}_{hi:.4f}'
        if folder.exists():raise FileExistsError(folder)
        shutil.copytree(SNAPSHOT,folder);project=folder/PROJECT_NAME
        tree=ET.parse(project)
        for b in tree.findall(f'{tag("Sweeps")}/{tag("Set")}'):
            if b.get('On')=='TRUE':
                for sw in b.findall(f'{tag("Frequencies")}/{tag("Sweep")}'):
                    if sw.get('On')=='TRUE':sw.set('Start',str(lo));sw.set('Stop',str(hi));sw.set('Target',str(target))
        write_xml(tree,project)
        check=ET.parse(project).getroot();assert canonical(check.find(tag('Geometry')))==geometry;assert canonical(check.find(tag('Control')))==control
        env=solver_environment();env.update(SON_THREAD_CONTROL_INI=str(threads),OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
        cmd=[str(BIN/'em.exe'),'-v','-AbsCacheMultiSweep','-AbsNoDiscrete','-SubFreqHz5400000000',str(project)]
        start=time.monotonic()
        with solver_lease(max_jobs=1,wait_seconds=5):
            with (folder/'probe.log').open('w') as log:r=_run(cmd,cwd=folder,stdout=log,env=env)
        text=(folder/'probe.log').read_text(encoding='utf-8-sig');assert not re.search('Matrix fill time:|Matrix solve time:',text)
        snrs=re.findall(r'SNR = ([-\d.]+) dB \(([\d.]+) dB Goal\).*?at ([\d.]+) GHz',text)
        item=dict(band_GHz=[lo,hi],target=target,elapsed_seconds=time.monotonic()-start,returncode=r.returncode,geometry_and_control_unchanged=True,zero_new_frequencies_confirmed='Total time for 0 frequencies:' in text,log=str(folder/'probe.log'),snr_check=snrs[-1] if snrs else None,tail=text[-3500:])
        if project.with_suffix('.snp').exists():
            net=read_snp(project.with_suffix('.snp'));mask=(net.f>=lo*1e9-1)&(net.f<=hi*1e9+1)
            item.update(total_export_samples=len(net.f),in_band_samples=int(mask.sum()),export_range_GHz=(net.f[[0,-1]]/1e9).tolist())
            if mask.any():item['in_band_range_GHz']=(net.f[mask][[0,-1]]/1e9).tolist()
        item['criterion_pass']=bool(snrs and float(snrs[-1][0])>=float(snrs[-1][1]))
        results.append(item);(root/'results.json').write_text(json.dumps(results,indent=2)+'\n',encoding='utf-8');print(json.dumps(item),flush=True)
    return results

if __name__=='__main__':main()
