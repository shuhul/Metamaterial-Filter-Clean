"""Independent ordinary field at the frequency used for the node-probe pilot."""
import json
import xml.etree.ElementTree as ET
import numpy as np
from scipy.signal import find_peaks
from geometry.taper.known_capacitor_control import NATIVE,OUT,W,opened
from geometry.prepare_sonnet import tag,write_xml
from geometry.taper.run_bridged_stage import digest


def run(case='frame700_wire88_left214'):
    source=NATIVE/case/'idc/generated_ncell_3.sonx'
    ordinary=W['read_network'](source.with_suffix('.snp'),2)
    peaks=find_peaks(abs(ordinary.s[:,1,0]),prominence=.025)[0]
    i=int(max(peaks,key=lambda k:abs(ordinary.s[k,1,0])))
    frequency=float(ordinary.f[i]/1e9)
    dest=NATIVE/case/'ordinary_peak_direct_check';dest.mkdir(exist_ok=True)
    out=OUT/case/'ordinary_peak_direct_check';out.mkdir(exist_ok=True)
    project=dest/'ordinary_peak.sonx'
    if not project.exists():
        tree=ET.parse(source);root=tree.getroot();sweeps=root.find(tag('Sweeps'))
        for q in list(sweeps):sweeps.remove(q)
        fs=ET.SubElement(ET.SubElement(sweeps,tag('Set'),On='TRUE'),tag('Frequencies'))
        ET.SubElement(fs,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(frequency),Stop=str(frequency),Step='1')
        for q in root.iter(tag('File')):
            if q.get('Type')=='TOUCH':
                q.set('Name',project.with_suffix('.snp').name);q.set('UseInterpolatedData','FALSE')
        write_xml(tree,project)
    if ET.tostring(ET.parse(source).find(tag('Geometry')))!=ET.tostring(ET.parse(project).find(tag('Geometry'))):
        raise ValueError('Ordinary spot check changed geometry/settings')
    marker=out/'completed.json'
    if marker.exists():
        if json.loads(marker.read_text())!=dict(project_sha256=digest(project),snp_sha256=digest(project.with_suffix('.snp'))):
            raise ValueError('Completed spot data changed')
        net=W['read_network'](project.with_suffix('.snp'),2)
    else:
        net=W['em'](project,2)
        marker.write_text(json.dumps(dict(project_sha256=digest(project),snp_sha256=digest(project.with_suffix('.snp'))),indent=2)+'\n')
    error=abs(net.s[0]-ordinary.s[i])
    result=dict(frequency_GHz=frequency,source_sha256=digest(source),ordinary_snapshot_sha256=digest(source.with_suffix('.snp')),
        ordinary_direct_S21_magnitude=float(abs(net.s[0,1,0])),ABS_S21_magnitude=float(abs(ordinary.s[i,1,0])),
        ABS_vs_direct_S21_error=float(error[1,0]),ABS_vs_direct_all_S_max=float(error.max()),probe_comparisons=[])
    for name in ('known_capacitor_control','known_capacitor_control_floating'):
        path=OUT/case/name/'pilot_data.snp'
        probe=W['read_network'](path,4)
        if not np.allclose(probe.f,net.f,rtol=0,atol=1):raise ValueError('Pilot frequencies differ')
        delta=abs(opened(probe.s)[0]-net.s[0])
        result['probe_comparisons'].append(dict(reference=name,S21_error_vs_ordinary_direct=float(delta[1,0]),all_S_max=float(delta.max())))
    (out/'comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result


if __name__=='__main__':run()
