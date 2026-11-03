"""One independent resonant-frequency check of the plain series operator."""
import json,shutil
import xml.etree.ElementTree as ET
import numpy as np
from geometry.taper.plain_series import ROOT,NATIVE,OUT,W,G,A,tag,write_xml,digest,analyze
from circuit_common import read_snp


def run(name='w500_top10_wire52_gap6'):
    p=NATIVE/name/'series_probes/series_probe.sonx';out=OUT/name/'series_probes'
    marker=out/'peak_refinement.json'
    if marker.exists():return json.loads(marker.read_text())
    old=np.load(out/'coupled_model.npz');freq=5.1339e9
    co=old['coefficients'];powers=old['powers']
    oldz=((np.array([freq])[:,None]/5e9)**powers@co).reshape(1,6,6)
    archive=out/'before_peak_refinement';archive.mkdir(exist_ok=True)
    for path in (p,p.with_suffix('.snp'),out/'coupled_model.npz',out/'analysis.json',out/'completed.json'):
        if not (archive/path.name).exists():shutil.copyfile(path,archive/path.name)
    tree=ET.parse(p);fs=tree.find(f"{tag('Sweeps')}/{tag('Set')}/{tag('Frequencies')}")
    if not any(abs(float(q.get('Start'))-freq/1e9)<1e-9 for q in fs):
        ET.SubElement(fs,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(freq/1e9),Stop=str(freq/1e9),Step='1')
    write_xml(tree,p);n=W['em'](p,6);i=np.argmin(abs(n.f-freq))
    shorted=G['close_wires'](n.s)[i]
    polarity=np.diag([1,1,-1,-1,-1,-1]);z=A['s_to_z'](polarity@n.s[[i]]@polarity)
    ref=read_snp(OUT/name/'direct/response.snp');j=np.argmin(abs(ref.f-freq))
    report=dict(frequency_GHz=freq/1e9,one_new_frequency=True,
        actual_probe_short_vs_ordinary_S21_error=float(abs(shorted[1,0]-ref.s[j,1,0])),
        actual_probe_short_vs_ordinary_max_all_S_error=float(abs(shorted-ref.s[j]).max()),
        interpolation_vs_new_probe_S21_error=float(abs(A['short_internal_z'](oldz)[0,1,0]-shorted[1,0])),
        interpolation_vs_new_probe_max_Z_error_ohm=float(abs(oldz-z).max()),
        project_sha256=digest(p),snp_sha256=digest(p.with_suffix('.snp')))
    marker.write_text(json.dumps(report,indent=2)+'\n')
    (out/'completed.json').write_text(json.dumps(dict(project_sha256=digest(p),snp_sha256=digest(p.with_suffix('.snp')),
        frequencies_GHz=(n.f/1e9).tolist(),held_out_GHz=5.2),indent=2)+'\n')
    print(json.dumps(report,indent=2),flush=True);analyze(name);return report


if __name__=='__main__':run()
