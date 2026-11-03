"""Check full-geometry capacitor/inductor interaction without changing metal.

First verify that ideal shorts on zero-width series probes reproduce the
original two-port simulation. Only then solve three low frequencies for a
whole-conductor capacitance and impedance diagnostic. No final acceptance.
"""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json,runpy
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np
from circuit_common import read_snp
from geometry.taper.run_bridged_stage import OUT,digest
from geometry.taper.bridged_fixtures import ROOT
from geometry.prepare_sonnet import tag,write_xml

G=runpy.run_path(str(ROOT/'investigations/taper2_fit_2026-10-22/gap_probes.py'))
A=runpy.run_path(str(ROOT/'investigations/taper2_fit_2026-10-22/analyze_series.py'))
W=G['W']


def run(name):
    case='04_t2_idc_restored'
    source=ROOT/'sonnet_files/taper_exploration/bridged_trials'/name/case/'generated_ncell_3.sonx'
    meta=json.loads((ROOT/'geometry/taper/output/bridged_trials'/name/case/'geometry.json').read_text())
    out=OUT/'trials'/name/'inplace_wire_check';out.mkdir(parents=True,exist_ok=True)
    probe=G['prepare'](source,meta,out,first_frequency=5.4)
    manifest=json.loads((out/'manifest.json').read_text())
    if manifest['source_sha256']!=digest(source):raise ValueError('Source changed after probe preparation')
    marker=out/'control_result.json'
    if not marker.exists():
        n=W['em'](probe,6);actual=G['close_wires'](n.s)
        original=read_snp(source.with_suffix('.snp'))
        indices=[int(np.argmin(abs(original.f-f))) for f in n.f]
        if np.max(abs(original.f[indices]-n.f))>1:raise ValueError('Control has no exact frequency reference')
        error=abs(actual-original.s[indices]);relative=error[:,1,0]/np.maximum(abs(original.s[indices,1,0]),1e-8)
        control=dict(frequencies_GHz=(n.f/1e9).tolist(),max_all_S_error=float(error.max()),
            relative_S21_error=float(relative.max()),physical_union_unchanged=True,
            passed=bool(error.max()<.005 and relative.max()<.01))
        marker.write_text(json.dumps(control,indent=2)+'\n')
        print(json.dumps(control),flush=True)
        for ext in ('.snp','.sonx'):probe.with_name('control'+ext).write_bytes(probe.with_suffix(ext).read_bytes())
        (out/'control_em.log').write_bytes((out/'em.log').read_bytes())
        resource=out/'em.log.resources.json'
        if resource.exists():(out/'control_em.log.resources.json').write_bytes(resource.read_bytes())
    if not json.loads(marker.read_text())['passed']:raise ValueError('Probe nonperturbation failed; no expanded run')
    completed=out/'low_frequency_completed.json'
    if not completed.exists():
        tree=ET.parse(probe);freq=tree.find(f"{tag('Sweeps')}/{tag('Set')}/{tag('Frequencies')}")
        for q in list(freq):freq.remove(q)
        for f in (.5,1,2,5.4):
            ET.SubElement(freq,tag('Sweep'),Type='SWEEP',On='TRUE',Adaptive='FALSE',Start=str(f),Stop=str(f),Step='1')
        write_xml(tree,probe);n=W['em'](probe,6)
        completed.write_text(json.dumps(dict(project_sha256=digest(probe),snp_sha256=digest(probe.with_suffix('.snp'))),indent=2)+'\n')
    else:
        saved=json.loads(completed.read_text())
        if saved['project_sha256']!=digest(probe) or saved['snp_sha256']!=digest(probe.with_suffix('.snp')):raise ValueError('Probe archive changed')
        n=W['read_network'](probe.with_suffix('.snp'),6)
    polarity=np.diag([1,1,-1,-1,-1,-1]);z=A['s_to_z'](polarity@n.s@polarity)
    low=n.f<3e9;c,coeff=A['extract_capacitance'](n.f[low],z[low])
    # Coefficient of omega in terminal Z includes distributed conductor and
    # fixture contributions. Do not mislabel it as an independently isolated L.
    basis=(n.f[low,None]/1e9)**np.array([-1,1,3])[None,:]
    allcoeff=np.linalg.solve(basis,z[low].imag.reshape(3,-1))
    l=allcoeff[1].reshape(6,6)/(2*np.pi*1e9)*1e9
    result=dict(conductor_order=['Input','T2','T1','Regular','Output'],
        C_matrix_fF=c.tolist(),C_ground_fF=c.sum(axis=1).tolist(),
        Cg_input_T2_fF=float(-c[0,1]),Cg_T2_T1_fF=float(-c[1,2]),
        terminal_linear_Z_coefficient_nH=l.tolist(),accepted=False,
        limitation='Static conductor C and terminal impedance coefficients, not direct equivalents of the ideal ladder L/C values. Three low points provide no independent extrapolation check.')
    (out/'diagnostic.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez(out/'operator.npz',f=n.f,z=z,static_C_fF=c,linear_Z_nH=l)
    print(json.dumps(result,indent=2),flush=True)
    return result


if __name__=='__main__':run('w52_frame800_wireleft_right548')
