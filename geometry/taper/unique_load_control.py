"""Validate uniquely named, ports-only design controls before expanding EM.

These files have unchanged planar metal and are never a final GDS design.
All-open and single-load controls must reproduce independently solved native
two-port projects. Three controls explore geometry directions, not extracted
production fit parameters.
"""
import json
import xml.etree.ElementTree as ET
import numpy as np
from circuit_common import read_snp
from geometry.prepare_sonnet import tag,write_xml
from geometry.taper.compact_stages import ROOT,NATIVE
from geometry.taper.response_cases import W
from geometry.taper.ideal_load_search import terminate,optimize
from geometry.taper.run_bridged_stage import digest


def check(*,two=False):
    word='two' if two else 'three'
    folder=ROOT/f'sonnet_files/taper_exploration/ideal_load_diagnostics/plain_{word}_unique_controls'
    p=folder/f'diagnostic_{word}_unique.sonx'
    out=ROOT/f'output_data/taper_exploration/ideal_load_diagnostics/plain_{word}_unique_controls'
    out.mkdir(parents=True,exist_ok=True)
    count=2 if two else 3
    n=W['read_network'](p.with_suffix('.snp'),2+2*count)
    references=[(NATIVE/'plain_w500_wire52_idc30_span172_l504/idc/generated_ncell_3.snp',[0]*count),
      (ROOT/'sonnet_files/taper_exploration/ideal_shunt_diagnostics/plain_w500_wire52_idc30_span172_l504_c420/diagnostic_ideal_shunt.snp',[420]+[0]*(count-1))]
    rows=[]
    for ref,loads in references:
        r=read_snp(ref);ix=np.array([np.argmin(abs(r.f-f)) for f in n.f])
        if max(abs(r.f[ix]-n.f))>1:raise ValueError('Control frequency absent from reference')
        e=abs(terminate(n.s,n.f,loads)-r.s[ix])
        rows.append(dict(known_caps_fF=loads,reference=str(ref),reference_sha256=digest(ref),
                         max_all_S_error=float(e.max()),max_S21_error=float(e[:,1,0].max())))
    report=dict(project_sha256=digest(p),operator_sha256=digest(p.with_suffix('.snp')),
                frequency_GHz=(n.f/1e9).tolist(),checks=rows,
                passed=all(q['max_all_S_error']<.002 and q['max_S21_error']<.001 for q in rows))
    (out/'control_checks.json').write_text(json.dumps(report,indent=2)+'\n')
    summary={k:v for k,v in report.items() if k!='frequency_GHz'}
    summary.update(sample_count=len(n.f),frequency_range_GHz=[float(n.f[0]/1e9),float(n.f[-1]/1e9)])
    print(json.dumps(summary,indent=2),flush=True)
    return p,out,report


def run(*,two=False):
    p,out,report=check(two=two)
    if not report['passed']:raise ValueError('Added ports perturb the native result; do not expand')
    for ext in ('.sonx','.snp'):(out/('two_point_control'+ext)).write_bytes(p.with_suffix(ext).read_bytes())
    tree=ET.parse(p);fs=tree.find(f"{tag('Sweeps')}/{tag('Set')}/{tag('Frequencies')}")
    for q in list(fs):fs.remove(q)
    ET.SubElement(fs,tag('Sweep'),Type='ABS_ENTRY',On='TRUE',Adaptive='TRUE',Start='4.65',Stop='5.4',Target='7500')
    for q in tree.iter(tag('File')):
        if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','TRUE')
    write_xml(tree,p);n=W['em'](p,6 if two else 8)
    _,_,full=check(two=two)
    if not full['passed']:raise ValueError('Full-band added-port reconstruction failed')
    return optimize(n,out)


def exploratory():
    """Uncertified search direction only; retain the failed validation gate.

    Added co-calibrated ports perturb the original narrow resonances. This
    operator cannot certify component values or accept a physical design.
    Its sole use is choosing a candidate for an ordinary two-port EM check.
    """
    p,out,report=check()
    if len(report['frequency_GHz'])>2:
        raise ValueError('Full operator already exists; optimize its saved SNP without rerunning EM')
    (out/'exploratory_intent.json').write_text(json.dumps(dict(
        diagnostic_only=True,model_verified=False,physical_design_accepted=False,
        original_validation=report,
        scope='Uncertified direction search; ordinary unprobed physical two-port validation required'),indent=2)+'\n')
    for ext in ('.sonx','.snp'):(out/('two_point_control'+ext)).write_bytes(p.with_suffix(ext).read_bytes())
    tree=ET.parse(p);fs=tree.find(f"{tag('Sweeps')}/{tag('Set')}/{tag('Frequencies')}")
    for q in list(fs):fs.remove(q)
    ET.SubElement(fs,tag('Sweep'),Type='ABS_ENTRY',On='TRUE',Adaptive='TRUE',Start='4.65',Stop='5.4',Target='7500')
    for q in tree.iter(tag('File')):
        if q.get('Type')=='TOUCH':q.set('UseInterpolatedData','TRUE')
    write_xml(tree,p);n=W['em'](p,8)
    check()  # Preserve the actual full-band error, including a failed gate.
    return optimize(n,out/'uncertified_direction_search')


if __name__=='__main__':run()
