"""Audit native ideal-C loading and its port perturbation before using it."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
os.environ['OMP_NUM_THREADS']='1'
import json
import xml.etree.ElementTree as ET
import numpy as np
from geometry.taper.compact_stages import NATIVE,OUT
from geometry.taper.response_cases import W
from geometry.taper.solver_guard import run as guarded
from geometry.taper.solver_resources import solver_environment
from geometry.prepare_sonnet import tag,write_xml,BIN
from geometry.taper.run_bridged_stage import digest
from circuit_common import read_snp


def export_operator(name):
    project=NATIVE/name/'direct/generated_ncell_3.sonx'
    out=OUT/name/'component_port_control';out.mkdir(parents=True,exist_ok=True)
    marker=out/'operator_completed.json';saved=out/'operator.snp'
    if marker.exists():return W['read_network'](saved,4),out
    original=project.read_bytes();original_s=project.with_suffix('.snp').read_bytes()
    try:
        tree=ET.parse(project);comp=next(q for q in tree.iter(tag('Component')) if q.get('MacroID')=='compactKnownCg')
        comp.set('Type','PORTS');comp.attrib.pop('Units',None);comp.remove(comp.find(tag('Value')))
        write_xml(tree,project);(out/'operator.sonx').write_bytes(project.read_bytes())
        with (out/'cache_export.log').open('w') as log:
            result=guarded([str(BIN/'em.exe'),'-v','-AbsNoDiscrete',str(project)],cwd=project.parent,stdout=log,env=solver_environment())
        result.check_returncode()
        text=(out/'cache_export.log').read_text()
        if 'Total time for 0 frequencies:' not in text:raise RuntimeError('Cache-only export did not confirm zero EM solves')
        saved.write_bytes(project.with_suffix('.snp').read_bytes())
        n=W['read_network'](saved,4)
    finally:
        project.write_bytes(original);project.with_suffix('.snp').write_bytes(original_s)
    actual=read_snp(project.with_suffix('.snp'))
    cap=W['terminate_capacitor'](n.s,n.f,.0567)
    if not np.array_equal(n.f,actual.f):raise ValueError('Operator grid differs from original')
    error=float(np.max(abs(cap-actual.s)))
    swap=np.array([[0,1],[1,0]],complex)[None,:,:];eye=np.eye(2)[None,:,:]
    short=n.s[:,:2,:2]+n.s[:,:2,2:]@swap@np.linalg.solve(eye-n.s[:,2:,2:]@swap,n.s[:,2:,:2])
    np.savez(out/'short_prediction.npz',f=n.f,s=short)
    marker.write_text(json.dumps(dict(source_project_sha256=digest(project),source_snp_sha256=digest(project.with_suffix('.snp')),
        operator_sha256=digest(saved),known_Cg_fF=56.7,native_cap_termination_max_S_error=error,
        zero_additional_EM_frequencies=True,short_nonperturbation_pending=True),indent=2)+'\n')
    print('Native known-C contraction max error:',error,flush=True)
    return n,out


def compare_bridge(name,physical_snp):
    _,out=export_operator(name);pred=np.load(out/'short_prediction.npz');n=W['read_network'](physical_snp,2)
    ix=np.array([np.argmin(abs(pred['f']-f)) for f in n.f])
    if np.max(abs(pred['f'][ix]-n.f))>1:raise ValueError('No exact shared frequency')
    d=abs(pred['s'][ix]-n.s)
    result=dict(frequencies_GHz=(n.f/1e9).tolist(),max_all_S_error=float(d.max()),
                S21_complex_error=d[:,1,0].tolist(),passed=bool(d.max()<.005),
                interpretation='A passing comparison supports this control; a mismatch can include both port effects and the finite L/C of the added 100 um physical bridge. It does not uniquely identify port error.')
    (out/'physical_bridge_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)
    return result
