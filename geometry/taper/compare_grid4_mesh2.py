"""Compare 2um and 4um meshes on exactly the same snapped T2 metal.

Only Box/NumCells changes. No geometry snapping, fitting, or cache copying.
"""
import os
os.environ['SONNET_MAX_JOBS']='1'
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
from hashlib import sha256
from copy import deepcopy
import json
import xml.etree.ElementTree as ET
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.signal import find_peaks,peak_widths
from geometry.prepare_sonnet import tag,write_xml
from geometry.taper.audit import signature
from geometry.taper import taper2_grid4 as layout
from geometry.taper.run_taper2_grid4 import OUT as COARSE_OUT
from geometry.taper.sonnet import simulate
from circuit_common import read_snp

NATIVE=layout.DEST.parent.parent/'best_003853_grid4_verify_mesh2/response'
OUT=COARSE_OUT.parent.parent/'best_003853_grid4_verify_mesh2/response'


def main():
    NATIVE.mkdir(parents=True,exist_ok=True);OUT.mkdir(parents=True,exist_ok=True)
    source=layout.DEST/'generated_ncell_3.sonx';project=NATIVE/source.name
    coarse_status=json.loads((COARSE_OUT/'run_status.json').read_text())
    if sha256(source.read_bytes()).hexdigest()!=coarse_status['project_sha256']:
        raise ValueError('Coarse source changed')
    if sha256(source.with_suffix('.snp').read_bytes()).hexdigest()!=coarse_status['snp_sha256']:
        raise ValueError('Coarse response changed')
    tree=ET.parse(source);root=tree.getroot();mesh=root.find(f'{tag("Geometry")}/{tag("Box")}/{tag("NumCells")}')
    old_attributes=mesh.attrib.copy()
    mesh.set('X',str(int(float(mesh.get('X'))*2)));mesh.set('Y',str(int(float(mesh.get('Y'))*2)))
    check=deepcopy(root);check.find(f'{tag("Geometry")}/{tag("Box")}/{tag("NumCells")}').attrib=old_attributes
    if signature(check)!=signature(ET.parse(source).getroot()):raise ValueError('Something besides mesh changed')
    if not project.exists():write_xml(tree,project)
    if signature(ET.parse(project).getroot())!=signature(root):raise ValueError('Prepared fine project differs')
    manifest=OUT/'mesh_only_proof.json'
    if not manifest.exists():manifest.write_text(json.dumps(dict(source=str(source),project=str(project),
        source_sha256=sha256(source.read_bytes()).hexdigest(),project_sha256=sha256(project.read_bytes()).hexdigest(),
        identical_metal_ports_materials_stack_controls_sweeps_output_settings=True,
        only_XML_change='Geometry/Box/NumCells X788->1576;Y375->750',
        fresh_fine_EM_cache=not (NATIVE/'sondata').exists(),
        raw_S21_comparison=True,phase_alignment_or_circuit_fitting=False),indent=2)+'\n')
    status=OUT/'run_status.json'
    if status.exists() and json.loads(status.read_text()).get('status')=='complete':
        state=json.loads(status.read_text())
        if sha256(project.with_suffix('.snp').read_bytes()).hexdigest()!=state['snp_sha256']:raise ValueError('Fine response changed')
        snp=project.with_suffix('.snp')
    else:snp=simulate(project,OUT)
    a=read_snp(source.with_suffix('.snp'));b=read_snp(snp)
    if len(a.f)!=len(b.f) or not np.allclose(a.f,b.f,atol=1,rtol=0):raise ValueError('Different frequency samples')
    f=a.f;x=a.s[:,1,0];y=b.s[:,1,0];error=abs(x-y)
    coarse_resource=json.loads((COARSE_OUT/'em.log.resources.json').read_text())
    fine_resource=json.loads((OUT/'em.log.resources.json').read_text())
    peaks=[]
    for name,s in [('mesh4',x),('mesh2',y)]:
        ix=find_peaks(abs(s),prominence=.05)[0]
        widths=peak_widths(abs(s),ix,rel_height=.5)[0]*(f[1]-f[0])/1e6
        peaks.append([dict(frequency_GHz=float(f[i]/1e9),magnitude=float(abs(s[i])),
                           half_prominence_width_MHz=float(w)) for i,w in zip(ix,widths)])
    result=dict(complex_S21_RMS=float(np.sqrt(np.mean(error**2))),complex_S21_max=float(error.max()),
        magnitude_RMS=float(np.sqrt(np.mean((abs(x)-abs(y))**2))),
        magnitude_max=float(np.max(abs(abs(x)-abs(y)))),
        coarse_solver_seconds=coarse_resource['elapsed_seconds'],fine_solver_seconds=fine_resource['elapsed_seconds'],
        same_geometry_speedup=fine_resource['elapsed_seconds']/coarse_resource['elapsed_seconds'],
        coarse_peaks=peaks[0],fine_peaks=peaks[1],
        peak_frequency_shifts_MHz=[1000*(p['frequency_GHz']-q['frequency_GHz']) for p,q in zip(*peaks)] if len(peaks[0])==len(peaks[1]) else None,
        samples=len(f),frequency_range_GHz=[float(f[0]/1e9),float(f[-1]/1e9)],
        phase_alignment_or_fitting=False,proof=str(manifest),
        limitation='Two mesh levels establish sensitivity, not proof that2um is fully converged.')
    (OUT/'mesh_comparison.json').write_text(json.dumps(result,indent=2)+'\n')
    np.savez(OUT/'mesh_comparison.npz',f=f,mesh4=x,mesh2=y)
    fig,(ax,err)=plt.subplots(2,1,figsize=(14,10),sharex=True,gridspec_kw={'height_ratios':[3,1]})
    ax.plot(f/1e9,abs(x),lw=2.6,color='#2768aa',label='4 x 4 um mesh')
    ax.plot(f/1e9,abs(y),'--',lw=2.3,color='#b72d25',label='2 x 2 um mesh; identical metal')
    ax.set_title('T2 Mesh Comparison: Identical Snapped Geometry',fontsize=23,pad=12)
    ax.set_ylabel(r'$|S_{21}|$',fontsize=23);ax.legend(fontsize=14,loc='upper right')
    text=f'Raw complex S21 RMS difference: {result["complex_S21_RMS"]:.5f}\nMagnitude RMS difference: {result["magnitude_RMS"]:.5f}\nNo phase alignment or circuit fit'
    ax.text(.025,.96,text,transform=ax.transAxes,va='top',fontsize=15,bbox=dict(boxstyle='round,pad=.5',fc='white',ec='gray'))
    ax.set_ylim(0,max(abs(x).max(),abs(y).max())*1.5)
    err.plot(f/1e9,error,color='#594477',lw=2);err.set_ylabel(r'$|\Delta S_{21}|$',fontsize=19)
    err.set_xlabel('Frequency (GHz)',fontsize=20)
    for q in (ax,err):q.grid(alpha=.2);q.tick_params(labelsize=15)
    fig.tight_layout()
    for ext in ('png','pdf'):fig.savefig(OUT/f'mesh_comparison.{ext}',dpi=190)
    plt.close(fig);print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':main()
