r"""Prepare the best measured T2 layout for coarse, 4 um mesh iterations.

Run from the project root: .venv\Scripts\python.exe -m geometry.taper.taper2_grid4
Default: generate and validate files only; do not run EM or fit old responses.
The 2 um source is immutable. Both future mesh comparisons should use the
same snapped geometry, so geometry changes are separated from mesh effects.
"""
from pathlib import Path
from hashlib import sha256
import json
import subprocess
import xml.etree.ElementTree as ET
import numpy as np
import gdstk
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from geometry.prepare_sonnet import tag, points, write_xml, BIN
from geometry.taper.audit import signature

ROOT = Path(__file__).resolve().parents[2]
SOURCE_CASE = 'idc_rail8_right30_dl256_pad148_w148_p200_i42n70_j50_r74tg2'
SOURCE = ROOT/'sonnet_files/taper_exploration/three_tapers_v3_T2_plain_2026-10-28'/SOURCE_CASE/'response/generated_ncell_3.sonx'
OLD_FIT = ROOT/'output_data/taper_exploration/three_tapers_v3_T2_plain_2026-10-28'/SOURCE_CASE/'response/extraction/fit.json'
DEST = ROOT/'sonnet_files/taper_exploration/three_tapers_v3_T2_grid4_2026-10-29/best_003853_grid4/response'
RUN_SIMULATION = False
RUN_FIT = False
GRID = 4
OLDER_REFERENCE = ROOT/'geometry/taper/taper2_fixed_older_reference.json'

# Choose the downward tie where nearest rounding would close an IDC tip gap
# or erase the thin ground strip. The T1 right edge also rounds downward so
# the widened strip remains exactly centered between T1 and the regular cap.
X_TIES = {982:980, 990:988, 1030:1028, 1990:1988, 2026:2024}


def fit_policy():
    """Use the user-adopted 1% older-circuit values, then freeze them."""
    ref=json.loads(OLDER_REFERENCE.read_text())
    source=Path(ref['source_fit'])
    if sha256(source.read_bytes()).hexdigest()!=ref['source_fit_sha256']:
        raise ValueError('Adopted older-circuit fit changed')
    return dict(free_capacitances=ref['free_capacitances'],
        names=['T2 C0','T2 input Cg','T2-to-T1 Cg'],fixed_L_nH=ref['fixed_L_nH'],
        initial_C_fF=[ref['T2_seed_C_fF'],ref['fixed_C_fF']['T1'],ref['fixed_C_fF']['regular'],ref['fixed_C_fF']['output_end']],
        initial_G_fF=ref['T2_seed_G_fF']+[ref['fixed_G_fF']['T1_to_regular'],ref['fixed_G_fF']['regular_to_output']],
        fixture_delay_free=True,delay_bounds_ps=[-100,100],metric='complex S21; real and imaginary residuals',
        fixed_older_values_source=str(OLDER_REFERENCE),
        fixed_older_values_sha256=sha256(OLDER_REFERENCE.read_bytes()).hexdigest(),
        inductance_note='2um-source calibrations provisionally fixed for coarse screening; not remeasured on snapped metal.',
        target_RMS_cutoff=.025)


def snap(xy):
    xy=np.asarray(xy,float);out=np.floor(xy/GRID+.5)*GRID
    for old,new in X_TIES.items():out[xy[:,0]==old,0]=new
    return out


def clean(xy):
    # Keep the port edge ordering; remove only zero-length consecutive edges.
    return xy[np.any(xy!=np.roll(xy,1,axis=0),axis=1)]


def contacts(polys):
    """Polygon incidence, including shared edges/points, independent of mesh."""
    expanded=[gdstk.offset([p],.001,precision=.0001) for p in polys]
    boxes=[np.array(p.bounding_box()) for p in polys];pairs=[]
    for i,a in enumerate(polys):
        for j in range(i):
            if np.any(boxes[i][0]>boxes[j][1]+.002) or np.any(boxes[j][0]>boxes[i][1]+.002):continue
            if gdstk.boolean(expanded[i],[polys[j]],'and',precision=.0001):pairs.append((j,i))
    return set(pairs)


def prepare():
    DEST.mkdir(parents=True,exist_ok=True)
    original_bytes=SOURCE.read_bytes();source_hash=sha256(original_bytes).hexdigest()
    tree=ET.parse(SOURCE);root=tree.getroot();geometry=root.find(tag('Geometry'))
    elements=list(geometry.iter(tag('PlanarPolygon')))
    original=[gdstk.Polygon(points(p)) for p in elements]
    # This tie map is scoped to this exact best-run geometry, not a general
    # arbitrary-GDS snapping routine.
    if original[1].bounding_box()!=((2026.,448.),(2028.,1066.)):
        raise ValueError('Source geometry changed; review the scoped tie choices')
    transformed=[];moves=[];port_count=0
    for element,old in zip(elements,original):
        raw=snap(old.points);xy=clean(raw);new=gdstk.Polygon(xy)
        if new.area()<=0 or len(gdstk.boolean([new],[],'or',precision=.001))!=1:
            raise ValueError('Snapping collapsed or split a polygon')
        moves.append(float(np.max(abs(raw-old.points))))
        # Reattach ports by their transformed edge, not an obsolete vertex index.
        for port in element.findall(f'{tag("Objects")}/{tag("Port")}'):
            edge=int(port.get('Vertex'));aa=raw[edge];bb=raw[(edge+1)%len(raw)]
            matches=[i for i,(a,b) in enumerate(zip(xy,np.roll(xy,-1,axis=0)))
                     if np.array_equal(a,aa) and np.array_equal(b,bb)]
            if len(matches)!=1 or np.array_equal(aa,bb):raise ValueError('Port edge collapsed')
            port.set('Vertex',str(matches[0]));center=(aa+bb)/2
            port.find(tag('Center')).text=f'({center[0]:g},{center[1]:g})';port_count+=1
        oldnodes=element.findall(tag('Points'))
        for q in oldnodes:element.remove(q)
        ET.SubElement(element,tag('Points')).text=' '.join(f'({x:g},{y:g})' for x,y in xy)
        transformed.append(new)
    if port_count!=2:raise ValueError('Expected exactly two physical box ports')
    if max(moves)>2.000001:raise ValueError('Metal moved farther than nearest 4um grid point')
    old_components=len(gdstk.boolean(original,[],'or',precision=.001))
    new_components=len(gdstk.boolean(transformed,[],'or',precision=.001))
    old_contacts=contacts(original);new_contacts=contacts(transformed)
    if old_components!=new_components or old_contacts!=new_contacts:
        raise ValueError(f'Connectivity changed: components {old_components}->{new_components}; '
                         f'new contacts {new_contacts-old_contacts}; lost contacts {old_contacts-new_contacts}')
    size=geometry.find(f'{tag("Box")}/{tag("Size")}')
    old_size=[float(size.get(a)) for a in ('X','Y')]
    width,height=(np.floor(np.array(old_size)/4+.5)*4).astype(int)
    size.set('X',str(width));size.set('Y',str(height))
    mesh=geometry.find(f'{tag("Box")}/{tag("NumCells")}')
    mesh.set('X',str(width//4));mesh.set('Y',str(height//4))
    project=DEST/'generated_ncell_3.sonx'
    for q in root.iter(tag('File')):
        if q.get('Type')=='TOUCH':q.set('Name',project.with_suffix('.snp').name)
    root.find(tag('OutputFiles')).set('Folder','.')
    if list(geometry.iter(tag('Component'))):raise ValueError('Ideal component in production layout')
    write_xml(tree,project)
    # Physical-coordinate GDS of precisely the native metal. Y points upward.
    lib=gdstk.Library(unit=1e-6,precision=1e-9);top=lib.new_cell('TOP')
    for p in transformed:
        xy=p.points.copy();xy[:,1]=height-xy[:,1]
        top.add(gdstk.Polygon(xy,layer=1,datatype=0))
    gds=DEST/'geometry_grid4.gds';lib.write_gds(str(gds))
    loaded=gdstk.read_gds(str(gds)).top_level()[0].get_polygons()
    back=[]
    for p in loaded:
        xy=p.points.copy();xy[:,1]=height-xy[:,1];back.append(gdstk.Polygon(xy))
    if gdstk.boolean(back,transformed,'xor',precision=.001):raise ValueError('GDS/native XOR failed')
    if not all(np.allclose(p.points/4,np.round(p.points/4)) for p in transformed):raise ValueError('Metal off grid')
    baseline=ET.fromstring(original_bytes);oldgeom=baseline.find(tag('Geometry'))
    for name in ('Units','Control','Sweeps'):
        if signature(root.find(tag(name)))!=signature(baseline.find(tag(name))):raise ValueError(name+' changed')
    for name in ('Dielectric','Conductor','DielectricMaterialModel','BoxCover','RefPlane'):
        if [signature(x) for x in geometry.iter(tag(name))]!=[signature(x) for x in oldgeom.iter(tag(name))]:
            raise ValueError(name+' changed')
    read=subprocess.run([str(BIN/'soncmd.exe'),'-Read',str(project)],cwd=DEST,capture_output=True,text=True)
    (DEST/'sonnet_read.log').write_text(read.stdout+read.stderr);read.check_returncode()
    policy=fit_policy()
    (DEST/'fit_policy.json').write_text(json.dumps(policy,indent=2)+'\n')
    report=dict(source=str(SOURCE),source_sha256=source_hash,source_target_RMS=.03853008915822824,
        project=str(project),gds=str(gds),mesh_um=[4,4],box_um=[int(width),int(height)],
        old_box_um=old_size,mesh_cells=[int(width//4),int(height//4)],
        max_vertex_displacement_per_axis_um=max(moves),x_tie_choices_um=X_TIES,
        metal_components_before=old_components,metal_components_after=new_components,
        polygon_contact_graph_preserved=True,polygon_count=len(elements),
        input_IDC_tip_gaps_um=4,centered_ground_strip_width_um=4,
        centered_ground_strip_center_x_um=2026,facing_cap_edges_x_um=[1988,2064],
        source_unchanged=SOURCE.read_bytes()==original_bytes,source_T2_pads_retained=True,
        materials_stack_sweeps_control_and_reference_planes_preserved=True,
        GDS_native_XOR_empty=True,Sonnet_read_exit_code=read.returncode,
        simulation_run=False,coarse_design_accepted=False,
        final_verification='Use identical snapped metal on a2x2um mesh; do not confuse mesh error with geometry change.')
    if not report['source_unchanged']:raise ValueError('Source was modified')
    for path in (project,gds,DEST/'fit_policy.json'):
        report.setdefault('files_sha256',{})[path.name]=sha256(path.read_bytes()).hexdigest()
    (DEST/'validation.json').write_text(json.dumps(report,indent=2)+'\n')
    fig,ax=plt.subplots(figsize=(17,6))
    for p in transformed:ax.add_patch(Polygon(p.points,fc='#287b53',ec='none'))
    ax.set(xlim=(360,width-360),ylim=(1100,420),aspect='equal',xlabel='Sonnet x (um)',ylabel='Sonnet y (um)',
           title='Best T2 Layout: Snapped to 4 x 4 um Grid')
    fig.tight_layout();fig.savefig(DEST/'geometry.png',dpi=180);plt.close(fig)
    fig,axes=plt.subplots(1,2,figsize=(12,7))
    for ax,limits,title in zip(axes,[(970,1050,1050,450),(1970,2080,1090,430)],['Input IDC: 4 um Tip Gaps','Centered 4 um Ground Strip']):
        for p in transformed:ax.add_patch(Polygon(p.points,fc='#287b53',ec='none'))
        ax.set(xlim=limits[:2],ylim=limits[2:],aspect='equal',title=title)
        ax.tick_params(labelsize=10)
    fig.tight_layout();fig.savefig(DEST/'details.png',dpi=180);plt.close(fig)
    print(json.dumps(report,indent=2),flush=True)
    return project,policy


def main():
    project,policy=prepare()
    if RUN_SIMULATION:
        from geometry.taper.sonnet import simulate
        simulate(project,DEST/'simulation')
    if RUN_FIT:
        from geometry.taper.progressive_fitting import fit
        if not project.with_suffix('.snp').exists():raise ValueError('Run the new 4um response before fitting')
        fit(project.with_suffix('.snp'),DEST/'three_cap_fit',policy['fixed_L_nH'],policy['initial_C_fF'],
            policy['initial_G_fF'],[tuple(q) for q in policy['free_capacitances']],title='T2: Three-Capacitance Fit')


if __name__=='__main__':main()
