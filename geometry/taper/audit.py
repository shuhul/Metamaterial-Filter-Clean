"""Verify a completed taper trial against its saved geometry and EM evidence."""
from pathlib import Path
from hashlib import sha256
import json
import xml.etree.ElementTree as ET
import gdstk
import numpy as np
from scipy.stats import linregress
from geometry.prepare_sonnet import tag,points,SOURCE
from geometry.taper.layout import validate
from circuit_common import read_snp,single_sweep_resonance
ROOT=Path(__file__).resolve().parents[2]


def signature(e):
    return e.tag,e.attrib,(e.text or '').strip(),[signature(c) for c in e]


def check_geometry(project,gds,metadata,*,geometry_validator=validate):
    _,original=geometry_validator(gds,metadata)
    expected=[]
    height=metadata['parameters']['l_sim'];top=metadata['parameters']['l_top_ground']
    for p in original:
        xy=p.points.copy();xy[:,1]+=height-top;xy=np.floor(xy/2+.5)*2;xy[:,1]=height-xy[:,1]
        expected.append(gdstk.Polygon(xy))
    r=ET.parse(project).getroot();actual=[gdstk.Polygon(points(p)) for p in r.iter(tag('PlanarPolygon'))]
    if gdstk.boolean(expected,actual,'xor',precision=.001):raise ValueError('Project metal differs from source GDS')
    if any(not np.allclose(p.points/2,np.round(p.points/2)) for p in actual):raise ValueError('Project is off grid')
    base=ET.parse(SOURCE).getroot()
    for key in ('Units','Control'):
        if signature(r.find(tag(key)))!=signature(base.find(tag(key))):raise ValueError(f'{key} changed')
    a=r.find(tag('Geometry'));b=base.find(tag('Geometry'))
    for key in ('Dielectric','Conductor','DielectricMaterialModel','BoxCover'):
        if [signature(x) for x in a.iter(tag(key))]!=[signature(x) for x in b.iter(tag(key))]:raise ValueError(f'Material/stack changed: {key}')
    size=a.find(f"{tag('Box')}/{tag('Size')}");mesh=a.find(f"{tag('Box')}/{tag('NumCells')}")
    if [float(size.get('X')),float(size.get('Y'))]!=[metadata['w_sim'],height]:raise ValueError('Box changed')
    if [float(mesh.get('X')),float(mesh.get('Y'))]!=[metadata['w_sim']/2,height/2]:raise ValueError('Mesh changed')
    if any(float(p.get('RefLength'))!=400 for p in a.iter(tag('RefPlane'))):raise ValueError('Reference planes changed')
    if len(list(a.iter(tag('Port'))))!=2:raise ValueError('Port count changed')


def audit(name):
    folder=ROOT/'output_data/taper_exploration'/name
    record=json.loads((folder/'trial.json').read_text());gds=Path(record['gds'])
    metadata=json.loads(gds.with_suffix('.json').read_text())
    if record['gds_sha256']!=sha256(gds.read_bytes()).hexdigest():raise ValueError('GDS changed')
    check_geometry(record['project'],gds,metadata)
    status=json.loads((folder/'run_status.json').read_text())
    if status['status']!='complete' or status['project_sha256']!=sha256(Path(record['project']).read_bytes()).hexdigest() or status['snp_sha256']!=sha256(Path(record['snp']).read_bytes()).hexdigest():raise ValueError('Mixed run identity failed')
    log=(folder/'em.log').read_text()
    if 'Em errors detected: 0    Em warnings detected: 0.' not in log or 'Em simulation completed' not in log:raise ValueError('Mixed run has warnings/errors or incomplete status')
    lp=Path(record['inductance_reference']);lr=json.loads(lp.read_text())
    if record['inductance_sha256']!=sha256(lp.read_bytes()).hexdigest():raise ValueError('Inductance record changed')
    iso_gds=ROOT/'geometry/taper/output'/f'{lp.parent.name}.gds';iso_meta=json.loads(iso_gds.with_suffix('.json').read_text())
    if iso_meta['taper']!=record['taper']:raise ValueError('Inductance geometry does not match taper')
    points_record=[]
    for entry in lr['sweep']:
        stem=lp.parent/f"cadd_{entry['C_add_fF']:g}ff"
        for suffix,key in (('.sonx','project_sha256'),('.snp','snp_sha256')):
            if sha256(stem.with_suffix(suffix).read_bytes()).hexdigest()!=entry[key]:raise ValueError('Archived inductance artifact changed')
        check_geometry(stem.with_suffix('.sonx'),iso_gds,iso_meta)
        log=stem.with_suffix('.log').read_text()
        if 'Em errors detected: 0    Em warnings detected: 0.' not in log:raise ValueError('Inductance run has warnings/errors')
        net=read_snp(stem.with_suffix('.snp'));f=single_sweep_resonance(net)
        if abs(f-entry['frequency_Hz'])>1:raise ValueError('Recorded resonance does not match Touchstone data')
        points_record.append((entry['C_add_fF']*1e-15,1/(2*np.pi*f)**2))
    regression=linregress(*np.array(points_record).T)
    if not np.isclose(regression.slope*1e9,lr['L0_nH'],rtol=1e-12):raise ValueError('Inductance regression mismatch')
    report={'trial':name,'gds_project_xor_empty':True,'mesh_um':[2,2],'reference_planes_um':[400,400],
            'materials_and_stack_preserved':True,'mixed_project_and_snp_hashes_match':True,
            'inductance_geometry_matches':True,'inductance_archive_hashes_match':True,'inductance_points_rechecked':len(points_record),
            'inductance_regression_recomputed':True,'all_native_errors_and_warnings':0,
            'code_sha256':{str(p.relative_to(ROOT)):sha256(p.read_bytes()).hexdigest() for p in [
                ROOT/'geometry/taper/layout.py',ROOT/'geometry/taper/fitting.py',ROOT/'geometry/taper/experiment.py',ROOT/'geometry/taper/sonnet.py',ROOT/'geometry/prepare_sonnet.py',ROOT/'extract_generated_inductance.py']}}
    (folder/'audit.json').write_text(json.dumps(report,indent=2)+'\n');return report
