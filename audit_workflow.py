"""Read-only audit of saved GDS, Sonnet projects, S-parameters and regression."""

from hashlib import sha256
from pathlib import Path
import json
import re
import xml.etree.ElementTree as ET

import gdstk
import numpy as np
from scipy.stats import linregress

from geometry.validation import require, validate_gds_metadata
from geometry.prepare_sonnet import points, tag
from circuit_common import read_snp, single_sweep_resonance

BASE=Path(__file__).resolve().parent


def main():
    report={'projects':{},'networks':{},'limitations':[
        '2 um mesh convergence is not established; no finer-mesh runs were performed.',
        'ABS interpolation accuracy near narrow peaks needs spot checks before claiming peak loss.',
        'End-cell loading, reference-plane sensitivity, and transfer of N=1 L0 to larger N are not independently established.',
        'Default multi-cell fitting still uses target L0 unless --l0-nh is specified.',
    ]}
    baseline=BASE/'sonnet_files'/'filter_v6_cells_n=4.sonx'
    b=ET.parse(baseline)
    def sig(e):
        return (e.tag,e.attrib,(e.text or '').strip(),[sig(x) for x in e])
    for n in (1,4):
        folder=BASE/'sonnet_files'/f'generated_ncell_{n}'
        gds=BASE/'geometry'/'output'/f'ncell_{n}.gds'
        metadata=json.loads(gds.with_suffix('.json').read_text())
        p,_=validate_gds_metadata(gds,metadata)
        tree=ET.parse(folder/f'generated_ncell_{n}.sonx'); root=tree.getroot()
        polygons=tree.findall(f'.//{tag("PlanarPolygon")}')
        actual=[gdstk.Polygon(points(q)) for q in polygons]
        snapped=gdstk.read_gds(str(folder/f'ncell_{n}_grid2.gds')).top_level()[0].get_polygons()
        target=[gdstk.Polygon([(x,p.l_sim-y) for x,y in q.points]) for q in snapped]
        require(not gdstk.boolean(actual,target,'xor',precision=.001),'GDS/Sonnet mismatch')
        require(all(np.allclose(points(q)/2,np.round(points(q)/2),rtol=0,atol=1e-9) for q in polygons),'Off-grid vertex')
        require(len(gdstk.boolean(actual,[],'or'))==1,'Disconnected ground/inductor network')
        ports=tree.findall(f'.//{tag("PlanarPolygon")}/{tag("Objects")}/{tag("Port")}')
        require(len(ports)==2,'Wrong external port count')
        for poly in polygons:
            xy=points(poly)
            for port in poly.findall(f'{tag("Objects")}/{tag("Port")}'):
                i=int(port.get('Vertex')); a,bb=xy[i],xy[(i+1)%len(xy)]
                require(a[0]==bb[0] and a[0] in (0,p.w_sim),'Port is not on a sidewall')
                require(sorted((a[1],bb[1]))==[456,466],'Port is not on the feed edge')
        for key in ('Control','Sweeps','Units'):
            require(sig(root.find(tag(key)))==sig(b.getroot().find(tag(key))),f'{key} differs from baseline')
        for key in ('Dielectric','Conductor','DielectricMaterialModel','BoxCover'):
            require([sig(e) for e in tree.findall(f'.//{tag(key)}')]==[sig(e) for e in b.findall(f'.//{tag(key)}')],f'{key} changed')
        for layer in tree.findall(f'.//{tag("PlanarTechLayer")}'):
            if layer.findall(f'.//{tag("PlanarPolygon")}'):
                require(layer.get('MaterialName')=='SUPERCONDUCTOR','Wrong metal')
                for key in ('PlanarMaterialModel','SubsectionMeshing','PlanarMeshing'):
                    require(sig(layer.find(tag(key)))==sig(b.find(f'.//{tag(key)}')),f'Metal/mesh setting {key} differs')
        size=tree.find(f'.//{tag("Box")}/{tag("Size")}'); mesh=tree.find(f'.//{tag("NumCells")}')
        require(float(size.get('X'))==p.w_sim and float(size.get('Y'))==1500,'Wrong box')
        require(float(size.get('X'))/float(mesh.get('X'))==2 and 1500/float(mesh.get('Y'))==2,'Wrong mesh')
        require(all(float(q.get('RefLength'))==400 for q in tree.findall(f'.//{tag("RefPlane")}')),'Wrong reference planes')
        saved=json.loads((folder/'validation.json').read_text())
        require(saved['baseline_sha256']==sha256(baseline.read_bytes()).hexdigest(),'Baseline changed')
        report['projects'][str(n)]={'box_um':[p.w_sim,p.l_sim],'polygons':len(polygons),'ports':2,'mesh_um':2,'geometry_xor_empty':True}

    folder=BASE/'output_data'/'generated_inductance'
    all_files=[baseline.with_suffix('.snp'),BASE/'sonnet_files'/'generated_ncell_4'/'generated_ncell_4.snp',*sorted(folder.glob('cadd_*.snp'))]
    for path in all_files:
        net=read_snp(path)
        singular=np.linalg.svd(net.s,compute_uv=False)[:,0]
        reciprocal=float(np.max(abs(net.s[:,0,1]-net.s[:,1,0])))
        require(np.max(singular)<=1+1e-6,'Passivity violation')
        require(reciprocal<1e-6,'Reciprocity violation')
        report['networks'][str(path.relative_to(BASE))]={'samples':len(net.f),'max_singular_value':float(max(singular)),
                                                        'reciprocity_max_error':reciprocal,'sha256':sha256(path.read_bytes()).hexdigest()}
    frequencies=[]
    previous_geometry=None
    for c in (20,30,40,50,.1):
        name=f'cadd_{c:g}ff'; path=folder/f'{name}.snp'
        text=path.read_text(); value=float(re.search(r'^!< C_add = (\S+)',text,re.M)[1])*1000
        require(abs(value-c)<1e-9,'Stale capacitor export')
        require('Em errors detected: 0    Em warnings detected: 0.' in (folder/f'{name}.log').read_text(),'Solver errors/warnings')
        net=read_snp(path); freq=single_sweep_resonance(net)
        t=ET.parse(folder/f'{name}.sonx')
        geometry=[points(q).tolist() for q in t.findall(f'.//{tag("PlanarPolygon")}')]
        if previous_geometry is not None: require(geometry==previous_geometry,'Geometry changed during C_add sweep')
        previous_geometry=geometry
        comp,=t.findall(f'.//{tag("Component")}')
        require(comp.findtext(tag('Value'))=='C_add','Wrong component expression')
        require([q.findtext(tag('Center')) for q in comp.findall(tag('Pin'))]==['(940,448)','(940,456)'],'Wrong capacitor pins')
        if c!=.1: frequencies.append(freq)
    f=np.array(frequencies); x=np.array([20,30,40,50])*1e-15
    fit=linregress(x,1/(2*np.pi*f)**2)
    saved=json.loads((folder/'results.json').read_text())
    require(abs(fit.slope*1e9-saved['L0_nH'])<1e-10,'Inductance not reproducible')
    report['inductance']={'L0_nH':fit.slope*1e9,'R_squared':fit.rvalue**2,'frequencies_GHz':(f/1e9).tolist()}
    out=BASE/'output_data'/'workflow_audit'; out.mkdir(exist_ok=True)
    (out/'checks.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS: 2 generated projects, {len(all_files)} networks, 5 C_add snapshots, and inductance regression.')
    print(f'Results: {out / "checks.json"}')


if __name__=='__main__': main()
