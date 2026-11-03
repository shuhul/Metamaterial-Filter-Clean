"""Independent checks of the new series-probe extraction and circuit reduction."""
import os
os.environ['OPENBLAS_NUM_THREADS']='1'
from pathlib import Path
import runpy
import unittest
import json
import tempfile
import xml.etree.ElementTree as ET
import numpy as np

A=runpy.run_path(str(Path(__file__).resolve().parents[1]/'investigations/taper2_fit_2026-10-22/analyze_series.py'))


def known_cap(c0=100,cg=459):
    c=np.diag([84.,c0,101.49,188.45,190.])
    for i,g in enumerate([cg,56.15,8.229,8.312]):
        c[i,i]+=g;c[i+1,i+1]+=g;c[i,i+1]-=g;c[i+1,i]-=g
    return c


def known_z(f,c):
    mapping=A['conductor_mapping']()
    ell=np.diag([0.,0.,.74,6.014,4.923,4.923])
    ell[2,3]=ell[3,2]=.24;ell[2,4]=ell[4,2]=.06
    ell[3,4]=ell[4,3]=.10
    w=2*np.pi*f
    return (mapping.T@np.linalg.inv(c*1e-15)@mapping)[None,:,:]/(1j*w[:,None,None])+1j*w[:,None,None]*ell*1e-9


class SeriesProbeTests(unittest.TestCase):
    def test_series_probe_tracks_offset_wire_and_preserves_metal(self):
        g=runpy.run_path(str(Path(__file__).resolve().parents[1]/'investigations/taper2_fit_2026-10-22/gap_probes.py'))
        tag=g['tag'];fn=g['prepare'];env=fn.__globals__
        saved={k:env[k] for k in ('ROOT','OUT','CASE')}
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);name='offset_test'
            source=root/'sonnet_files/taper_exploration'/name/'generated_ncell_3.sonx'
            source.parent.mkdir(parents=True)
            meta_path=root/'geometry/taper/output'/f'{name}.json';meta_path.parent.mkdir(parents=True)
            meta_path.write_text(json.dumps({'cell_centers':[0,20,60,100,140],'w_sim':200,
                'second':{'x_meander':-4,'l_ground_rise':0}}))
            r=ET.Element(tag('Project'));layer=ET.SubElement(r,tag('PlanarTechLayer'))
            objects=ET.SubElement(layer,tag('Objects'))
            for i,x in enumerate((16,60,100,140)):
                q=ET.SubElement(objects,tag('PlanarPolygon'),Id=str(i+1))
                ET.SubElement(q,tag('Points')).text=f'({x-2},400) ({x+2},400) ({x+2},1100) ({x-2},1100)'
            ET.SubElement(r,tag('Sweeps'));ET.ElementTree(r).write(source)
            env.update(ROOT=root,OUT=root/'results',CASE=name)
            try:
                p=fn();manifest=json.loads((p.parent/'manifest.json').read_text())
                self.assertEqual(manifest['probes'][0]['x_um'],16)
                self.assertTrue(manifest['physical_metal_union_XOR_empty'])
                self.assertEqual(len(manifest['probes']),4)
                meta=json.loads(meta_path.read_text());meta['second']['x_meander']=0
                meta['series_probe_T2_x']=16;meta_path.write_text(json.dumps(meta))
                env['OUT']=root/'explicit_bottom_terminal'
                explicit=fn();other=json.loads((explicit.parent/'manifest.json').read_text())
                self.assertEqual(other['probes'][0]['x_um'],16)
            finally:env.update(saved)

    def test_static_cap_extraction_recovers_independent_known_network(self):
        f=np.array([.5,1.,2.])*1e9;c=known_cap()
        extracted,_=A['extract_capacitance'](f,known_z(f,c))
        np.testing.assert_allclose(extracted,c,rtol=1e-10,atol=1e-9)

    def test_shorting_gap_ports_matches_independent_nodal_LC_analysis(self):
        f=np.linspace(4.5e9,5.5e9,151);w=2*np.pi*f;c=known_cap()
        ell=np.diag([.74,6.014,4.923,4.923])
        ell[0,1]=ell[1,0]=.24;ell[0,2]=ell[2,0]=.06;ell[1,2]=ell[2,1]=.10
        y=1j*w[:,None,None]*c*1e-15
        y[:,1:,1:]+=np.linalg.inv(ell*1e-9)[None,:,:]/(1j*w[:,None,None])
        pp=[0,4];ii=[1,2,3]
        effective=y[:,pp,:][:,:,pp]-y[:,pp,:][:,:,ii]@np.linalg.solve(y[:,ii,:][:,:,ii],y[:,ii,:][:,:,pp])
        eye=np.tile(np.eye(2),(len(f),1,1));expected=(eye-50*effective)@np.linalg.inv(eye+50*effective)
        np.testing.assert_allclose(A['short_internal_z'](known_z(f,c)),expected,rtol=1e-8,atol=1e-9)

    def test_two_capacitances_recovered_with_mutual_coupling_fixed(self):
        f=np.linspace(4.5e9,5.5e9,301);reference=known_cap()
        truth=[101.5,454.41]
        observed=A['short_internal_z'](known_z(f,known_cap(*truth)))[:,1,0]
        fit,_=A['fit_two_capacitances'](f,observed,known_z(f,reference),reference)
        self.assertEqual(fit['free_parameters'],['C0_T2_fF','Cg_input_T2_fF'])
        np.testing.assert_allclose(fit['x'],truth,atol=1e-5)
        self.assertLess(fit['S21_complex_RMSE'],1e-8)


if __name__=='__main__':unittest.main()
