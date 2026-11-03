"""Offline checks for new taper geometry and nonuniform circuit physics."""
import json
from pathlib import Path
import tempfile
import unittest
import gdstk
import numpy as np
from geometry.taper.layout import TaperParams,build,write,validate
from geometry.taper.fitting import chain
from circuit_models import s21_n_cells

class TaperTests(unittest.TestCase):
    def test_uniform_limit_and_energy_conservation(self):
        omega=2*np.pi*np.linspace(4e9,6e9,151)
        s=chain(omega,[4.9194]*2,[190.4]*2,[8.2865]*3)
        np.testing.assert_allclose(s[:,1,0],s21_n_cells(omega,190.4,8.2865,4.9194e-9,2),rtol=1e-10,atol=1e-12)
        s=chain(omega,[6.06,4.9194],[100,190.4],[56.7,8.2865,8.2865])
        np.testing.assert_allclose(np.swapaxes(s.conj(),1,2)@s,np.broadcast_to(np.eye(2),s.shape),atol=1e-10)
    def test_continuous_top_boundary_and_grounded_meanders(self):
        top,m=build(TaperParams());polys=top.get_polygons();centers=m['cell_centers']
        probes=[(x,1) for x in np.linspace(2,m['w_sim']-2,60)]
        self.assertTrue(all(gdstk.inside(probes,polys)))
        self.assertFalse(any(gdstk.inside([(x,-2) for x in centers],polys)))
        self.assertTrue(all(gdstk.inside([(x,-618) for x in centers],polys)))
    def test_rebuilt_geometry_rejects_stale_metadata(self):
        with tempfile.TemporaryDirectory() as d:
            path=write(TaperParams(),Path(d)/'taper');m=json.loads(path.with_suffix('.json').read_text())
            validate(path,m)
            m['taper']['w']=332
            with self.assertRaises(ValueError):validate(path,m)
    def test_extracts_capacitances_without_absorbing_fixture_phase(self):
        import contextlib,io
        import skrf as rf
        from geometry.taper.fitting import fit
        freq=np.linspace(4e9,5.8e9,901);omega=2*np.pi*freq
        s=chain(omega,[6.06,4.922723945068658],[100,188.8781163861302],[56.7,8.2865,8.2865])
        s=s*np.exp(-1j*omega*4e-12)[:,None,None]
        with tempfile.TemporaryDirectory() as d,contextlib.redirect_stdout(io.StringIO()):
            path=Path(d)/'synthetic'
            net=rf.Network(frequency=rf.Frequency.from_f(freq,unit='Hz'),s=s,z0=50)
            net.write_touchstone(str(path))
            r=fit(path.with_suffix('.s2p'),Path(d)/'fit',l_taper=6.06,l_regular=4.922723945068658,c_regular=188.8781163861302,free_regular=True)['complex_with_delay']
            for key,target in [('C_taper_fF',100),('Cg_input_fF',56.7),('Cg_junction_fF',8.2865),('effective_fixture_delay_ps',4)]:
                self.assertAlmostEqual(r[key],target,places=4)
    def test_cache_only_rejects_an_older_export(self):
        from unittest.mock import patch
        import shutil,subprocess
        from geometry.taper import sonnet
        root=Path(__file__).resolve().parents[1]
        source=root/'sonnet_files/taper_exploration/final_designs/finger536_t20_right310_input6_j56_final/generated_ncell_2.sonx'
        with tempfile.TemporaryDirectory() as d:
            project=Path(d)/'case.sonx';shutil.copyfile(source,project)
            shutil.copyfile(source.with_suffix('.snp'),project.with_suffix('.snp'))
            sonnet.configure(project,4,5.6,3200)
            def cached(*args,**kwargs):
                kwargs['stdout'].write('Total time for 0 frequencies: 0 seconds.\nEm errors detected: 0\nEm simulation completed\n')
                return subprocess.CompletedProcess([],0)
            with patch.object(sonnet,'guarded_solver',side_effect=cached):
                with self.assertRaisesRegex(RuntimeError,'export does not match'):
                    sonnet.simulate(project,Path(d)/'result',cache_only=True)
            self.assertEqual(json.loads((Path(d)/'result/run_status.json').read_text())['status'],'failed')
    def test_input_finger_matches_native_benchmark_topology(self):
        top,m=build(TaperParams(),input_style='finger')
        self.assertEqual(m['expected_metal_components'],2)
        x=m['cell_centers'][0]
        self.assertFalse(gdstk.inside([(x,-100)],top.get_polygons())[0])
        self.assertTrue(gdstk.inside([(x+152,-100)],top.get_polygons())[0])
        self.assertEqual(len(gdstk.boolean(top.get_polygons(),[],'or',precision=.001)),2)
    def test_asymmetric_cap_does_not_touch_bottom(self):
        top,m=build(TaperParams(),isolation=True)
        polys=top.get_polygons();center=m['cell_centers'][1]
        self.assertFalse(gdstk.inside([(center+162,-615)],polys)[0])
        self.assertEqual(len(gdstk.boolean(polys,[],'or',precision=.001)),2)

if __name__=='__main__':unittest.main()
