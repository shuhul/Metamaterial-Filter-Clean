"""Protect the verified downstream geometry and the two-variable fit contract."""
from dataclasses import replace
from pathlib import Path
from tempfile import TemporaryDirectory
import json
import os
import sys
import unittest
from unittest.mock import patch
import numpy as np
import gdstk
from geometry.taper.run_locked_taper2 import GEOMETRY
from geometry.taper.locked_taper2 import T1,FIXED,require_locked,locked_validate,response,fit
from geometry.taper.second_layout import write,require_same_isolation_geometry
from geometry.taper.solver_guard import run as guarded_solver
from geometry.taper import locked_taper2
from geometry.taper import solver_resources


class LockedTaperTests(unittest.TestCase):
    def test_compact_idc_preserves_reference_box_and_T1_position(self):
        from geometry.taper.locked_taper2 import VERIFIED
        from geometry.taper.second_layout import build
        original=json.loads(Path(VERIFIED['artifacts']['fixture_gds']).with_suffix('.json').read_text())
        top,meta=build(GEOMETRY)
        self.assertEqual(meta['w_sim'],original['w_sim'])
        self.assertEqual(meta['cell_centers'][2:],original['cell_centers'][1:])
        # The extension starts 4 um beyond the input trunk: it never shorts
        # the two IDC terminals, despite their aligned top edges.
        self.assertEqual(GEOMETRY.idc_span-GEOMETRY.w_top_left_extension-20,4)
        self.assertEqual(len(gdstk.boolean(top.get_polygons(),[],'or',precision=.001)),2)

    def test_straight_mode_has_no_horizontal_kinks(self):
        from geometry.taper.layout import make_taper
        p=GEOMETRY.taper;x=GEOMETRY.x_meander
        metal=make_taper(p,meander_offset_x=x,straight_inductor=True).get_polygons()
        self.assertTrue(all(gdstk.inside([(x,y) for y in (-34,-44,-54,-200,-600)],metal)))
        self.assertFalse(any(gdstk.inside([(x+offset,y) for offset in (-4,4) for y in (-34,-44,-54,-200)],metal)))

    def test_admittance_diagnostics_recover_synthetic_model(self):
        from geometry.taper.diagnose_locked_response import required_shunt,required_junction
        f=np.linspace(4.6e9,5.5e9,451);w=2*np.pi*f;l=.74
        signal=response(f,[100,459],l)[:,1,0]
        np.testing.assert_allclose(required_shunt(f,signal,459),1j*w*100e-15+1/(1j*w*l*1e-9),atol=1e-10)
        np.testing.assert_allclose(required_junction(f,signal,100,459,l),1/(1j*w*FIXED['Cg_input_fF']*1e-15),atol=1e-9)

    def test_mutual_model_reduces_to_chain_and_remains_reciprocal(self):
        from geometry.taper.diagnose_mutual_inductance import coupled_response
        f=np.linspace(4.5e9,5.5e9,251)
        np.testing.assert_allclose(coupled_response(f,[100,459,0],.74),response(f,[100,459],.74),atol=1e-9)
        s=coupled_response(f,[100,459,.02],.74)
        np.testing.assert_allclose(s[:,0,1],s[:,1,0],atol=1e-9)
        # The added coupling is lossless; it must not manufacture gain/loss.
        identity=s@np.swapaxes(s.conj(),1,2)
        np.testing.assert_allclose(identity,np.tile(np.eye(2),(len(f),1,1)),atol=1e-9)

    def test_existing_case_geometry_is_not_rewritten(self):
        with TemporaryDirectory() as directory:
            root=Path(directory);name='immutable_case'
            gds=write(GEOMETRY,root/'geometry/taper/output'/name)
            original=gds.read_bytes()
            with patch.object(locked_taper2,'ROOT',root),patch.object(locked_taper2,'preview'),\
                 patch.object(locked_taper2,'write',side_effect=AssertionError('Must not rewrite existing geometry')):
                locked_taper2.prepare_case(name,GEOMETRY)
                self.assertEqual(gds.read_bytes(),original)
                changed=replace(GEOMETRY,taper=replace(GEOMETRY.taper,l_left=GEOMETRY.taper.l_left-2))
                with self.assertRaisesRegex(ValueError,'new case name'):
                    locked_taper2.prepare_case(name,changed)

    def test_guard_stops_our_child_if_thread_cap_is_ignored(self):
        with TemporaryDirectory() as directory:
            with patch.object(solver_resources,'LEASE_PATH',Path(directory)/'test.lock'), \
                 (Path(directory)/'child.log').open('w') as log:
                with self.assertRaisesRegex(RuntimeError,'four-thread cap'):
                    guarded_solver([sys.executable,'-c',
                        "import time; print('Threads: using 36 thread(s)',flush=True); time.sleep(30)"],
                        cwd=directory,stdout=log,env=os.environ.copy())

    def test_downstream_metal_is_exactly_verified(self):
        with TemporaryDirectory() as directory:
            gds=write(GEOMETRY,Path(directory)/'candidate')
            locked_validate(gds,json.loads(gds.with_suffix('.json').read_text()))

    def test_rejects_T1_and_junction_changes(self):
        for p in (replace(GEOMETRY,taper1=replace(T1,l_right=300)),
                  replace(GEOMETRY,d_t1_regular=48),replace(GEOMETRY,w_upper_plate=100)):
            with self.assertRaises(ValueError):require_locked(p)

    def test_response_phase_can_tune_T2_gap_without_changing_T1_metal(self):
        changed=replace(GEOMETRY,d_t2_t1=2,w_top_left_extension=0)
        with self.assertRaises(ValueError):require_locked(changed)
        require_locked(changed,allow_t2_junction=True)
        with TemporaryDirectory() as directory:
            gds=write(changed,Path(directory)/'response_case')
            locked_validate(gds,json.loads(gds.with_suffix('.json').read_text()),allow_t2_junction=True)

    def test_changed_ground_pedestal_cannot_reuse_L(self):
        with TemporaryDirectory() as directory:
            gds=write(GEOMETRY,Path(directory)/'cal',isolation=True)
            require_same_isolation_geometry(GEOMETRY,gds)
            with self.assertRaisesRegex(ValueError,'different isolation metal'):
                require_same_isolation_geometry(replace(GEOMETRY,l_ground_rise=128),gds)

    def test_synthetic_complex_fit_recovers_only_two_parameters(self):
        frequency=np.linspace(3.8e9,9e9,601);expected=[104.,470.]
        s=response(frequency,expected,.74)
        with self.assertRaises(ValueError):response(frequency,[100,459,5],.74)
        with TemporaryDirectory() as directory:
            path=Path(directory)/'synthetic.snp'
            lines=['# GHz S RI R 50']
            for f,point in zip(frequency,s):
                row=[f/1e9]
                for i,j in ((0,0),(1,0),(0,1),(1,1)):row.extend([point[i,j].real,point[i,j].imag])
                lines.append(' '.join(f'{v:.16g}' for v in row))
            path.write_text('\n'.join(lines)+'\n')
            r=fit(path,Path(directory),.74)
            np.testing.assert_allclose(r['x'],expected,rtol=1e-7)
            self.assertEqual(r['free_parameters'],['C0_t2_fF','Cg_t2_fF'])
            self.assertEqual(r['fixed_model'],FIXED)
            self.assertLess(r['s21_complex_rmse'],1e-10)


if __name__=='__main__':unittest.main()
