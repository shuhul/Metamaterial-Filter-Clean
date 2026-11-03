"""Offline regression tests for geometry, saved imports, and extraction guards."""

from dataclasses import replace
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import gdstk
import numpy as np

from geometry import CapParams, MeanderParams, NCellParams, make_cap, make_meander, make_ncell, place
from geometry.validation import validate_gds_metadata, validate_grid_params, require
from circuit_common import single_sweep_resonance, run_simulation
from circuit_models import s21_single, s21_n_cells

BASE = Path(__file__).resolve().parents[1]


class WorkflowTests(unittest.TestCase):
    def test_cap_opening_and_thickness(self):
        cap = make_cap(CapParams())
        self.assertEqual(cap.bounding_box(), ((0,0),(302,600)))
        self.assertEqual(gdstk.inside([(3,300),(299,300),(151,593),(151,300),(151,1)],cap.polygons),
                         (True,True,True,False,False))
        self.assertAlmostEqual(cap.polygons[0].area(), 11260)

    def test_meander_turns_and_gaps_both_parities(self):
        for n in (3,71,72,73):
            m=MeanderParams(n=n)
            c=make_meander(m)
            self.assertEqual(len(c.polygons),1)
            self.assertEqual(c.bounding_box(),((0,0),(190,596)))
            for i in range(n-1):
                y=m.l-m.l_top-m.t-i*(m.t+m.d)-m.d/2
                x=m.w-m.t/2 if i%2==0 else m.t/2
                self.assertEqual(gdstk.inside([(m.w/2,y),(x,y)],c.polygons),(False,True))

    def test_array_ports_clearance_ground(self):
        for n in (1,4,10):
            p=NCellParams(n_cells=n); c=make_ncell(p)
            self.assertEqual(c.bounding_box(),p.box)
            self.assertEqual(len(c.references),n+4)
            ground=c.references[-2].get_polygons(); feeds=c.references[-1].get_polygons()
            self.assertFalse(gdstk.boolean(ground,feeds,'and'))
            for i,ref in enumerate(c.references[:p.total_cells]):
                cap=ref.cell.references[0].get_polygons()
                for q in cap: q.translate(*ref.origin)
                self.assertFalse(gdstk.boolean(cap,ground,'and'))
                self.assertEqual(bool(gdstk.boolean(cap,feeds,'and')), i in (0,p.total_cells-1))
            self.assertEqual(len(gdstk.boolean(c.get_polygons(),[],'or')),1)

    def test_placement(self):
        c=make_cap(CapParams())
        self.assertEqual(place(c,origin=(0,-8),anchor='top').bounding_box(),((-151,-608),(151,-8)))
        self.assertEqual(tuple(place(c,anchor=(151,608)).origin),(-151,-608))

    def test_grid_rejects_collapsed_clearance(self):
        validate_grid_params(NCellParams())
        for p in (NCellParams(t_ground=1),NCellParams(d_cell=4,t_ground=2),NCellParams(w_side=400)):
            with self.assertRaises(ValueError): validate_grid_params(p)
        with self.assertRaises(ValueError): MeanderParams(n=74)
        with self.assertRaises(ValueError): require(False,'always active')

    def test_saved_gds_matches_metadata(self):
        for n in (1,4):
            p=BASE/'geometry'/'output'/f'ncell_{n}.gds'
            metadata=json.loads(p.with_suffix('.json').read_text())
            validate_gds_metadata(p,metadata)
            metadata['parameters']['unit']['cap']['t_top']=16
            with self.assertRaises(ValueError): validate_gds_metadata(p,metadata)

    def test_n1_cascade_equals_closed_form(self):
        w=2*np.pi*np.linspace(3e9,7e9,2001)
        np.testing.assert_allclose(s21_n_cells(w,190.4,8.2865,4.9194e-9,1),
                                  s21_single(w,190.4e-15,8.2865e-15,4.9194e-9),atol=1e-11)

    def test_peak_tracking_rejects_boundary_and_ambiguity(self):
        for values in ([1,.2,.1,0,0],[0,1,0,1,0],[0,0,0,0,0]):
            s=np.zeros((5,2,2),complex); s[:,1,0]=values
            with self.assertRaises(ValueError): single_sweep_resonance(SimpleNamespace(f=np.arange(5),s=s))
        s=np.zeros((5,2,2),complex); s[:,1,0]=[0,.1,1,.1,0]
        self.assertEqual(single_sweep_resonance(SimpleNamespace(f=np.arange(5),s=s)),2)

    def test_solver_runs_in_project_directory(self):
        with tempfile.TemporaryDirectory() as directory:
            project=Path(directory)/'example.sonx'
            with patch('geometry.taper.solver_guard.run',return_value=SimpleNamespace(returncode=0)) as run:
                self.assertFalse(run_simulation(project))
                self.assertEqual(run.call_args.kwargs['cwd'],project.parent)

    def test_inductance_preparation_without_writing_or_running(self):
        from geometry import prepare_inductance as prep
        from geometry.prepare_sonnet import tag
        with patch.object(prep,'write_xml') as write, patch.object(prep,'run'):
            prep.main()
            tree=write.call_args.args[0]
            component,=tree.findall(f'.//{tag("Component")}')
            self.assertEqual([p.findtext(tag('Center')) for p in component.findall(tag('Pin'))],
                             ['(940,448)','(940,456)'])
            self.assertEqual(component.findtext(tag('Value')),'C_add')


if __name__ == '__main__':
    unittest.main()
