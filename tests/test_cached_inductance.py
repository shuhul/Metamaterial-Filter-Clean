"""Cache-only calibration must explicitly forbid direct EM computations."""
from pathlib import Path
from tempfile import TemporaryDirectory
from types import SimpleNamespace
import subprocess
import unittest
from unittest.mock import patch
import numpy as np
import extract_generated_inductance as extraction


class CachedInductanceTests(unittest.TestCase):
    def exercise(self,frequencies):
        with TemporaryDirectory() as directory:
            root=Path(directory);project=root/'test.sonx';project.write_text('archivable project')
            signal=np.zeros((5,2,2),complex);signal[:,1,0]=[0,.1,1,.1,0]
            net=SimpleNamespace(f=np.linspace(4e9,6e9,5),s=signal)
            def cached(command,**kwargs):
                self.assertIn('-AbsNoDiscrete',command)
                kwargs['stdout'].write(f'Em errors detected: 0\nTotal time for {frequencies} frequencies: 1 second.\nEm simulation completed\n')
                project.with_suffix('.snp').write_text('!< C_add = 0.8\n')
                return subprocess.CompletedProcess(command,0)
            with patch.object(extraction,'write_sonx'),patch.object(extraction,'read_snp',return_value=net),\
                 patch('geometry.taper.solver_guard.run',side_effect=cached):
                return extraction.simulate(.8,'cached',project_file=project,output_dir=root,cache_only=True)

    def test_cache_flag_and_zero_frequency_confirmation(self):
        _,frequency=self.exercise(0)
        self.assertEqual(frequency,5e9)

    def test_missing_cache_confirmation_rejected(self):
        with self.assertRaisesRegex(RuntimeError,'zero new EM frequencies'):self.exercise(1)


if __name__=='__main__':unittest.main()
