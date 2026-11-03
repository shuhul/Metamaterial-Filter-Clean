"""A crashed native solver must not leave a misleading 'running' marker."""
import json
from pathlib import Path
import subprocess
import tempfile
import unittest
from unittest.mock import patch
from geometry.taper.sonnet import simulate


class FailedRunStatusTests(unittest.TestCase):
    def test_native_crash_is_recorded_and_propagated(self):
        with tempfile.TemporaryDirectory() as folder:
            root=Path(folder);project=root/'diagnostic.sonx';project.write_text('<not-read/>')
            output=root/'results'
            result=subprocess.CompletedProcess(['em.exe'],3221225477)
            with patch('geometry.taper.sonnet.guarded_solver',return_value=result):
                with self.assertRaises(subprocess.CalledProcessError):simulate(project,output)
            status=json.loads((output/'run_status.json').read_text())
            self.assertEqual(status['status'],'failed')
            self.assertEqual(status['returncode'],3221225477)
            self.assertIn('ended_utc',status)
            self.assertFalse((output/'response.snp').exists())
