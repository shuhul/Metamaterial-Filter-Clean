"""Verify that an isolated L screen cannot launch a full response."""
import json,tempfile,unittest
from pathlib import Path
from dataclasses import asdict
from unittest.mock import patch
import numpy as np
from geometry.taper import parallel_trials as p

class StageTests(unittest.TestCase):
    def test_inductance_only_records_out_of_tolerance_without_response(self):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'job.json'
            path.write_text(json.dumps(dict(case='screen',settings=asdict(p.w.SETTINGS),inductance_only=True)))
            with patch.object(p.w,'measure',return_value=dict(L0_nH=4.6,target_error_percent=np.float64(6.2))), \
                 patch.object(p.w,'response') as response:
                p.worker(path)
            response.assert_not_called()
            result=json.loads(path.with_suffix('.result.json').read_text())
            self.assertEqual(result['stage'],'inductance_only')
            self.assertFalse(result['inductance_pass'])

    def test_full_response_stays_gated_by_inductance(self):
        with tempfile.TemporaryDirectory() as td:
            path=Path(td)/'job.json'
            path.write_text(json.dumps(dict(case='response',settings=asdict(p.w.SETTINGS))))
            with patch.object(p.w,'measure',return_value=dict(L0_nH=4.6,target_error_percent=6.2)), \
                 patch.object(p.w,'response') as response:
                with self.assertRaisesRegex(ValueError,'Independent L failed'):
                    p.worker(path)
            response.assert_not_called()

if __name__=='__main__':unittest.main()
