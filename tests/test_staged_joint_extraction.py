"""Recover deliberately off-target capacitors from two independent cascades."""
import json
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch
import unittest
import numpy as np
import skrf
from geometry.taper import staged_joint_extraction as joint
from geometry.taper.staged_fitting import circuit,LT1,LR
from geometry.taper.corrected_iteration import REGULAR_C,REGULAR_G


class JointExtractionTests(unittest.TestCase):
    def test_two_stage_recovery_with_fixed_inductances_and_wrong_seeds(self):
        truth=np.array([174.,723.,121.,34.,9.4])
        with TemporaryDirectory() as directory:
            root=Path(directory);case='synthetic'
            for j,stage in enumerate(('direct','idc')):
                folder=root/case/stage/'fit';folder.mkdir(parents=True)
                ls=[.748+0.001*j,LT1,LR,LR]
                f=np.linspace(4.65e9,5.4e9,1501)
                s=circuit(f,ls,[truth[0] if j else 0,truth[2],REGULAR_C,REGULAR_C],
                    [truth[1] if j else None,truth[3],truth[4],REGULAR_G],left_ps=4.1+3.2*j)
                np.savez(folder/'response.npz',f=f,actual=s)
                source=folder/'source.s2p'
                skrf.Network(frequency=skrf.Frequency.from_f(f,unit='hz'),s=s,z0=50).write_touchstone(str(source))
                # These start values deliberately differ from the truth.
                r=dict(L_nH=ls,C_fF=[50 if j else 0,130,REGULAR_C,REGULAR_C],
                    G_fF=[500 if j else None,28,8.6,REGULAR_G],delay_ps=5,
                    provenance={'source':str(source),'sha256':joint.digest(source)})
                (folder/'fit.json').write_text(json.dumps(r))
            with patch.object(joint,'OUT',root):result=joint.run(case)
            actual=np.array([result['fitted'][k] for k in
                ('C0_T2_fF','Cg_input_T2_fF','C0_T1_fF','Cg_T2_T1_fF','Cg_T1_regular_fF')])
            np.testing.assert_allclose(actual,truth,rtol=.01,atol=.02)
            self.assertLess(max(result['complex_RMS']),1e-5)
            self.assertEqual(result['fixed_inductances_nH'][0],[.748,LT1,LR,LR])
            self.assertEqual(result['fixed_inductances_nH'][1],[.749,LT1,LR,LR])
            with patch.object(joint,'OUT',root):full=joint.run(case,reflections=True,window_GHz=(4.8,5.25))
            actual=np.array([full['fitted'][k] for k in
                ('C0_T2_fF','Cg_input_T2_fF','C0_T1_fF','Cg_T2_T1_fF','Cg_T1_regular_fF')])
            np.testing.assert_allclose(actual,truth,rtol=.01,atol=.02)
            self.assertLess(abs(full['common_right_delay_ps']),.001)
            self.assertLess(max(q for stage in full['frequency_holdout'] for q in stage.values()),1e-5)
            control=root/case/'known_capacitor_control';control.mkdir()
            (control/'open_validation.json').write_text(json.dumps(dict(passed=True,scope='Synthetic test')))
            loads=[]
            for cap in (0,400,800,1200):
                s=circuit(f,[.749,LT1,LR,LR],[truth[0]+cap,truth[2],REGULAR_C,REGULAR_C],
                    [truth[1],truth[3],truth[4],REGULAR_G],left_ps=7.3)
                path=control/f'load_{cap}.s2p'
                skrf.Network(frequency=skrf.Frequency.from_f(f,unit='hz'),s=s,z0=50).write_touchstone(str(path))
                loads.append(dict(known_C_fF=cap,path=str(path),sha256=joint.digest(path)))
            (control/'loads.json').write_text(json.dumps(loads))
            with patch.object(joint,'OUT',root):loaded=joint.run(case,known_loads=True)
            actual=np.array([loaded['fitted'][k] for k in
                ('C0_T2_fF','Cg_input_T2_fF','C0_T1_fF','Cg_T2_T1_fF','Cg_T1_regular_fF')])
            np.testing.assert_allclose(actual,truth,rtol=.01,atol=.02)
            self.assertLess(loaded['held_out_800fF']['complex_RMS'],1e-5)


if __name__=='__main__':unittest.main()
