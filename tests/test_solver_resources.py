"""Resource controls tested with idle Python children; never launch Sonnet."""
from pathlib import Path
import json
import os
import sys
import tempfile
import unittest
from unittest.mock import patch
from geometry.taper import solver_guard as guard
from geometry.taper import solver_resources as resources


class SolverResourceTests(unittest.TestCase):
    def test_two_job_budget_excludes_third_and_serial_callers(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'slots.lock'
            with resources.solver_lease(path,wait_seconds=0,max_jobs=2):
                with self.assertRaises(RuntimeError):
                    with resources.solver_lease(path,wait_seconds=0,max_jobs=1):pass
                with resources.solver_lease(path,wait_seconds=0,max_jobs=2):
                    with self.assertRaises(RuntimeError):
                        with resources.solver_lease(path,wait_seconds=0,max_jobs=2):pass
            with resources.solver_lease(path,wait_seconds=0,max_jobs=1):
                with self.assertRaises(RuntimeError):
                    with resources.solver_lease(path,wait_seconds=0,max_jobs=2):pass

    def test_completed_EM_caches_are_retained_without_overriding_explicit_choices(self):
        original=['em.exe','-v','test.sonx']
        self.assertEqual(guard._cache_command(original),['em.exe','-AbsCacheMultiSweep','-v','test.sonx'])
        self.assertEqual(original,['em.exe','-v','test.sonx'])
        explicit=['em.exe','-AbsCacheStopRestart','test.sonx']
        self.assertEqual(guard._cache_command(explicit),explicit)
        self.assertEqual(guard._cache_command(['em.exe','-AbsNoDiscrete','test.sonx']),
                         ['em.exe','-AbsCacheMultiSweep','-AbsNoDiscrete','test.sonx'])
        self.assertEqual(guard._cache_command([sys.executable,'-c','pass']),[sys.executable,'-c','pass'])

    def test_environment_enforces_child_limits_without_global_changes(self):
        before=dict(os.environ)
        env=resources.solver_environment({'MKL_NUM_THREADS':'36','OMP_NUM_THREADS':'36','EXAMPLE':'kept'})
        self.assertEqual(env['MKL_NUM_THREADS'],'4')
        self.assertEqual(env['OMP_NUM_THREADS'],'4')
        self.assertEqual(env['OPENBLAS_NUM_THREADS'],'1')
        self.assertEqual(env['EXAMPLE'],'kept')
        self.assertTrue(Path(env['SON_THREAD_CONTROL_INI']).is_file())
        self.assertEqual(dict(os.environ),before)

    def test_exclusive_lock_and_exception_cleanup(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'slot.lock'
            with resources.solver_lease(path,wait_seconds=0):
                with self.assertRaisesRegex(RuntimeError,'resource lock'):
                    with resources.solver_lease(path,wait_seconds=0):pass
            with self.assertRaisesRegex(ValueError,'test interruption'):
                with resources.solver_lease(path,wait_seconds=0):raise ValueError('test interruption')
            with resources.solver_lease(path,wait_seconds=0):pass

    def test_host_saturation_does_not_stop_healthy_child(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory);calls=iter((0,100,200,300,400,500,600))
            def busy_system():
                n=next(calls,700)
                return (0,n,n)  # Synthetic 100% host use, no real CPU stress.
            code="import os,time; print(os.environ['OMP_NUM_THREADS'],flush=True); time.sleep(1.2)"
            with patch.object(resources,'LEASE_PATH',root/'slot.lock'), \
                 patch.object(resources,'_system_times',side_effect=busy_system), \
                 (root/'idle.log').open('w') as log:
                result=guard.run([sys.executable,'-c',code],cwd=root,stdout=log,env={})
            self.assertEqual(result.returncode,0)
            self.assertEqual((root/'idle.log').read_text().strip(),'4')
            record=json.loads((root/'idle.log.resources.json').read_text())
            self.assertEqual(record['state'],'completed')
            self.assertEqual(record['samples'][0]['host_busy_percent'],100)
            self.assertEqual(record['thread_environment']['MKL_NUM_THREADS'],'4')

    def test_only_our_child_is_stopped_on_measured_cpu_fault(self):
        with tempfile.TemporaryDirectory() as directory:
            root=Path(directory)
            with patch.object(resources,'LEASE_PATH',root/'slot.lock'), \
                 patch.object(guard,'_cpu_seconds',side_effect=[0,100,200,300,400,500]), \
                 (root/'fault.log').open('w') as log:
                with self.assertRaisesRegex(RuntimeError,'six core equivalents'):
                    guard.run([sys.executable,'-c','import time; time.sleep(30)'],cwd=root,stdout=log)
            record=json.loads((root/'fault.log.resources.json').read_text())
            self.assertEqual(record['state'],'stopped')
            self.assertIsNotNone(record['returncode'])
            # OS lock released after the fault; next job can obtain the slot.
            with resources.solver_lease(root/'slot.lock',wait_seconds=0):pass


if __name__=='__main__':unittest.main()
