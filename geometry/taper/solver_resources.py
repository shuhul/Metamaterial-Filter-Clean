"""Child-only thread settings, OS-owned solver lock and CPU evidence."""
from contextlib import contextmanager
import ctypes
from ctypes import wintypes
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import time

LEASE_PATH = Path(__file__).resolve().parents[2]/'output_data/solver_runtime/active.lock'


def solver_environment(base=None):
    env = os.environ.copy()
    if base is not None:
        env.update(base)
    env.update(SON_THREAD_CONTROL_INI=str(Path(__file__).with_name('sonnet_threads.ini')),
               MKL_NUM_THREADS='4', OMP_NUM_THREADS='4', OPENBLAS_NUM_THREADS='1',
               NUMEXPR_NUM_THREADS='1', KMP_BLOCKTIME='0')
    return env


@contextmanager
def solver_lease(path=None, wait_seconds=600, max_jobs=None):
    """Queue concurrent workflow launches. OS releases the lock after a crash.

    An old lock file is harmless: file existence is never used as liveness.
    The lock does not cover independent GUI jobs or other project folders.
    """
    # Serial callers reserve both slots. Explicit two-job callers reserve one.
    # This keeps serial and parallel workflows under the same global budget.
    count=int(os.environ.get('SONNET_MAX_JOBS','1')) if max_jobs is None else max_jobs
    if count not in (1,2):raise ValueError('Only one or two four-thread Sonnet jobs are allowed')
    path=Path(path or LEASE_PATH);path.parent.mkdir(parents=True,exist_ok=True)
    paths=[path,path.with_name(path.stem+'.slot1'+path.suffix)]
    def unlock(h):
        h.seek(0)
        if os.name=='nt':
            import msvcrt
            msvcrt.locking(h.fileno(),msvcrt.LK_UNLCK,1)
        else:
            import fcntl
            fcntl.flock(h,fcntl.LOCK_UN)
        h.close()
    start=time.monotonic();announced=False;held=[]
    while not held:
        for candidate in ([paths] if count==1 else [[p] for p in paths]):
            acquired=[]
            try:
                for p in candidate:
                    h=p.open('a+b');h.seek(0,2)
                    if not h.tell():h.write(b'\0');h.flush()
                    h.seek(0)
                    try:
                        if os.name=='nt':
                            import msvcrt
                            msvcrt.locking(h.fileno(),msvcrt.LK_NBLCK,1)
                        else:
                            import fcntl
                            fcntl.flock(h,fcntl.LOCK_EX|fcntl.LOCK_NB)
                    except BaseException:h.close();raise
                    acquired.append(h)
                held=acquired;break
            except OSError:
                for h in reversed(acquired):unlock(h)
            except BaseException:
                for h in reversed(acquired):unlock(h)
                raise
        if not held:
            if time.monotonic()-start>=wait_seconds:
                raise RuntimeError('Another workflow solver holds the resource lock; no new solver launched')
            if not announced:
                print('Sonnet resource slots occupied; waiting.',flush=True);announced=True
            time.sleep(1)
    try:yield
    finally:
        for h in reversed(held):unlock(h)


def _utc():
    return datetime.now(timezone.utc).isoformat()


def _system_times():
    if os.name != 'nt':
        return None
    values = [wintypes.FILETIME() for _ in range(3)]
    ok = ctypes.windll.kernel32.GetSystemTimes(*(ctypes.byref(v) for v in values))
    return tuple((v.dwHighDateTime << 32) | v.dwLowDateTime for v in values) if ok else None


class ResourceRecord:
    """CPU samples are evidence, never a reason to stop a healthy child."""
    def __init__(self, process, log, cwd, env):
        self.path = Path(str(log.name)+'.resources.json')
        self.process = process
        self.start = time.monotonic()
        self.previous = (self.start, None, _system_times())
        self.last_saved = self.start; self.last_warning = self.start-60
        self.record = {'started_utc': _utc(), 'state': 'running',
                       'launcher_pid': os.getpid(), 'child_pid': process.pid,
                       'working_directory': str(cwd), 'thread_limit': 4,
                       'logical_processors': os.cpu_count(),
                       'max_concurrent_jobs':int(env.get('SONNET_MAX_JOBS','1')),
                       'thread_environment': {key: env[key] for key in (
                           'SON_THREAD_CONTROL_INI','MKL_NUM_THREADS','OMP_NUM_THREADS',
                           'OPENBLAS_NUM_THREADS','NUMEXPR_NUM_THREADS','KMP_BLOCKTIME')},
                       'policy': 'Host saturation warns; only child thread/CPU faults stop the child.',
                       'samples': [], 'peak_child_core_equivalents': 0.0}
        self.save()

    def save(self):
        self.path.write_text(json.dumps(self.record, indent=2)+'\n')

    def sample(self, now, cpu, declared_threads):
        system = _system_times()
        point = {'utc': _utc(), 'elapsed_seconds': now-self.start,
                 'child_cpu_seconds': cpu, 'child_core_equivalents': None,
                 'host_busy_percent': None}
        before, old_cpu, old_system = self.previous
        if cpu is not None and old_cpu is not None and now > before:
            point['child_core_equivalents'] = (cpu-old_cpu)/(now-before)
            self.record['peak_child_core_equivalents'] = max(
                self.record['peak_child_core_equivalents'], point['child_core_equivalents'])
        if system is not None and old_system is not None:
            idle, kernel, user = [b-a for a,b in zip(old_system, system)]
            if kernel+user > 0:
                point['host_busy_percent'] = 100*(1-idle/(kernel+user))
        self.record['declared_threads'] = declared_threads
        busy = point['host_busy_percent']
        if busy is not None and busy >= 90 and now-self.last_warning >= 60:
            print(f'Host CPU {busy:.0f}%: continuing our four-thread job, PID {self.process.pid}; '
                  'resource measurements saved beside its log.', flush=True)
            self.last_warning = now
        if now-self.last_saved >= 10 or not self.record['samples']:
            self.record['samples'].append(point)
            self.save(); self.last_saved = now
        self.previous = (now, cpu, system)

    def finish(self, error=None):
        self.record.update(state='stopped' if error else (
            'completed' if self.process.returncode == 0 else 'failed'),
            ended_utc=_utc(), elapsed_seconds=time.monotonic()-self.start,
            returncode=self.process.returncode, error=error,
            final_child_cpu_seconds=self.previous[1])
        self.save()
