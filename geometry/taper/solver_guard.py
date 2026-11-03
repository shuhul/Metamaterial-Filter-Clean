"""Run only our own child EM process; enforce the declared four-thread cap."""
import ctypes
from ctypes import wintypes
from pathlib import Path
import os
import re
import subprocess
import time
from geometry.taper.solver_resources import solver_environment, solver_lease, ResourceRecord


def run(command,*,cwd,stdout,env=None):
    command=_cache_command(command)
    child_env=solver_environment(env)
    slots=int(child_env.get('SONNET_MAX_JOBS','1'))
    # A queued parallel trial may wait for a healthy full sweep to finish.
    # Waiting reserves no solver threads and never raises the two-job cap.
    with solver_lease(max_jobs=slots,wait_seconds=1800 if slots==2 else 600):
        return _run(command,cwd=cwd,stdout=stdout,env=child_env)


def _cache_command(command):
    """Retain completed EM field caches, respecting an explicit caller choice.

    The inherited template uses STOP_RESTART, which deletes its ABS cache on
    completion. The CLI override avoids changing archived native projects.
    """
    result=list(command)
    if Path(result[0]).name.lower() in ('em','em.exe') and not any(
            str(q).lower().startswith('-abscache') for q in result[1:]):
        result.insert(1,'-AbsCacheMultiSweep')
    return result


def _cpu_seconds(process):
    if os.name!='nt':return None
    creation,exit_time,kernel,user=[wintypes.FILETIME() for _ in range(4)]
    ok=ctypes.windll.kernel32.GetProcessTimes(wintypes.HANDLE(int(process._handle)),
            ctypes.byref(creation),ctypes.byref(exit_time),ctypes.byref(kernel),ctypes.byref(user))
    if not ok:return None
    return ((kernel.dwHighDateTime<<32)+kernel.dwLowDateTime+(user.dwHighDateTime<<32)+user.dwLowDateTime)/1e7


def _run(command,*,cwd,stdout,env):
    process=subprocess.Popen(command,cwd=cwd,stdout=stdout,stderr=subprocess.STDOUT,env=env)
    previous=None;overload=0;resources=None;error_message=None
    try:
        resources=ResourceRecord(process,stdout,cwd,env)
        while process.poll() is None:
            time.sleep(1)
            try:text=Path(stdout.name).read_text(errors='replace')
            except OSError:text=''
            declared=re.findall(r'using\s+(\d+)\s+thread',text,re.I)
            reason=None
            if declared and max(map(int,declared))>4:reason='Sonnet ignored the four-thread cap'
            now=time.monotonic();cpu=_cpu_seconds(process)
            resources.sample(now,cpu,max(map(int,declared)) if declared else None)
            if cpu is not None and previous:
                cores=(cpu-previous[1])/(now-previous[0])
                overload=overload+1 if cores>6 else 0
                if overload>=3:reason='Solver CPU exceeded six core equivalents for three consecutive checks'
            if cpu is not None:previous=(now,cpu)
            if reason:
                process.kill();process.wait()
                print(f'Stopped our solver PID {process.pid}: {reason}',flush=True)
                raise RuntimeError(reason+'; child solver stopped, no next simulation launched')
        return subprocess.CompletedProcess(command,process.returncode)
    except BaseException as error:
        error_message=str(error)
        if process.poll() is None:process.kill();process.wait()
        raise
    finally:
        if resources is not None:resources.finish(error_message)
