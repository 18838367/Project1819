import numpy as np 
import subprocess

def count_jobs(job_name):
    """Returns how many jobs with self.jobs_name are currently queued or running"""

    try:
        out, err, code = exec_command("squeue")
    except OSError:
        raise RuntimeError("Couldn't run qstat, is it installed?")

    if code:
        raise RuntimeError("qstat failed with code %d: stdout: %s, stderr: %s" % (code, out, err))

    lines_with_jobname = [l for l in out.splitlines() if job_name in l]
    return len(lines_with_jobname)

def exec_command(cmd, shell=False, **kwargs):
    """Executes `cmd` and returns the stdout, stderr and exit code"""
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         shell=shell, **kwargs)
    out, err = p.communicate()
    return out, err, p.poll()
