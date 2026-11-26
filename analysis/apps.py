"""App registry: re-uses the run arguments of scripts/perf_final.sh (single source of truth for inputs and sizes)."""
import os, re, shlex
HERE = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CPUSET = "0-63,96-159"; NTHR = 128
def load():
    sub = {"$HERE": HERE, "$TO": f"{HERE}/omp/test", "$TC": f"{HERE}/opencl/test", "$D": f"{HERE}/data/inputs"}
    apps = {}
    for line in open(os.path.join(HERE, "scripts", "perf_final.sh")):
        if not line.startswith("reg "): continue
        for k, v in sub.items(): line = line.replace(k, v)
        t = shlex.split(line)
        apps[t[1]] = dict(name=t[1], omp_dir=t[2], ocl_dir=t[3], omp_args=shlex.split(t[4]), ocl_args=t[5])
    return apps
APPS = load()

ACC_APPS = {"swat"}  # apps that use dwarf_acc_* spans instead of region_begin/end (see common/dwarf_timer.h)
