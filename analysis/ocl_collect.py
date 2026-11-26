#!/usr/bin/env python3
"""Collect an evidence bundle for one OpenCL (app, runtime) using the LD_PRELOAD interposer analysis/clshim.c. Works for PoCL, Intel's CPU runtime and NVIDIA.
   python3 -I analysis/ocl_collect.py --app crc --config ocl-intel --out DIR [--reps 5 --warmup 1]
   configs: ocl-nvidia (GPU 0), ocl-pocl (CPU, 128 threads), ocl-intel (CPU, 128 threads). Bundle schema = gpu_collect.py (so compare.py works on both);
   extra keys: idle_ms_per_rep (device idle time between consecutive launches on the same queue) and runtime (the interposer's raw class totals)."""
import argparse, glob, json, os, re, subprocess, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import apps as A, budget as BG
HERE = A.HERE; SHIM = os.path.join(HERE, "analysis", "libclshim.so"); MARK = os.path.join(HERE, "analysis", "libmark.so")
VENDOR = {"ocl-nvidia": ("nvidia", "gpu"), "ocl-pocl": ("pocl", "cpu"), "ocl-intel": ("intel", "cpu")}
MAP = {"launch": ["enqueue_ndrange", "set_kernel_arg"], "copy_h2d": ["write_buffer"], "copy_d2h": ["read_buffer"], "copy_other": ["copy_buffer", "fill_buffer", "map_buffer", "unmap"],
       "sync": ["finish", "wait_events", "flush"], "alloc": ["create_buffer"], "free": ["release_mem"], "init": ["init"]}
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--app", required=True); ap.add_argument("--config", required=True, choices=VENDOR); ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=5); ap.add_argument("--warmup", type=int, default=1); ap.add_argument("--threads", type=int, default=A.NTHR, help="CPU worker threads (affinity set = first N cpus of the pinned set)"); ap.add_argument("--tag", default="", help="suffix for output file names"); ap.add_argument("--build", action="store_true", help="(re)build libclshim.so first")
    o = ap.parse_args(); o.out = os.path.abspath(o.out); os.makedirs(o.out, exist_ok=True)
    if o.build or not os.path.exists(SHIM):
        inc = [f"-I{d}" for d in (os.path.expanduser("~/.local/apps/opencl-hudson/include"), "/usr/include") if os.path.isdir(os.path.join(d, "CL"))][:1]
        subprocess.run(["gcc", "-shared", "-fPIC", "-O2", "-Wall", "-o", SHIM, os.path.join(HERE, "analysis", "clshim.c")] + inc + ["-ldl", "-lpthread"], check=True)
    vend, dev = VENDOR[o.config]; a = A.APPS[o.app]; tag = f"{o.app}.{o.config}{o.tag}"; base = os.path.join(o.out, tag + ".shim")
    for f in glob.glob(base + ".*"): os.remove(f)
    env = dict(os.environ); env.update(LD_PRELOAD=SHIM + ":" + MARK, CLSHIM_OUT=base, MARK_OUT=base + ".marks", DWARF_WARMUP=str(o.warmup), DWARF_REPEATS=str(o.reps))
    for f in glob.glob(base + ".marks.*"): os.remove(f)
    if o.app in A.ACC_APPS: env["MARK_ACC"] = "1"
    if dev == "gpu": env["CUDA_VISIBLE_DEVICES"] = "0"
    else: env.update(POCL_MAX_PTHREAD_COUNT=str(o.threads), OMP_NUM_THREADS=str(o.threads))
    cpus = [c for part in A.CPUSET.split(",") for c in (range(int(part.split("-")[0]), int(part.split("-")[1]) + 1) if "-" in part else [int(part)])]
    cmd = ["taskset", "-c", ",".join(map(str, cpus[:o.threads])), "make", "-s", "-C", f"{HERE}/opencl/{a['ocl_dir']}", "run", f"OCL_VENDOR={vend}", "GPU_INDEX=0", f"DWARF_REPEATS={o.reps}", f"DWARF_WARMUP={o.warmup}", f"RUN_ARGS={a['ocl_args']}"]
    t0 = time.time(); r = subprocess.run(cmd, env=env, capture_output=True, text=True); log = r.stdout + r.stderr; open(os.path.join(o.out, tag + ".log"), "w").write(log)
    m = re.search(r"DWARF_RESULT,.*", log)
    if not m: sys.exit(f"no DWARF_RESULT for {tag}; see {o.out}/{tag}.log")
    res = dict(kv.split("=", 1) for kv in m.group(0).split(",")[1:] if "=" in kv)
    files = sorted([f for f in glob.glob(base + ".*") if not f.endswith(".tl")], key=os.path.getsize)
    if not files: sys.exit(f"interposer wrote no output for {tag} (is the app linked dynamically against libOpenCL?)")
    sh = json.load(open(files[-1])); nreps = o.reps + o.warmup
    per = lambda ns: ns / 1e6 / nreps
    buckets = {b: per(sum(sh["api"][c]["total_ns"] for c in cs)) for b, cs in MAP.items()}
    counts = {b: sum(sh["api"][c]["count"] for c in cs) / nreps for b, cs in MAP.items()}
    kernels = {}
    for k in sh["kernels"]:
        d = kernels.setdefault(k["name"], dict(raw_name=k["name"], instances=0, total_ns=0.0, regs=0, grids=[], blocks=[], q2s_ns=0.0, idle_ns=0.0))
        d["instances"] += k["count"]; d["total_ns"] += k["dev_ns"]; d["q2s_ns"] += k["queued_to_start_ns"]; d["idle_ns"] += k["idle_before_ns"]
        g = [x for x in k["global"]]; l = [x for x in k["local"]]
        if g not in d["grids"]: d["grids"].append(g)
        if l not in d["blocks"]: d["blocks"].append(l)
    for d in kernels.values():
        d["grids"] = d["grids"][:4]; d["blocks"] = d["blocks"][:4]; d["ms_per_rep"] = per(d["total_ns"]); d["inst_per_rep"] = d["instances"] / nreps; d["idle_ms_per_rep"] = per(d["idle_ns"])
    tl = files[-1] + ".tl"; ivs = {}; kiv = []; names = {}
    if os.path.exists(tl):
        for line in open(tl):
            t = line.split()
            if line.startswith("# classes:"): names = {int(x.split("=")[0]): x.split("=")[1] for x in t[2:]}
            elif t[0] == "C":
                b = next((bb for bb, cs in MAP.items() if names.get(int(t[1])) in cs), "other"); ivs.setdefault(b, []).append((int(t[2]), int(t[3])))
            elif t[0] == "E": kiv.append((int(t[1]), int(t[2])))
    win = []; mk = sorted(glob.glob(base + ".marks.*"), key=os.path.getsize)
    if mk:
        pend = None; spans = []
        for line in open(mk[-1]):
            t = line.split()
            if t[0] == "B": pend = int(t[1])
            elif t[0] == "E" and pend is not None:
                if int(t[2]) >= 0: win.append((pend, int(t[1])))
                pend = None
            elif t[0] == "S" and pend is not None: spans.append((pend, int(t[1]))); pend = None
            elif t[0] == "K":
                if int(t[2]) >= 0: win += spans
                spans = []
    part = BG.partition(ivs, kiv, [], window_iv=win or None); nreps_p = o.reps if win else nreps
    bundle = dict(app=o.app, config=o.config, date=time.strftime("%F %T"), reps=o.reps, warmup=o.warmup, wall_ms=float(res["median_ms"]), buckets_ms_per_rep=buckets, bucket_counts_per_rep=counts,
                  kernel_ms_per_rep=per(sum(d["total_ns"] for d in kernels.values())), launches_per_rep=sum(d["instances"] for d in kernels.values()) / nreps,
                  idle_ms_per_rep=per(sum(d["idle_ns"] for d in kernels.values())), api={}, kernels=kernels, memcpy={}, lowered={}, runtime=sh["api"],
                  toolchain=res.get("toolchain", ""), valid=res.get("valid", ""), events_without_profile=sh.get("events_without_profile", 0),
                  partition_ms_per_rep={k: v / 1e6 / nreps_p for k, v in part["buckets_ns"].items()}, window_ms_per_rep=part["window_ns"] / 1e6 / nreps_p, window_is_timed_region=bool(win), init_in_window_ms=part["init_in_window_ns"] / 1e6)
    json.dump(bundle, open(os.path.join(o.out, tag + ".bundle.json"), "w"), indent=1)
    print(f"{tag}: wall={bundle['wall_ms']:.3f} ms  kernel(dev)={bundle['kernel_ms_per_rep']:.3f}  device_idle={bundle['idle_ms_per_rep']:.3f} ms/rep  launches/rep={bundle['launches_per_rep']:.0f}  valid={bundle['valid']}  ({time.time()-t0:.0f}s)")
main()
