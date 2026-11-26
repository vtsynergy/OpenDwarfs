#!/usr/bin/env python3
"""PoCL-side code-generation knobs on kmeans and lud (1 worker thread). python3 -I analysis/pocl_knobs.py OUTDIR [--apps kmeans,lud]
Prediction registered before the run (text below is copied into the result file):
  If PoCL's slowness on kmeans/lud (Intel executes 5.5x / 19.6x fewer instructions with wide vectors) is due to a conservative vector-width choice or cost model, forcing wider vectors
  (POCL_VECTORIZER_PREFER_VECTOR_WIDTH=512, POCL_VECTORIZER_FORCE_VECTOR_WIDTH=8/16) should raise the share of 256/512-bit FP micro-ops, cut instructions executed and time, moving PoCL toward Intel.
  If the loop vectorizer cannot vectorize these work-item loops (legality / dependences), the remarks will say why and the width knobs will not help.
  Controls: FORCE_VECTOR_WIDTH=1 (vectorization off) must not be faster than the default; the default is collected twice per app."""
import json, os, re, shutil, statistics as st, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cpu_lowered as C, apps as A
out = os.path.abspath(sys.argv[1]); os.makedirs(out, exist_ok=True)
apps = (sys.argv[sys.argv.index("--apps") + 1] if "--apps" in sys.argv else "kmeans,lud").split(","); N = 2
V = [("default", {}), ("default (repeat)", {}), ("prefer 512-bit", {"POCL_VECTORIZER_PREFER_VECTOR_WIDTH": "512"}),
     ("force VF 8", {"POCL_VECTORIZER_FORCE_VECTOR_WIDTH": "8"}), ("force VF 16 + prefer 512-bit", {"POCL_VECTORIZER_FORCE_VECTOR_WIDTH": "16", "POCL_VECTORIZER_PREFER_VECTOR_WIDTH": "512"}),
     ("force VF 1 (vectorizer off)", {"POCL_VECTORIZER_FORCE_VECTOR_WIDTH": "1"}), ("method loops (no vectorizer)", {"POCL_WORK_GROUP_METHOD": "loops"}),
     ("method cbs", {"POCL_WORK_GROUP_METHOD": "cbs"}), ("target cpu znver4", {"POCL_LLVM_CPU_NAME": "znver4"}), ("no wg specialization", {"POCL_WORK_GROUP_SPECIALIZATION": "0"})]
def summarize(d):
    w = {k: d.get("fp_ops_retired_by_width." + k + "_uops_retired", 0) for k in ("scalar", "pack_128", "pack_256", "pack_512")}; t = max(sum(w.values()), 1)
    return dict(ms=d.get("median_ms"), instr=d.get("instructions", 0), cycles=d.get("cycles", 0), ipc=d.get("ipc", 0), fpi=d.get("fp_ret_sse_avx_ops.all", 0) / max(d.get("instructions", 1), 1),
                w=[100 * w[k] / t for k in ("scalar", "pack_128", "pack_256", "pack_512")], valid=d.get("valid"))
res = {}
for app in apps:
    for name, env in V:
        cache = os.path.join(out, f"cache-{app}-{re.sub('[^a-z0-9]+', '_', name.lower())}"); shutil.rmtree(cache, ignore_errors=True)
        e = dict(env); e.update(POCL_CACHE_DIR=cache, POCL_LEAVE_KERNEL_COMPILER_TEMP_FILES="1")
        C.run_app(app, "pocl", e, 1, 0, out)                                    # populate the (per-variant) kernel cache
        runs = [summarize(C.run_app(app, "pocl", e, 3, 1, out, perf=True)) for _ in range(N)]
        ks = C.pocl_kernels(cache); feat = {k: dict(instr=v["features"].get("instr", 0), vec=round(v["features"]["vector_pct"]), width=v["features"]["max_width"]) for k, v in ks.items()}
        res[(app, name)] = dict(runs=runs, static=feat)
    # vectorizer remarks (default build) for the explanation of what is / is not vectorized
    cache = os.path.join(out, f"cache-{app}-remarks"); shutil.rmtree(cache, ignore_errors=True)
    a = A.APPS[app]; env = dict(os.environ, POCL_CACHE_DIR=cache, POCL_VECTORIZER_REMARKS="1", POCL_MAX_PTHREAD_COUNT="1", OMP_NUM_THREADS="1", DWARF_REPEATS="1", DWARF_WARMUP="0")
    r = subprocess.run(["taskset", "-c", "0", "make", "-s", "-C", f"{A.HERE}/opencl/{a['ocl_dir']}", "run", "OCL_VENDOR=pocl", f"RUN_ARGS={a['ocl_args']}"], env=env, capture_output=True, text=True)
    open(os.path.join(out, f"remarks-{app}.txt"), "w").write(r.stdout + r.stderr)
json.dump({f"{a}|{n}": v for (a, n), v in res.items()}, open(os.path.join(out, "pocl_knobs.json"), "w"), indent=1)
L = ["# PoCL code-generation knobs on kmeans and lud (1 worker thread, exploratory)", "", *[l.strip() for l in __doc__.split("Prediction registered before the run (text below is copied into the result file):")[1].strip().splitlines()], "",
     f"Each row: median of {N} collections (1 warm-up + 3 timed repetitions, counters over the timed repetitions only); fresh kernel cache per variant. Intel reference (default, same protocol, earlier session): kmeans 1411 ms, lud 637 ms.", ""]
for app in apps:
    L += [f"## {app}", "", "| variant | ms | user instr | IPC | FLOP/instr | FP uops scalar/128/256/512 % | static kernels: instr / vector % / max width |", "|---|---:|---:|---:|---:|---|---|"]
    for name, _ in V:
        rr = res[(app, name)]["runs"]; m = lambda k: st.median(x[k] for x in rr); w = [st.median(x["w"][i] for x in rr) for i in range(4)]
        sf = "; ".join(f"{k}: {v['instr']} / {v['vec']}% / {v['width']}" for k, v in res[(app, name)]["static"].items())
        L.append(f"| {name} | {m('ms'):.1f} | {m('instr'):.3g} | {m('ipc'):.2f} | {m('fpi'):.2f} | {w[0]:.0f}/{w[1]:.0f}/{w[2]:.0f}/{w[3]:.0f} | {sf} |")
    L.append("")
open(os.path.join(out, "pocl_knobs.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
