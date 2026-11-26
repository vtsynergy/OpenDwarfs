#!/usr/bin/env python3
"""Collect an evidence bundle for one (app, config) on the GPU: DWARF wall time, Nsight-derived time-budget buckets, per-kernel facts and lowered-code features.
   python3 -I analysis/gpu_collect.py --app crc --config omp-clang --out DIR [--reps 5 --warmup 1]"""
import argparse, glob, json, os, re, sqlite3, statistics, subprocess, sys, shlex, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import apps as A, lowered_gpu as L, budget as BG
NSYS = "/opt/nvidia/hpc_sdk/Linux_x86_64/25.7/compilers/bin/nsys"
NCU = "/opt/nvidia/hpc_sdk/Linux_x86_64/25.7/compilers/bin/ncu"
NCU_METRICS = "gpu__time_duration.sum,launch__block_size,launch__grid_size,launch__registers_per_thread,launch__shared_mem_per_block_static,sm__warps_active.avg.pct_of_peak_sustained_active,smsp__inst_executed.sum,smsp__thread_inst_executed_per_inst_executed.ratio,smsp__inst_executed_op_global_ld.sum,smsp__inst_executed_op_global_st.sum,smsp__inst_executed_op_shared_ld.sum,smsp__inst_executed_op_shared_st.sum,smsp__inst_executed_op_local_ld.sum,smsp__inst_executed_op_local_st.sum,dram__bytes_read.sum,dram__bytes_write.sum,sm__throughput.avg.pct_of_peak_sustained_elapsed,gpu__dram_throughput.avg.pct_of_peak_sustained_elapsed"
BUCKETS = [("launch", r"^(cu|cuda)LaunchKernel"), ("copy_h2d", r"Memcpy.*HtoD|MemcpyHtoD|cudaMemcpy.*Host.*Device"), ("copy_d2h", r"MemcpyDtoH"),
           ("copy_other", r"^cu(da)?Memcpy|^cuMemcpy"), ("sync", r"Synchronize|StreamWaitEvent|EventQuery|StreamQuery"),
           ("alloc", r"MemAlloc|cudaMalloc|MemHostAlloc|MemCreate|MemMap|MemAddressReserve"), ("free", r"MemFree|cudaFree|MemRelease|MemUnmap|MemAddressFree"),
           ("event", r"^cu(da)?Event(Record|Create|Destroy)"), ("init", r"^cu(Init|Ctx|DevicePrimary|Module|Library|Device|Get|Func|Occup|Pointer|Stream(Create|Destroy|Is)|Mem(HostRegister|HostUnregister|GetInfo))")]
def bucket(name):
    for b, rx in BUCKETS:
        if re.search(rx, name): return b
    return "other"
def build_cmd(app, cfg, reps, warm):
    a = A.APPS[app]; e = dict(CUDA_VISIBLE_DEVICES="0", OMP_THREAD_LIMIT="2048", OMP_NUM_THREADS=str(A.NTHR))
    if cfg.startswith("omp-"):
        e["OMP_TARGET_OFFLOAD"] = "MANDATORY"
        cmd = ["make", "-s", "-C", f"{A.HERE}/omp/{a['omp_dir']}", "run", "DEVICE=h100", f"TOOLCHAIN={cfg[4:]}", f"RUN_REPEATS={reps}", f"RUN_WARMUP={warm}"] + a["omp_args"]
    else:
        cmd = ["make", "-s", "-C", f"{A.HERE}/opencl/{a['ocl_dir']}", "run", "OCL_VENDOR=nvidia", "GPU_INDEX=0", f"DWARF_REPEATS={reps}", f"DWARF_WARMUP={warm}", f"RUN_ARGS={a['ocl_args']}"]
    return e, cmd
def q(db, sql):
    try: return list(db.execute(sql))
    except sqlite3.OperationalError: return []
def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--app", required=True); ap.add_argument("--config", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--reps", type=int, default=5); ap.add_argument("--warmup", type=int, default=1); ap.add_argument("--no-capture", action="store_true", help="profile the whole process instead of only the timed repetitions"); ap.add_argument("--ncu", action="store_true", help="also collect per-kernel Nsight Compute counters (1 profiled launch per kernel)")
    o = ap.parse_args(); os.makedirs(o.out, exist_ok=True); tag = f"{o.app}.{o.config}"
    env, cmd = build_cmd(o.app, o.config, o.reps, o.warmup)
    full = dict(os.environ); full.update(env)
    cap = not o.no_capture
    if cap:
        full.update(LD_PRELOAD=os.path.join(A.HERE, "analysis", "libmark.so"), MARK_PROFILER="1", MARK_OUT="/dev/null/mark")
        if o.app in A.ACC_APPS: full["MARK_ACC"] = "1"
    rep = os.path.join(o.out, tag)
    pre = ["taskset", "-c", A.CPUSET, NSYS, "profile", "-f", "true", "--trace=cuda", "--sample=none", "--cpuctxsw=none", "-o", rep, "--stats=false"] + (["--capture-range=cudaProfilerApi", "--capture-range-end=stop"] if cap else [])
    t0 = time.time(); r = subprocess.run(pre + cmd, env=full, capture_output=True, text=True); log = r.stdout + r.stderr
    open(rep + ".log", "w").write(log)
    m = re.search(r"DWARF_RESULT,.*", log)
    if not m: sys.exit(f"no DWARF_RESULT for {tag}; see {rep}.log")
    res = dict(kv.split("=", 1) for kv in m.group(0).split(",")[1:] if "=" in kv)
    subprocess.run([NSYS, "export", "--type", "sqlite", "-f", "true", "-o", rep + ".sqlite", rep + ".nsys-rep"], capture_output=True)
    db = sqlite3.connect(rep + ".sqlite"); nreps = o.reps if cap else o.reps + o.warmup
    S = dict(q(db, "select id,value from StringIds"))
    api = {}; bk = {}; ivs = {}
    for st, en, nid in q(db, "select start,end,nameId from CUPTI_ACTIVITY_KIND_RUNTIME"):
        n = S.get(nid, "?"); n = re.sub(r"_v\d+$|_ptsz$|_ptds$", "", n); b = bucket(n)
        ivs.setdefault(b, []).append((st, en)); a = api.setdefault(n, dict(bucket=b, count=0, total_ns=0)); a["count"] += 1; a["total_ns"] += en - st
        bk.setdefault(b, dict(count=0, total_ns=0)); bk[b]["count"] += 1; bk[b]["total_ns"] += en - st
    kern = {}
    for st, en, dn, mn, sn, rg, gx, gy, gz, bx, by, bz, ssm, dsm, lm in q(db, "select start,end,demangledName,mangledName,shortName,registersPerThread,gridX,gridY,gridZ,blockX,blockY,blockZ,staticSharedMemory,dynamicSharedMemory,localMemoryPerThread from CUPTI_ACTIVITY_KIND_KERNEL"):
        raw = S.get(mn) or S.get(sn) or S.get(dn) or "?"; key = L.kernel_key(raw)
        k = kern.setdefault(key, dict(raw_name=raw, instances=0, total_ns=0, regs=rg, grids=set(), blocks=set(), static_smem=ssm, local_per_thread=lm))
        k["instances"] += 1; k["total_ns"] += en - st; k["grids"].add((gx, gy, gz)); k["blocks"].add((bx, by, bz))
    for k in kern.values(): k["grids"] = sorted(k["grids"])[:4]; k["blocks"] = sorted(k["blocks"])[:4]
    mem = {}
    for st, en, nb, ck in q(db, "select start,end,bytes,copyKind from CUPTI_ACTIVITY_KIND_MEMCPY"):
        d = mem.setdefault(ck, dict(count=0, total_ns=0, bytes=0)); d["count"] += 1; d["total_ns"] += en - st; d["bytes"] += nb
    kiv = [(r[0], r[1]) for r in q(db, "select start,end from CUPTI_ACTIVITY_KIND_KERNEL")]; civ = [(r[0], r[1]) for r in q(db, "select start,end from CUPTI_ACTIVITY_KIND_MEMCPY")]
    part = BG.partition(ivs, kiv, civ)
    lowered = {}
    if o.config.startswith("omp-"):
        suf = "_h100" if o.config == "omp-clang" else "_h100_nvhpc"
        a = A.APPS[o.app]; bins = glob.glob(f"{A.HERE}/omp/{a['omp_dir']}/build/*{suf}")
        if bins: lowered = L.extract(bins[0])
    per = lambda ns: ns / 1e6 / nreps
    bundle = dict(app=o.app, config=o.config, date=time.strftime("%F %T"), reps=o.reps, warmup=o.warmup, wall_ms=float(res["median_ms"]),
                  buckets_ms_per_rep={b: per(v["total_ns"]) for b, v in bk.items()}, bucket_counts_per_rep={b: v["count"] / nreps for b, v in bk.items()},
                  kernel_ms_per_rep=per(sum(k["total_ns"] for k in kern.values())), launches_per_rep=sum(k["instances"] for k in kern.values()) / nreps,
                  api={n: dict(a, total_ms_per_rep=per(a["total_ns"]), count_per_rep=a["count"] / nreps) for n, a in api.items()},
                  kernels={k: dict(v, ms_per_rep=per(v["total_ns"]), inst_per_rep=v["instances"] / nreps) for k, v in kern.items()},
                  memcpy={str(k): dict(v, ms_per_rep=per(v["total_ns"]), count_per_rep=v["count"] / nreps) for k, v in mem.items()},
                  lowered=lowered, toolchain=res.get("toolchain", ""), valid=res.get("valid", ""),
                  partition_ms_per_rep={k: v / 1e6 / nreps for k, v in part["buckets_ns"].items()}, window_ms_per_rep=part["window_ns"] / 1e6 / nreps, init_in_window_ms=part["init_in_window_ns"] / 1e6)
    if o.ncu:
        import csv, io
        bundle["ncu"] = {}
        for key, kk in kern.items():
            skip = min(3, max(kk["instances"] - 1, 0))
            r2 = subprocess.run(["taskset", "-c", A.CPUSET, NCU, "--kernel-name", "regex:" + re.escape(key), "--launch-skip", str(skip), "--launch-count", "1", "--metrics", NCU_METRICS, "--csv"] + cmd,
                                env=full, capture_output=True, text=True, timeout=600)
            rows = list(csv.DictReader(io.StringIO("\n".join(l for l in r2.stdout.splitlines() if l.startswith('"')))))
            m = {}
            for x in rows:
                try: m[x["Metric Name"]] = float(x["Metric Value"].replace(",", ""))
                except Exception: pass
            bundle["ncu"][key] = m
    json.dump(bundle, open(rep + ".bundle.json", "w"), indent=1)
    print(f"{tag}: wall={bundle['wall_ms']:.2f} ms  kernel={bundle['kernel_ms_per_rep']:.2f} ms/rep  launches/rep={bundle['launches_per_rep']:.0f}  ({time.time()-t0:.0f}s) -> {rep}.bundle.json")
main()
