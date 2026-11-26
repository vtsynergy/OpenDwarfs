#!/usr/bin/env python3
"""Lowered-code comparison of the OpenCL CPU runtimes (PoCL vs Intel oneAPI) for one app, single worker thread.

  python3 -I analysis/cpu_lowered.py --app kmeans --out DIR [--reps 3] [--env KEY=VAL ...]

For each runtime it (1) runs the app once with the runtime's code-dump facility (PoCL: isolated kernel cache with the compiled work-group objects; Intel: CL_CONFIG_DUMP_ASM
and DUMP_IR_*), (2) extracts, per kernel, static features of the generated x86 code (instruction count, vector/scalar mix, register width, packed vs scalar FP, FMA, gathers/scatters,
mask-register use, loops, calls, stack traffic) and (3) measures dynamic counters with `perf stat` on a second, cache-warm run (user-space instructions, cycles, IPC) plus the wall/kernel
time the app reports. Output: DIR/<app>.lowered.json and DIR/<app>.lowered.md. Features are static counts of the emitted code, not execution counts; dynamic totals are per whole process."""
import argparse, glob, json, os, re, shutil, subprocess, sys, collections
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import apps as A
HERE = A.HERE; OMP = os.path.expanduser("~/.local/apps/llvm-offload-hudson/bin")
OBJDUMP = os.path.join(OMP, "llvm-objdump")

def run_app(app, vendor, extra_env, reps, warm, cwd_out, perf=False, extra_args="", events=None):
    a = A.APPS[app]
    env = dict(os.environ); env.update(extra_env); env.update(POCL_MAX_PTHREAD_COUNT="1", OMP_NUM_THREADS="1")
    cmd = ["taskset", "-c", "0", "make", "-s", "-C", f"{HERE}/opencl/{a['ocl_dir']}", "run", f"OCL_VENDOR={vendor}", f"DWARF_REPEATS={reps}", f"DWARF_WARMUP={warm}", f"RUN_ARGS={a['ocl_args']} {extra_args}".strip()]
    fifos = None
    if perf:
        import tempfile; td = tempfile.mkdtemp(prefix="perfctl"); fifos = (os.path.join(td, "ctl"), os.path.join(td, "ack")); os.mkfifo(fifos[0]); os.mkfifo(fifos[1])
        env.update(LD_PRELOAD=os.path.join(HERE, "analysis", "libmark.so"), MARK_OUT="/dev/null/mark", MARK_PERF_CTL=fifos[0], MARK_PERF_ACK=fifos[1], DWARF_WARMUP=str(warm), DWARF_REPEATS=str(reps))
        if app in A.ACC_APPS: env["MARK_ACC"] = "1"
    if perf: cmd = ["perf", "stat", "-D", "-1", "--control", f"fifo:{fifos[0]},{fifos[1]}", "-x,", "-e", events or "instructions:u,cycles:u,branches:u,branch-misses:u,fp_ret_sse_avx_ops.all:u,fp_ops_retired_by_width.scalar_uops_retired:u,fp_ops_retired_by_width.pack_128_uops_retired:u,fp_ops_retired_by_width.pack_256_uops_retired:u,fp_ops_retired_by_width.pack_512_uops_retired:u,cycles:k,instructions:k", "--"] + cmd
    r = subprocess.run(cmd, env=env, capture_output=True, text=True, cwd=cwd_out); out = r.stdout + r.stderr
    m = re.search(r"median_ms=([0-9.]+)", out); res = dict(median_ms=float(m.group(1)) if m else None, valid=("valid=1" in out))
    if perf:
        for line in out.splitlines():
            f = line.split(",")
            if len(f) > 3 and f[0].replace(".", "").isdigit() and f[2][-2:] in (":u", ":k"): res[f[2][:-2] + ("_kernel" if f[2].endswith(":k") else "")] = float(f[0])
        if res.get("cycles"): res["ipc"] = res.get("instructions", 0) / res["cycles"]
    return res

# ---------------- static feature extraction (shared by both runtimes: AT&T x86 text) ----------------
VEC_REG = re.compile(r"%([xyz])mm\d+")
def features(lines):
    f = collections.Counter(); width = collections.Counter(); labels = {}; backs = 0
    for i, l in enumerate(lines):
        m = re.match(r"^\s*([.\w$@]+):\s*$", l) or re.match(r"^[0-9a-f]+ <([^>]+)>:", l)
        if m: labels[m.group(1)] = i
    for i, l in enumerate(lines):
        am = re.match(r"^\s*([0-9a-f]+):\s", l); cur_addr = int(am.group(1), 16) if am else None
        s = l.strip()
        if not s or s.startswith(("#", ".", ";")) or s.endswith(":"): continue
        s = re.sub(r"^[0-9a-f]+:\s*", "", s)                       # objdump address prefix
        mn = s.split(None, 1)[0]; ops = s[len(mn):]
        if mn in ("nop", "nopw", "nopl", "int3", "endbr64", "xchg") or mn.startswith("data16") : continue
        f["instr"] += 1
        regs = VEC_REG.findall(ops); is_vec_op = bool(regs) and mn.startswith(("v", "p", "mov", "add", "sub", "mul", "div", "xor", "and", "or", "cvt", "sqrt", "max", "min", "cmp", "ucomis", "comis", "unpck", "shuf", "pxor"))
        if mn.startswith("v") or mn.startswith(("p", "movap", "movup", "movdq", "movs", "addp", "mulp")):
            pass
        if regs:
            w = {"x": 128, "y": 256, "z": 512}[max(regs, key=lambda c: "xyz".index(c))]
            scalar_fp = re.search(r"(ss|sd|sh)$", mn) and mn.startswith(("v", "add", "sub", "mul", "div", "mov", "cvt", "sqrt", "max", "min", "comis", "ucomis", "fmadd"))
            scalar_mv = mn in ("vmovd", "vmovq", "movd", "movq")
            if scalar_fp or scalar_mv or mn.startswith(("vfmadd", "vfmsub", "vfnmadd", "vfnmsub")) and re.search(r"s[sd]$", mn):
                f["fp_scalar"] += 1; f["vecreg_scalar_ops"] += 1
            else:
                f["vector"] += 1; width[w] += 1
                if re.search(r"p[sd]$|p[sd]\b", mn) or re.search(r"p[sd][a-z]*$", mn) and mn.startswith("v") and not mn.startswith("vp"): f["fp_packed"] += 1
                if mn.startswith("vp") or mn.startswith("p"): f["int_vector"] += 1
        if re.match(r"v?f(n)?m(add|sub)\d*[sp][sd]", mn): f["fma"] += 1
        if "gather" in mn or "scatter" in mn: f["gather_scatter"] += 1
        if re.match(r"k[a-z]+[bwdq]$", mn) or "{%k" in ops or re.search(r"%k[0-7]", ops): f["mask_ops"] += 1
        if "blend" in mn: f["blend"] += 1
        if mn.startswith(("vmov", "mov")) and "(" in ops: f["mem_ops"] += 1
        elif "(" in ops and not mn.startswith(("lea", "nop")): f["mem_ops"] += 1
        if "(%rsp)" in ops or "(%rbp)" in ops: f["stack_refs"] += 1
        if mn.startswith("call"): f["calls"] += 1
        if mn.startswith("j") and mn != "jmp" or mn == "jmp":
            f["branches"] += 1
            t = re.search(r"([.\w$@]+)\s*$", ops.strip());
            if t and t.group(1) in labels and labels[t.group(1)] <= i: backs += 1
            else:
                ta = re.search(r"(?:0x)?([0-9a-f]{2,})\s*(?:<[^>]*>)?\s*$", ops.strip())
                if cur_addr is not None and ta and mn.startswith("j") and int(ta.group(1), 16) <= cur_addr: backs += 1
    f["loops_backward_branches"] = backs
    d = dict(f); d["width_hist"] = dict(width); d["max_width"] = max(width) if width else 0
    tot = max(d.get("instr", 0), 1); d["vector_pct"] = 100.0 * d.get("vector", 0) / tot
    return d

def split_functions(text):
    """AT&T asm text -> {function: [lines]} using '.type name,@function' / '<name>:' headers"""
    funcs = collections.OrderedDict(); cur = None
    for l in text.splitlines():
        m = re.match(r"^[0-9a-f]+ <([^>]+)>:", l) or re.match(r"^([A-Za-z_][\w$.@]*):\s*$", l)
        if m and not m.group(1).startswith((".L", "LBB", ".LBB")): cur = m.group(1); funcs[cur] = []; continue
        if cur: funcs[cur].append(l)
    return funcs

def pocl_kernels(cache):
    out = {}
    for p in glob.glob(os.path.join(cache, "*", "*", "*", "*", "*.so.o")):
        kern = os.path.basename(p)[:-5]; wg = os.path.basename(os.path.dirname(p))
        txt = subprocess.run([OBJDUMP, "-d", "--no-show-raw-insn", p], capture_output=True, text=True).stdout
        funcs = split_functions(txt); lines = [l for fl in funcs.values() for l in fl]
        out[kern] = dict(file=p, wg=wg, features=features(lines), functions={k: len(v) for k, v in funcs.items()})
    return out
def intel_kernels(dumpdir, names):
    out = {}
    for p in sorted(glob.glob(os.path.join(dumpdir, "*.s"))):
        funcs = split_functions(open(p).read())
        for k, v in funcs.items():
            if k in names: out[k] = dict(file=p, wg="runtime", features=features(v), functions={k: len(v)})
    return out

def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--app", required=True); ap.add_argument("--out", required=True); ap.add_argument("--reps", type=int, default=3)
    ap.add_argument("--env", action="append", default=[], help="extra KEY=VAL for the INTEL run only (e.g. CL_CONFIG_USE_VECTORIZER=false)"); ap.add_argument("--tag", default="")
    o = ap.parse_args(); out = os.path.abspath(o.out); os.makedirs(out, exist_ok=True); app = o.app; tag = app + o.tag
    # PoCL
    pd = os.path.join(out, tag + ".pocl"); shutil.rmtree(pd, ignore_errors=True); os.makedirs(pd)
    penv = dict(POCL_CACHE_DIR=os.path.join(pd, "cache"), POCL_LEAVE_KERNEL_COMPILER_TEMP_FILES="1")
    r1 = run_app(app, "pocl", penv, 1, 0, pd); rp = run_app(app, "pocl", penv, o.reps, 1, pd, perf=True)   # 2nd run: cache warm, no JIT in the counters
    pk = pocl_kernels(penv["POCL_CACHE_DIR"])
    # Intel
    idr = os.path.join(out, tag + ".intel"); shutil.rmtree(idr, ignore_errors=True); os.makedirs(idr)
    ienv = dict(CL_CONFIG_DUMP_ASM="true", CL_CONFIG_DUMP_IR_AFTER_OPTIMIZER="true", CL_CONFIG_DUMP_FILE_NAME_PREFIX=os.path.join(idr, "k"))
    for kv in o.env: k, v = kv.split("=", 1); ienv[k] = v
    r2 = run_app(app, "intel", ienv, 1, 0, idr)
    ik = intel_kernels(idr, set(pk))
    ienv_run = {k: v for k, v in ienv.items() if not k.startswith("CL_CONFIG_DUMP")}
    ip = run_app(app, "intel", ienv_run, o.reps, 1, idr, perf=True)
    res = dict(app=app, env_intel=o.env, pocl=dict(kernels=pk, dyn=rp), intel=dict(kernels=ik, dyn=ip), kernel_names=sorted(set(pk) | set(ik)))
    json.dump(res, open(os.path.join(out, tag + ".lowered.json"), "w"), indent=1)
    L = [f"# {app}: lowered code, PoCL vs Intel OpenCL (1 worker thread{', Intel env: ' + ' '.join(o.env) if o.env else ''})", "",
         "Static features count emitted instructions (not executed ones). 'vector' = packed vector instructions; scalar FP in vector registers is counted as scalar.", "",
         "| kernel | runtime | instr | vector % | max width | packed FP | scalar FP | FMA | int-vector | gather/scatter | mask ops | blends | mem ops | backward branches | calls | stack refs |", "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for k in res["kernel_names"]:
        for name, kk in (("PoCL", pk), ("Intel", ik)):
            if k in kk:
                f = kk[k]["features"]; g = f.get
                L.append(f"| {k} | {name} | {g('instr',0)} | {f['vector_pct']:.0f} | {f['max_width']} | {g('fp_packed',0)} | {g('fp_scalar',0)} | {g('fma',0)} | {g('int_vector',0)} | {g('gather_scatter',0)} | {g('mask_ops',0)} | {g('blend',0)} | {g('mem_ops',0)} | {g('loops_backward_branches',0)} | {g('calls',0)} | {g('stack_refs',0)} |")
            else: L.append(f"| {k} | {name} | (no code found) |" + " |" * 14)
    L += ["", "## Dynamic (counters enabled only between the timer marks of the timed repetitions; 1 thread, cache-warm run)", "", "| runtime | median ms (app timer) | user instructions | user cycles | kernel-mode cycles % of all | IPC | branches | branch-miss % | FLOPs / instr | FP uops: scalar % | 128-bit % | 256-bit % | 512-bit % |", "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for name, d in (("PoCL", rp), ("Intel", ip)):
        w = {k: d.get("fp_ops_retired_by_width." + k + "_uops_retired", 0) for k in ("scalar", "pack_128", "pack_256", "pack_512")}; tw = max(sum(w.values()), 1)
        L.append(f"| {name} | {d.get('median_ms')} | {d.get('instructions',0):.3g} | {d.get('cycles',0):.3g} | {100*d.get('cycles_kernel',0)/max(d.get('cycles',0)+d.get('cycles_kernel',0),1):.0f} | {d.get('ipc',0):.2f} | {d.get('branches',0):.3g} | {100*d.get('branch-misses',0)/max(d.get('branches',1),1):.1f} | {d.get('fp_ret_sse_avx_ops.all',0)/max(d.get('instructions',1),1):.3f} | {100*w['scalar']/tw:.0f} | {100*w['pack_128']/tw:.0f} | {100*w['pack_256']/tw:.0f} | {100*w['pack_512']/tw:.0f} |")
    open(os.path.join(out, tag + ".lowered.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
if __name__ == "__main__":
    main()
