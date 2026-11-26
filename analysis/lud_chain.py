#!/usr/bin/env python3
"""lud dependency-chain test: in a SCRATCH copy of the lud OpenCL build, replace the 16-term serial accumulation of lud_internal by K independent partial sums (K = 1, 2, 4, 8, 16; K = 1 is the original
kernel) and time PoCL, Intel and Intel with the vectorizer off at 1 worker thread. The production kernel is never modified.   python3 -I analysis/lud_chain.py OUTDIR
Prediction registered before the run: lud_internal's per-work-item `sum += a*b` over 16 iterations is a serial dependent chain. 2.86e9 FMAs x 4-cycle latency / 3.76 GHz = 3.0 s per repetition, close to PoCL's measured
2.7 s for lud_internal. If PoCL is latency-bound by that chain, K independent sums should shorten the chain to 16/K and PoCL's time should fall roughly as 1/K until throughput-bound (clearly below 1.5 s at K >= 4);
Intel (vectorized across work-items, i.e. already independent chains) should change little. If PoCL's time does not move, the hypothesis is wrong. Results change the floating-point summation order, so the validity flag is recorded."""
import json, os, re, shutil, statistics as st, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import apps as A
HERE = A.HERE; APP = f"{HERE}/opencl/dense-linear-algebra/lud"; out = os.path.abspath(sys.argv[1]); os.makedirs(out, exist_ok=True)
Ks = [1, 2, 4, 8, 16]; N = 2
subprocess.run(["make", "-s", "-C", APP, "build"], check=True, capture_output=True)
src = open(f"{APP}/lud_kernel.cl").read()
OLD = re.search(r"\tsum = 0;\n\tfor \(i=0; i < BLOCK_SIZE; i\+\+\)\n\t\tsum \+= m\[\(global_row_id\+get_local_id\(1\)\)\*matrix_dim\+offset\+i\] \* m\[\(offset\+i\)\*matrix_dim\+global_col_id\+get_local_id\(0\)\];\n", src)
assert OLD, "lud_internal accumulation not found in lud_kernel.cl"
def variant(K):
    if K == 1: return src
    terms = "\n".join(f"\t\ts{j} += m[(global_row_id+get_local_id(1))*matrix_dim+offset+i+{j}] * m[(offset+i+{j})*matrix_dim+global_col_id+get_local_id(0)];" for j in range(K))
    tree = "s0"
    parts = [f"s{j}" for j in range(K)]
    while len(parts) > 1: parts = [f"({parts[i]}+{parts[i+1]})" for i in range(0, len(parts), 2)]
    new = "\t" + "".join(f"float s{j} = 0;" for j in range(K)) + f"\n\tfor (i=0; i < BLOCK_SIZE; i+={K}) {{\n{terms}\n\t}}\n\tsum = {parts[0]};\n"
    return src[:OLD.start()] + new + src[OLD.end():]
RT = [("PoCL", "pocl", {}), ("Intel", "intel", {}), ("Intel vectorizer off", "intel", {"CL_CONFIG_USE_VECTORIZER": "false"})]
inp = f"{HERE}/opencl/test/dense-linear-algebra/lud/2048.dat"
def recipe(vendor):
    r = subprocess.run(["make", "-n", "-C", APP, "run", f"OCL_VENDOR={vendor}", "DWARF_REPEATS=3", "DWARF_WARMUP=1", f"RUN_ARGS=-i {inp}"], capture_output=True, text=True).stdout
    return [l for l in r.splitlines() if l.startswith("cd ")][-1]
res = {}
for K in Ks:
    sd = os.path.join(out, f"K{K}"); shutil.rmtree(sd, ignore_errors=True); os.makedirs(sd); shutil.copy(f"{APP}/build/lud", sd); open(os.path.join(sd, "lud_kernel.cl"), "w").write(variant(K))
    for name, vend, env in RT:
        cmd = recipe(vend).replace("cd build", f"cd {sd}", 1); runs = []
        for i in range(N):
            e = dict(os.environ); e.update(env); e.update(POCL_MAX_PTHREAD_COUNT="1", OMP_NUM_THREADS="1")
            if vend == "pocl": e["POCL_CACHE_DIR"] = os.path.join(sd, "pocl-cache")
            r = subprocess.run(["taskset", "-c", "0", "bash", "-c", cmd], env=e, capture_output=True, text=True); t = r.stdout + r.stderr
            m = re.search(r"median_ms=([0-9.]+)", t); runs.append(dict(ms=float(m.group(1)) if m else None, valid="valid=1" in t, text=t[-300:] if not m else ""))
        res[(K, name)] = runs
json.dump({f"{k}|{n}": v for (k, n), v in res.items()}, open(os.path.join(out, "lud_chain.json"), "w"), indent=1)
L = ["# lud dependency-chain test (1 worker thread; scratch kernel with K independent partial sums; median of 2 runs, 1 warm-up + 3 timed repetitions)", ""] + [l.strip() for l in __doc__.split("Prediction registered before the run:")[1].strip().splitlines()] + ["",
     "| K (independent sums) | chain length | PoCL ms | Intel ms | Intel vec-off ms | PoCL/Intel | valid (all runs) |", "|---|---:|---:|---:|---:|---:|---|"]
for K in Ks:
    g = lambda n: st.median(r["ms"] for r in res[(K, n)] if r["ms"]); v = all(r["valid"] for n, _, _ in RT for r in res[(K, n)])
    L.append(f"| {K} | {16//K} | {g('PoCL'):.0f} | {g('Intel'):.0f} | {g('Intel vectorizer off'):.0f} | {g('PoCL')/g('Intel'):.2f} | {v} |")
open(os.path.join(out, "lud_chain.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
