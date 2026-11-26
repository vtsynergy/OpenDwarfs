#!/usr/bin/env python3
"""Memory-behaviour comparison for lud, PoCL vs Intel (and PoCL variants), 1 worker thread, AMD Zen 4 counters over the timed repetitions only.
  python3 -I analysis/lud_memory.py OUTDIR
Prediction registered before the run: if Intel's lud advantage (590 ms vs 2820 ms per repetition) comes from better memory behaviour, Intel shows fewer L1D fills / loads per FLOP and a smaller
'load not complete' stall share than PoCL; PoCL's variants (forced vector widths, vectorizer off) change the instruction count by up to 6x but time by < 7%, so their memory counters should stay about equal
(explaining the flat time). If the memory counters of PoCL and Intel are about equal, memory behaviour is NOT the cause."""
import json, os, statistics as st, subprocess, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cpu_lowered as C, apps as A
out = os.path.abspath(sys.argv[1]); os.makedirs(out, exist_ok=True); N = 2; REPS = 3
G = {"g1": "cycles:u,instructions:u,ls_dispatch.ld_dispatch:u,ls_dispatch.store_dispatch:u,ls_any_fills_from_sys.all:u,ls_stlf:u,ls_bad_status2.stli_other:u",
     "g2": "ls_any_fills_from_sys.local_l2:u,ls_any_fills_from_sys.local_ccx:u,ls_any_fills_from_sys.dram_io_near:u,ls_any_fills_from_sys.dram_io_far:u,ls_hw_pf_dc_fills.all:u,ls_dmnd_fills_from_sys.all:u",
     "g3": "ex_no_retire.load_not_complete:u,ex_no_retire.not_complete:u,ex_no_retire.empty:u,ls_l1_d_tlb_miss.all:u,ls_misal_loads.ma64:u,cycles:u",
     "g4": "fp_ret_sse_avx_ops.all:u,cycles:u,instructions:u"}
V = [("PoCL", "pocl", {}), ("PoCL (repeat)", "pocl", {}), ("PoCL VF16+512", "pocl", {"POCL_VECTORIZER_FORCE_VECTOR_WIDTH": "16", "POCL_VECTORIZER_PREFER_VECTOR_WIDTH": "512"}),
     ("PoCL vectorizer off", "pocl", {"POCL_VECTORIZER_FORCE_VECTOR_WIDTH": "1"}), ("Intel", "intel", {}), ("Intel vectorizer off", "intel", {"CL_CONFIG_USE_VECTORIZER": "false"})]
res = {}
for name, vend, env in V:
    cache = os.path.join(out, "cache-" + name.replace(" ", "_").replace("+", "_"))
    e = dict(env)
    if vend == "pocl": e["POCL_CACHE_DIR"] = cache
    C.run_app("lud", vend, e, 1, 0, out)                       # populate caches
    ms = []; ev = {}
    for gname, g in G.items():
        runs = [C.run_app("lud", vend, e, REPS, 1, out, perf=True, events=g) for _ in range(N)]
        ms += [r["median_ms"] for r in runs]
        for k in {k for r in runs for k in r if k not in ("median_ms", "valid", "ipc")}:
            ev.setdefault(k, []).append(st.median(r.get(k, 0) for r in runs))
    res[name] = dict(ms=st.median(ms), ev={k: st.median(v) for k, v in ev.items()})
json.dump(res, open(os.path.join(out, "lud_memory.json"), "w"), indent=1)
# per-kernel event-timed kernel time (interposer, 1 thread)
pk = {}
for cfg in ("ocl-pocl", "ocl-intel"):
    subprocess.run([sys.executable, "-I", os.path.join(A.HERE, "analysis", "ocl_collect.py"), "--app", "lud", "--config", cfg, "--threads", "1", "--tag", ".t1", "--out", out], capture_output=True, text=True, env=dict(os.environ))
    b = json.load(open(os.path.join(out, f"lud.{cfg}.t1.bundle.json"))); pk[cfg] = {k: (v["ms_per_rep"], v["inst_per_rep"]) for k, v in b["kernels"].items()}; pk[cfg]["_wall"] = (b["wall_ms"], 0)
FL = 2.0 / 3.0 * 2048 ** 3                                   # useful FLOPs of one LU factorisation (2048 x 2048), one repetition
L = ["# lud memory behaviour: PoCL vs Intel (1 worker thread, Zen 4 counters over 3 timed repetitions, median of 2 collections per counter group)", ""] + [l.strip() for l in __doc__.split("Prediction registered before the run:")[1].strip().splitlines()] + ["",
     "## Event-timed kernel time per repetition (interposer, 1 thread)", "", "| kernel | PoCL ms | Intel ms | PoCL/Intel | launches/rep |", "|---|---:|---:|---:|---:|"]
for k in sorted(set(pk["ocl-pocl"]) | set(pk["ocl-intel"])):
    a = pk["ocl-pocl"].get(k, (0, 0)); b = pk["ocl-intel"].get(k, (0, 0)); L.append(f"| {k} | {a[0]:.1f} | {b[0]:.1f} | {a[0]/max(b[0],1e-9):.2f} | {a[1]:.0f} |")
L += ["", f"## Counters per repetition (totals / {REPS}); FLOPs per repetition = 2/3 n^3 = {FL:.3g}", "",
      "| variant | ms | instr | loads | stores | L1D fills | fills/FLOP x1e3 | fills from L2 % | L3 % | DRAM % | HW-prefetch fills % | load-not-complete stall cycles % | no-retire empty % | L1 DTLB miss | store-fwd | STLI stall | misaligned (64B) loads |", "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
for name, _, _ in V:
    r = res[name]; e = r["ev"]; g = lambda k: e.get(k, 0) / REPS; fills = max(g("ls_any_fills_from_sys.all"), 1); cyc = max(e.get("cycles", 1) / REPS, 1)
    sh = lambda k: 100 * g(k) / fills
    L.append(f"| {name} | {r['ms']:.0f} | {g('instructions'):.3g} | {g('ls_dispatch.ld_dispatch'):.3g} | {g('ls_dispatch.store_dispatch'):.3g} | {fills:.3g} | {fills/FL*1000:.2f} | {sh('ls_any_fills_from_sys.local_l2'):.0f} | {sh('ls_any_fills_from_sys.local_ccx'):.0f} | {sh('ls_any_fills_from_sys.dram_io_near')+sh('ls_any_fills_from_sys.dram_io_far'):.0f} | {sh('ls_hw_pf_dc_fills.all'):.0f} | {100*g('ex_no_retire.load_not_complete')/cyc:.0f} | {100*g('ex_no_retire.empty')/cyc:.0f} | {g('ls_l1_d_tlb_miss.all'):.3g} | {g('ls_stlf'):.3g} | {g('ls_bad_status2.stli_other'):.3g} | {g('ls_misal_loads.ma64'):.3g} |")
open(os.path.join(out, "lud_memory.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
