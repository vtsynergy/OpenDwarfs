#!/usr/bin/env python3
"""Compare two evidence bundles (same app, two toolchains) -> gap decomposition + ranked, evidence-backed findings (a rule engine).
   python3 -I analysis/compare.py BASE.bundle.json OTHER.bundle.json [--md out.md]
   Convention: gap = BASE wall - OTHER wall (positive = OTHER is faster). Every finding states its evidence and the intervention that would confirm it;
   a finding is a *hypothesis* until that intervention has been run (see analysis/README.md)."""
import argparse, json, sys
HOST_B = ["launch", "copy_h2d", "copy_d2h", "copy_other", "sync", "alloc", "free", "event", "other"]
def load(p): return json.load(open(p))
def g(d, *k, default=0.0):
    for x in k:
        if not isinstance(d, dict) or x not in d: return default
        d = d[x]
    return d
def compare(A, B):
    F = []; bA, bB = A["buckets_ms_per_rep"], B["buckets_ms_per_rep"]; cA, cB = A["bucket_counts_per_rep"], B["bucket_counts_per_rep"]
    gap = A["wall_ms"] - B["wall_ms"]; nl = max(A["launches_per_rep"], B["launches_per_rep"], 1)
    delta = {b: bA.get(b, 0) - bB.get(b, 0) for b in HOST_B}
    pA, pB = A.get("partition_ms_per_rep"), B.get("partition_ms_per_rep"); additive = bool(pA and pB)
    if additive:   # additive accounting: every instant of the window belongs to exactly one state (analysis/budget.py)
        delta = {b: pA.get(b, 0) - pB.get(b, 0) for b in set(pA) | set(pB)}
        for b in HOST_B: delta.setdefault(b, 0.0)
    kA, kB = A["kernel_ms_per_rep"], B["kernel_ms_per_rep"]
    # --- launch-count / structure
    if abs(A["launches_per_rep"] - B["launches_per_rep"]) > 0.1 * nl:
        F.append(dict(id="L1", title="different number of kernel launches per repetition", ms=None, conf="measured",
                      evidence=f"{A['launches_per_rep']:.0f} vs {B['launches_per_rep']:.0f} launches/rep: the two builds do not decompose the work the same way",
                      check="compare source-level launch structure / LAUNCH-CONFIGS.md"))
    # --- kernel time
    kr = kA / kB if kB else float("inf")
    if abs(kr - 1) <= 0.10:
        F.append(dict(id="K0", title="device kernel time is not the cause", ms=0.0, conf="measured",
                      evidence=f"summed kernel time per rep {kA:.2f} vs {kB:.2f} ms (ratio {kr:.2f}); the {abs(gap):.1f} ms wall gap is outside the kernels", check="-"))
    else:
        rows = []
        for k in sorted(set(A["kernels"]) & set(B["kernels"])):
            a, b = A["kernels"][k], B["kernels"][k]
            rows.append((k, a["ms_per_rep"], b["ms_per_rep"], a["regs"], b["regs"], a["blocks"][:1], b["blocks"][:1], a["grids"][:1], b["grids"][:1]))
        ev = "; ".join(f"{k}: {x:.2f} vs {y:.2f} ms (regs {ra} vs {rb}, block {ba} vs {bb}, grid {ga} vs {gb})" for k, x, y, ra, rb, ba, bb, ga, gb in rows)
        F.append(dict(id="K1", title="device kernel time differs", ms=(delta["device_kernel"] if additive else kA - kB), conf="measured", evidence=f"summed kernel time {kA:.2f} vs {kB:.2f} ms/rep: " + ev, check="compare ncu counters per kernel (below); pin identical launch config and re-measure"))
    # --- launch configuration mismatch (per matched kernel)
    for k in sorted(set(A["kernels"]) & set(B["kernels"])):
        a, b = A["kernels"][k], B["kernels"][k]
        if a["blocks"] != b["blocks"] or a["grids"] != b["grids"]:
            F.append(dict(id="C1", title=f"launch configuration differs for kernel '{k}'", ms=None, conf="measured",
                          evidence=f"grid {a['grids'][:2]} vs {b['grids'][:2]}, block {a['blocks'][:2]} vs {b['blocks'][:2]}", check="pin num_teams/thread_limit and confirm in the trace that both launch the same shape"))
    # --- allocation churn
    apl, bpl = (cA.get("alloc", 0) + cA.get("free", 0)) / nl, (cB.get("alloc", 0) + cB.get("free", 0)) / nl
    if max(apl, bpl) >= 0.5 and min(apl, bpl) < 0.05:
        lo = "BASE" if apl < bpl else "OTHER"
        F.append(dict(id="A1", title="device allocation/free around every launch in one toolchain only", ms=delta["alloc"] + delta["free"], conf="measured",
                      evidence=f"alloc+free per launch: {apl:.2f} (BASE) vs {bpl:.2f} (OTHER); time {bA.get('alloc',0)+bA.get('free',0):.2f} vs {bB.get('alloc',0)+bB.get('free',0):.2f} ms/rep",
                      check="hoist the mapping (target enter data) or pool allocations; re-measure the wall gap"))
    # --- event / sync overhead
    if max(cA.get("event", 0), cB.get("event", 0)) / nl >= 1 and min(cA.get("event", 0), cB.get("event", 0)) / nl < 0.1:
        F.append(dict(id="E1", title="CUDA event bookkeeping around every launch in one toolchain only", ms=delta["event"], conf="measured",
                      evidence=f"event calls per launch {cA.get('event',0)/nl:.1f} vs {cB.get('event',0)/nl:.1f}; time {bA.get('event',0):.2f} vs {bB.get('event',0):.2f} ms/rep", check="runtime knob / different mapping path"))
    sA, sB = bA.get("sync", 0) / nl * 1000, bB.get("sync", 0) / nl * 1000
    if max(sA, sB) > 5 and max(sA, sB) / max(min(sA, sB), 0.1) >= 3:
        F.append(dict(id="S1", title="much more host time blocked in synchronization in one toolchain", ms=delta["sync"], conf="measured",
                      evidence=f"sync {sA:.1f} vs {sB:.1f} us/launch ({bA.get('sync',0):.2f} vs {bB.get('sync',0):.2f} ms/rep)", check="count and placement of stream syncs; async copies"))
    # --- copy path
    def d2h(X): return g(X, "api", "cuMemcpyDtoHAsync", "total_ms_per_rep") * 1000 / max(g(X, "api", "cuMemcpyDtoHAsync", "count_per_rep"), 1e-9)
    dA, dB = d2h(A), d2h(B)
    devcopy = lambda X: sum(v["ms_per_rep"] for kk, v in X["memcpy"].items() if kk == "2") / max(sum(v["count_per_rep"] for kk, v in X["memcpy"].items() if kk == "2"), 1e-9) * 1000
    if dA and dB and max(dA, dB) / min(dA, dB) >= 1.5:
        pin = lambda X: any("AllocHost" in n or "HostAlloc" in n for n in X["api"])
        F.append(dict(id="D1", title="device-to-host copy call is slower in one toolchain although the device-side copy time is similar", ms=delta["copy_d2h"], conf="measured; mechanism hypothesis",
                      evidence=f"D2H API call {dA:.0f} vs {dB:.0f} us avg; device-side memcpy {devcopy(A):.1f} vs {devcopy(B):.1f} us; pinned host allocation seen: BASE={pin(A)} OTHER={pin(B)} (pinned staging would explain faster async copies)",
                      check="allocate pinned host staging / compare cudaMemcpyAsync source memory kind in trace"))
    # --- GPU utilisation
    ua, ub = kA / A["wall_ms"], kB / B["wall_ms"]
    F.append(dict(id="U1", title="device busy fraction of the timed region", ms=None, conf="measured",
                  evidence=f"kernel time / wall = {ua:.0%} (BASE) vs {ub:.0%} (OTHER); {A['launches_per_rep']:.0f} launches/rep of {kA/max(A['launches_per_rep'],1)*1000:.1f} us (BASE) vs {kB/max(B['launches_per_rep'],1)*1000:.1f} us (OTHER)", check="-"))
    # --- OpenCL-runtime rules (bundles from ocl_collect.py carry idle_ms_per_rep)
    if "idle_ms_per_rep" in A and "idle_ms_per_rep" in B:
        ia, ib = A["idle_ms_per_rep"] / nl * 1000, B["idle_ms_per_rep"] / nl * 1000
        if max(ia, ib) > 3 and max(ia, ib) / max(min(ia, ib), 0.5) >= 2:
            F.append(dict(id="Q1", title="device idle time between consecutive launches differs", ms=A["idle_ms_per_rep"] - B["idle_ms_per_rep"], conf="measured (event profiling)",
                          evidence=f"{ia:.1f} vs {ib:.1f} us idle per launch ({A['idle_ms_per_rep']:.2f} vs {B['idle_ms_per_rep']:.2f} ms/rep): time the device waits between one command finishing and the next kernel starting (host work, transfers, runtime scheduling)",
                          check="vary worker-thread count / batch launches; compare against host time between launches"))
    la, lb = bA.get("launch", 0) / nl * 1000, bB.get("launch", 0) / nl * 1000
    if max(la, lb) > 3 and max(la, lb) / max(min(la, lb), 0.5) >= 2:
        F.append(dict(id="H1", title="host time spent in launch calls (enqueue + argument setting) differs per launch", ms=delta["launch"], conf="measured",
                      evidence=f"{la:.1f} vs {lb:.1f} us per launch ({bA.get('launch',0):.2f} vs {bB.get('launch',0):.2f} ms/rep)", check="batch arguments / reduce launches; compare API call counts per launch"))
    wa, wb = A["wall_ms"] / max(kA, 1e-9), B["wall_ms"] / max(kB, 1e-9)
    F.append(dict(id="O1", title="overhead factor: wall time / device kernel time", ms=None, conf="measured",
                  evidence=f"{wa:.1f}x (BASE) vs {wb:.1f}x (OTHER); non-kernel time per rep {A['wall_ms']-kA:.2f} vs {B['wall_ms']-kB:.2f} ms", check="-"))
    # --- ncu drill-down
    nc = []
    for k in sorted(set(A.get("ncu", {})) & set(B.get("ncu", {}))):
        a, b = A["ncu"][k], B["ncu"][k]
        keys = ["gpu__time_duration.sum", "launch__block_size", "launch__grid_size", "launch__registers_per_thread", "sm__warps_active.avg.pct_of_peak_sustained_active", "smsp__inst_executed.sum",
                "smsp__inst_executed_op_global_ld.sum", "smsp__inst_executed_op_local_ld.sum", "dram__bytes_read.sum", "sm__throughput.avg.pct_of_peak_sustained_elapsed"]
        nc.append((k, [(m, a.get(m), b.get(m)) for m in keys]))
        ir = a.get("smsp__inst_executed.sum", 0) / max(b.get("smsp__inst_executed.sum", 1), 1)
        if abs(ir - 1) > 0.15:
            F.append(dict(id="K2", title=f"executed-instruction count differs for kernel '{k}' ({ir:.2f}x)", ms=None, conf="measured (1 profiled launch)",
                          evidence=f"inst executed {a.get('smsp__inst_executed.sum')} vs {b.get('smsp__inst_executed.sum')}; block {a.get('launch__block_size')} vs {b.get('launch__block_size')}; regs {a.get('launch__registers_per_thread')} vs {b.get('launch__registers_per_thread')}",
                          check="compare SASS; rule out launch-shape effects (different block size changes executed instruction counts)"))
    if additive:
        wgap = A["window_ms_per_rep"] - B["window_ms_per_rep"]
        F.append(dict(id="HO1", title="host time outside traced runtime calls (device idle) differs", ms=delta.get("host_other", 0.0), conf="measured (residual state of the partition)",
                      evidence=f"{pA.get('host_other',0):.2f} vs {pB.get('host_other',0):.2f} ms/rep; includes untimed host code that is the same in both runs unless the toolchain changes host-side work", check="-"))
        F.append(dict(id="W1", title="accounting check", ms=None, conf="measured",
                      evidence=f"sum of bucket deltas {sum(delta.values()):+.2f} ms = window gap {wgap:+.2f} ms (exact by construction); measured wall gap {gap:+.2f} ms; window/wall coverage {A['window_ms_per_rep']/A['wall_ms']:.2f} (BASE) vs {B['window_ms_per_rep']/B['wall_ms']:.2f} (OTHER)", check="coverage far from 1 means untimed host code is inside the window"))
    explained = sum(f["ms"] for f in F if f.get("ms") and f["id"] in ("A1", "E1", "S1", "D1", "K1", "H1") and f["ms"] > 0)
    return dict(gap_ms=gap, delta_ms=delta, kernel_ms=(kA, kB), findings=F, explained_ms=explained, ncu=nc, additive=additive)
def render(A, B, R):
    L = [f"# {A['app']}: {A['config']} (BASE) vs {B['config']} (OTHER)", "",
         f"Wall per rep: {A['wall_ms']:.2f} vs {B['wall_ms']:.2f} ms -> gap {R['gap_ms']:+.2f} ms (BASE - OTHER); kernel time/rep {R['kernel_ms'][0]:.2f} vs {R['kernel_ms'][1]:.2f} ms; launches/rep {A['launches_per_rep']:.0f} vs {B['launches_per_rep']:.0f}.", "",
         "## Host-timeline buckets (ms per repetition, from the Nsight CUDA API trace)", "", "| bucket | BASE | OTHER | delta | calls/rep BASE | OTHER |", "|---|---:|---:|---:|---:|---:|"]
    if R.get("additive"):
        L[-4:] = ["## Additive time budget (ms per repetition; every instant of the steady-state window is in exactly one state, device-busy first)", "",
                  "| state | BASE | OTHER | delta | share of window gap |", "|---|---:|---:|---:|---:|"]
        pa, pb = A["partition_ms_per_rep"], B["partition_ms_per_rep"]; wg = A["window_ms_per_rep"] - B["window_ms_per_rep"]
        for b in ["device_kernel", "device_copy", "sync", "copy_d2h", "copy_h2d", "copy_other", "alloc", "free", "event", "launch", "other", "host_other"]:
            x, y = pa.get(b, 0), pb.get(b, 0)
            if x or y: L.append(f"| {b} | {x:.2f} | {y:.2f} | {x-y:+.2f} | {(x-y)/wg*100 if wg else 0:+.0f}% |")
        L.append(f"| **window (sum)** | {A['window_ms_per_rep']:.2f} | {B['window_ms_per_rep']:.2f} | {wg:+.2f} | 100% |"); L += ["", "(The old per-API-class table is not additive: host waits include the kernel time they wait for.)", ""]
        return_early = True
    else: return_early = False
    for b in ([] if return_early else HOST_B):
        x, y = A["buckets_ms_per_rep"].get(b, 0), B["buckets_ms_per_rep"].get(b, 0)
        if x or y: L.append(f"| {b} | {x:.2f} | {y:.2f} | {x-y:+.2f} | {A['bucket_counts_per_rep'].get(b,0):.0f} | {B['bucket_counts_per_rep'].get(b,0):.0f} |")
    tot = sum(A["buckets_ms_per_rep"].get(b, 0) for b in HOST_B) - sum(B["buckets_ms_per_rep"].get(b, 0) for b in HOST_B)
    if return_early: L += ["## Findings (hypotheses until their check has been run; contributions are additive deltas of the budget above)", ""]
    else: L += ["", f"Sum of host API bucket deltas: {tot:+.2f} ms of the {R['gap_ms']:+.2f} ms gap (the rest is host time outside the traced calls). NOTE: sync/copy buckets include time spent waiting for kernels, so bucket deltas overlap with kernel time and are not additive.", "", "## Findings (hypotheses until their check has been run)", ""]
    for f in R["findings"]:
        L.append(f"- **{f['id']} {f['title']}** [{f['conf']}]" + (f" - contributes {f['ms']:+.2f} ms/rep" if f.get("ms") not in (None,) else "")); L.append(f"  - evidence: {f['evidence']}"); L.append(f"  - confirm by: {f['check']}")
    if R["ncu"]:
        L += ["", "## Per-kernel Nsight Compute counters (one profiled launch each)", ""]
        for k, rows in R["ncu"]:
            L += [f"**{k}**", "", "| metric | BASE | OTHER |", "|---|---:|---:|"] + [f"| {m} | {a} | {b} |" for m, a, b in rows] + [""]
    return "\n".join(L)
if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("base"); ap.add_argument("other"); ap.add_argument("--md")
    o = ap.parse_args(); A, B = load(o.base), load(o.other); R = compare(A, B); md = render(A, B, R)
    if o.md: open(o.md, "w").write(md + "\n"); json.dump(R, open(o.md.replace(".md", ".findings.json"), "w"), indent=1, default=str)
    print(md)
