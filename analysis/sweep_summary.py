#!/usr/bin/env python3
"""Summarise the CPU thread sweep (PROTOCOL-final.md Addendum 3): reads omp/results/perf-sweep-t<N>/raw.csv and writes sweep-summary.md / sweep-results.csv.
   python3 -I analysis/sweep_summary.py [RESULTS_DIR]   (default: omp/results)
Analysis as pre-registered: (a) speedup vs omp-clang at the same thread count; (b) each configuration's best-thread cell vs its own 128-thread cell, and the best-thread speedup vs
omp-clang at ITS best thread count; (c) apps where a configuration's best thread count is below 128."""
import csv, glob, math, os, re, statistics as st, sys, collections
d = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "omp", "results")
APPS = "spmv bfs crc nw srad cfd kmeans nqueens hmm swat tdm gem fft lud".split()
CFGS = ["omp-clang", "omp-nvhpc", "opencl", "opencl-intel"]; NAME = {"omp-clang": "omp-clang", "omp-nvhpc": "omp-nvhpc", "opencl": "PoCL", "opencl-intel": "Intel OpenCL"}
cell = collections.defaultdict(list); status = collections.defaultdict(collections.Counter); T = set()
for f in glob.glob(os.path.join(d, "perf-sweep-t*", "raw.csv")):
    t = int(re.search(r"perf-sweep-t(\d+)", f).group(1)); T.add(t)
    for r in csv.DictReader(open(f)):
        k = (r["app"], r["config"], t); status[k][r["status"]] += 1
        if r["status"] == "OK": cell[k].append(float(r["median_ms"]))
T = sorted(T); ms = {k: st.median(v) for k, v in cell.items()}
gm = lambda v: math.exp(sum(map(math.log, v)) / len(v)) if v else float("nan")
fmt = lambda v: "-" if v is None else (f"{v:.3f}" if v < 1 else f"{v:.1f}" if v < 100 else f"{v:.0f}")
L = [f"# CPU thread sweep (threads {T}); median of launch medians, ms", "", "Exploratory; 3 launches x (1 warm-up + 5 reps) per cell; incomplete cells shown as '-' (see status counts in raw.csv).", ""]
rows = [["app", "config", "threads", "median_ms", "speedup_vs_omp_clang_same_threads", "status"]]
for c in CFGS:
    L += [f"## {NAME[c]}: ms by thread count (best marked *)", "", "| app | " + " | ".join(f"t={t}" for t in T) + " | best t | best / t=128 |", "|---|" + "---:|" * (len(T) + 2)]
    for a in APPS:
        v = {t: ms.get((a, c, t)) for t in T}; ok = {t: x for t, x in v.items() if x}
        bt = min(ok, key=ok.get) if ok else None
        r128 = (ok[bt] / ok[128]) if bt and 128 in ok else None
        L.append(f"| {a} | " + " | ".join((fmt(v[t]) + ("*" if t == bt else "")) for t in T) + f" | {bt if bt else '-'} | {('%.2f' % r128) if r128 else '-'} |")
        for t in T:
            base = ms.get((a, "omp-clang", t)); x = v[t]
            rows.append([a, c, t, f"{x:.6f}" if x else "", f"{base/x:.4f}" if (x and base) else "", ",".join(f"{k}x{n}" for k, n in sorted(status.get((a, c, t), {}).items())) or "NOT_RUN"])
    L.append("")
L += ["## Speedup vs omp-clang at the same thread count (geomean over apps with both cells)", "", "| config | " + " | ".join(f"t={t}" for t in T) + " |", "|---|" + "---:|" * len(T)]
for c in CFGS[1:]:
    g = []
    for t in T:
        v = [ms[(a, "omp-clang", t)] / ms[(a, c, t)] for a in APPS if (a, c, t) in ms and (a, "omp-clang", t) in ms]
        g.append(f"{gm(v):.2f} (n={len(v)})" if v else "-")
    L.append(f"| {NAME[c]} | " + " | ".join(g) + " |")
L += ["", "## Best-thread comparison (each configuration at ITS best thread count, apps with all four cells)", "", "| config | geomean speedup vs omp-clang at its own best t | geomean of (best-t time / 128-thread time) | apps whose best t < 128 |", "|---|---:|---:|---:|"]
for c in CFGS:
    sp, rr, nb = [], [], 0
    for a in APPS:
        ok = {t: ms[(a, c, t)] for t in T if (a, c, t) in ms}; okb = {t: ms[(a, "omp-clang", t)] for t in T if (a, "omp-clang", t) in ms}
        if not ok or not okb: continue
        bt = min(ok, key=ok.get); sp.append(min(okb.values()) / ok[bt])
        if 128 in ok: rr.append(ok[bt] / ok[128]); nb += bt < 128
    L.append(f"| {NAME[c]} | {gm(sp):.2f} (n={len(sp)}) | {gm(rr):.2f} (n={len(rr)}) | {nb} |")
L += ["", "(The first column compares each configuration's best cell with omp-clang's best cell for the same app, i.e. best-vs-best.)"]
out = os.path.join(d, "sweep-summary.md"); open(out, "w").write("\n".join(L) + "\n")
with open(os.path.join(d, "sweep-results.csv"), "w", newline="") as f: csv.writer(f).writerows(rows)
print("\n".join(L))
