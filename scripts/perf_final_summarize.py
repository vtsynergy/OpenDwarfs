#!/usr/bin/env python3
"""Final summary: results.csv (one row per app x device x config) and FINAL-REPORT.md (speedups pinned to omp-clang/LLVM,
overall geomeans, the pinned launch configurations and the environment). Usage: perf_final_summarize.py OUTDIR  (python3 -I)."""
import csv, sys, os, math, statistics as st, collections
out = sys.argv[1]
raw = os.path.join(out, "raw.csv")
APPS = "spmv bfs crc nw srad cfd kmeans nqueens hmm swat tdm gem fft lud".split()
CFGS = ["omp-clang", "omp-nvhpc", "opencl"]
EXTRA = "opencl-intel"
rows = list(csv.DictReader(open(raw)))
if any(r["config"] == EXTRA for r in rows): CFGS = CFGS + [EXTRA]
med = collections.defaultdict(list); stat = collections.defaultdict(collections.Counter)
for r in rows:
    k = (r["app"], r["device"], r["config"]); stat[k][r["status"]] += 1
    if r["status"] == "OK":
        med[k].append(float(r["median_ms"]))
def cell(a, d, c):
    v = med.get((a, d, c))
    if not v: return None
    m = st.median(v)
    return dict(ms=m, lo=min(v), hi=max(v), spread=(max(v) - min(v)) / m if m else 0.0, n=len(v))
def meaningful(x, base):  # pre-registered rule vs the baseline
    if not x or not base: return ""
    ratio = x["ms"] / base["ms"]
    if abs(ratio - 1) <= 0.10: return "no"
    if not (x["hi"] < base["lo"] or base["hi"] < x["lo"]): return "no"
    if x["spread"] >= 0.15 or base["spread"] >= 0.15: return "no (noisy)"
    return "yes"
with open(os.path.join(out, "results.csv"), "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["app", "device", "config", "status", "launches_ok", "median_ms", "launch_spread", "baseline_omp_clang_ms", "speedup_vs_omp_clang", "meaningful_vs_baseline"])
    for d in ("cpu", "gpu"):
        for a in APPS:
            base = cell(a, d, "omp-clang")
            for c in [k for k in CFGS if not (k == EXTRA and d == "gpu")]:
                x = cell(a, d, c)
                s = ",".join(f"{k}x{v}" for k, v in sorted(stat.get((a, d, c), {}).items())) or "NOT_RUN"
                if x:
                    w.writerow([a, d, c, s, x["n"], f"{x['ms']:.6f}", f"{x['spread']:.3f}", f"{base['ms']:.6f}" if base else "", f"{base['ms']/x['ms']:.4f}" if base else "", "baseline" if c == "omp-clang" else meaningful(x, base)])
                else:
                    w.writerow([a, d, c, s, 0, "", "", f"{base['ms']:.6f}" if base else "", "", ""])
def sp(a, d, c):
    b, x = cell(a, d, "omp-clang"), cell(a, d, c)
    return b["ms"] / x["ms"] if b and x else None
def fmt(v): return "-" if v is None else (f"{v:.2f}" if v >= 0.1 else f"{v:.3f}")
def gm(d, c):
    v = [sp(a, d, c) for a in APPS if sp(a, d, c)]
    return (math.exp(sum(map(math.log, v)) / len(v)), len(v)) if v else (None, 0)
L = ["# Final results: speedups pinned to omp-clang (LLVM) = 1.00", "",
     "Speedup = omp-clang median / config median (above 1 = faster than the LLVM OpenMP baseline). One input per app; median of the launch medians; cross-model rows mix implementation and toolchain.", ""]
env = os.path.join(out, "environment.txt")
if os.path.exists(env): L += ["## Environment", "```", open(env).read().strip(), "```", ""]
for d, o in (("cpu", "PoCL"), ("gpu", "OpenCL NVIDIA")):
    cf = [k for k in CFGS if not (k == EXTRA and d == "gpu")]
    nm = {"omp-clang": "omp-clang", "omp-nvhpc": "omp-nvhpc", "opencl": o, EXTRA: "Intel OpenCL"}
    L += [f"## {d.upper()}  (speedup vs omp-clang)", "",
          "| app | " + " | ".join(f"{nm[k]} ms" for k in cf) + " | " + " | ".join(f"{nm[k]} speedup" for k in cf[1:]) + " | status notes |",
          "|---|" + "---:|" * (len(cf) + len(cf) - 1) + "---|"]
    for a in APPS:
        c = [cell(a, d, k) for k in cf]
        notes = "; ".join(f"{k}: {','.join(sorted(stat[(a, d, k)]))}" for k, x in zip(cf, c) if x is None and stat.get((a, d, k))) or ("not run" if not any(stat.get((a, d, k)) for k in cf) else "")
        L.append(f"| {a} | " + " | ".join(fmt(x['ms']) if x else '-' for x in c) + " | " + " | ".join(fmt(sp(a, d, k)) for k in cf[1:]) + f" | {notes} |")
    gs = [gm(d, k) for k in cf[1:]]
    L += ["| **geomean (all apps with both cells)** | " + " | ".join([""] * len(cf)) + " | " + " | ".join(f"**{fmt(g)}** (n={n})" for g, n in gs) + " | |", ""]
lc = os.path.join(out, "LAUNCH-CONFIGS.md")
if os.path.exists(lc): L += ["## Pinned launch configurations", "", open(lc).read(), ""]
open(os.path.join(out, "FINAL-REPORT.md"), "w").write("\n".join(L))
print("wrote", os.path.join(out, "results.csv"), "and FINAL-REPORT.md")
