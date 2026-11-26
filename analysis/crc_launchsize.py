#!/usr/bin/env python3
"""crc launch-size test: vary the pages per kernel launch (crc -p N; 16384 pages in total, so launches per repetition = 16384/N) at ONE worker thread and measure
PoCL, Intel and Intel with the implicit vectorizer off. Prediction registered before the run (see output header). python3 -I analysis/crc_launchsize.py OUTDIR"""
import json, os, statistics as st, sys, time
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cpu_lowered as C
out = os.path.abspath(sys.argv[1]); os.makedirs(out, exist_ok=True)
P = [16, 64, 256, 1024, 4096, 16384]; N = 3
V = {"PoCL": ("pocl", {}), "Intel": ("intel", {}), "Intel (vectorizer off)": ("intel", {"CL_CONFIG_USE_VECTORIZER": "false"})}
res = {}
J = os.path.join(out, 'crc_launchsize.json')
if '--render-only' in sys.argv and os.path.exists(J):
    res = {(int(k.split('|')[0]), k.split('|')[1]): v for k, v in json.load(open(J)).items()}
else:
  for p in P:
    for name, (vend, env) in V.items():
        runs = []
        for i in range(N):
            e = dict(env)
            if vend == "pocl": e.update(POCL_CACHE_DIR=os.path.join(out, "pocl-cache"))
            runs.append(C.run_app("crc", vend, e, 3, 1, out, perf=True, extra_args=f"-p {p}"))
        res[(p, name)] = runs
  json.dump({f"{p}|{n}": v for (p, n), v in res.items()}, open(os.path.join(out, "crc_launchsize.json"), "w"), indent=1)
med = lambda p, n, k: st.median(r.get(k, 0) or 0 for r in res[(p, n)])
POCL, INTEL, OFF = "PoCL", "Intel", "Intel (vectorizer off)"
L = ["# crc launch-size test (1 worker thread; 16384 pages in total; launches per repetition = 16384 / pages-per-launch)", "",
     "Prediction registered before the run: if Intel's extra time on crc is per-launch runtime overhead, the Intel-minus-PoCL time difference shrinks roughly in proportion to the number of launches",
     "(about 7.7 us per launch at the default of 1024 launches) and is about zero or negative at <= 64 launches; if it is code-based, it stays roughly constant.", "",
     f"Median of {N} collections per cell (each: 1 warm-up + 3 timed repetitions; counters only over the timed repetitions). All differences are runtime minus PoCL (positive = slower than PoCL).", "",
     "| pages / launch | launches / rep | PoCL ms | Intel ms | Intel (vectorizer off) ms | Intel - PoCL ms | (Intel - PoCL) us per launch | Intel(vec off) - PoCL ms | Intel user cycles per launch | PoCL user cycles per launch |", "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
pts = []
for p in P:
    la = 16384 / p; tp, ti, to = (med(p, n, "median_ms") for n in (POCL, INTEL, OFF)); nt = 3 * la
    L.append(f"| {p} | {la:.0f} | {tp:.2f} | {ti:.2f} | {to:.2f} | {ti-tp:+.2f} | {(ti-tp)/la*1000:+.1f} | {to-tp:+.2f} | {med(p, INTEL, 'cycles')/nt:.3g} | {med(p, POCL, 'cycles')/nt:.3g} |"); pts.append((la, ti - tp, to - tp))
def fit(xy):
    n = len(xy); mx = sum(x for x, _ in xy) / n; my = sum(y for _, y in xy) / n
    sxx = sum((x - mx) ** 2 for x, _ in xy); b = sum((x - mx) * (y - my) for x, y in xy) / sxx; return my - b * mx, b
a0, b0 = fit([(x, y) for x, y, _ in pts]); a1, b1 = fit([(x, z) for x, _, z in pts])
a2, b2 = fit([(x, y) for x, y, _ in pts if x >= 4]); a3, b3 = fit([(x, z) for x, _, z in pts if x >= 4])
L += ["", "Least-squares fit of the difference = a + b x launches:",
      f"- Intel - PoCL, all 6 points: a = {a0:+.2f} ms, b = {b0*1000:+.2f} us per launch; excluding the single-launch point: a = {a2:+.2f} ms, b = {b2*1000:+.2f} us per launch.",
      f"- Intel (vectorizer off) - PoCL, all 6 points: a = {a1:+.2f} ms, b = {b1*1000:+.2f} us per launch; excluding the single-launch point: a = {a3:+.2f} ms, b = {b3*1000:+.2f} us per launch.",
      "(6 points, 2 parameters; descriptive, no significance claimed.)"]
open(os.path.join(out, "crc_launchsize.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
