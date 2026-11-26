#!/usr/bin/env python3
"""Launch-configuration audit: for every app, compare the (grid, block) actually launched by omp-clang and omp-nvhpc (from Nsight kernel traces).
   Both compile the same source with the same pinned clauses (common/dwarf_launch.h), so any difference is a toolchain interpretation of the clauses.
   python3 -I analysis/launch_audit.py BUNDLE_DIR"""
import glob, json, os, sys
d = sys.argv[1]; apps = sorted({os.path.basename(f).split(".")[0] for f in glob.glob(d + "/*.bundle.json")})
out = ["| app | kernel | launches/rep | clang grid x block | nvhpc grid x block | same? |", "|---|---|---:|---|---|---|"]
bad = []; n = 0
for a in apps:
    try: c = json.load(open(f"{d}/{a}.omp-clang.bundle.json")); v = json.load(open(f"{d}/{a}.omp-nvhpc.bundle.json"))
    except Exception: continue
    for k in sorted(set(c["kernels"]) | set(v["kernels"])):
        kc, kv = c["kernels"].get(k), v["kernels"].get(k)
        fmt = lambda x: "-" if not x else f"{x['grids'][0]} x {x['blocks'][0]}" + (f" (+{len(x['grids'])-1} more shapes)" if len(x['grids']) > 1 else "")
        same = bool(kc and kv and kc["grids"] == kv["grids"] and kc["blocks"] == kv["blocks"])
        n += 1
        if not same: bad.append((a, k))
        out.append(f"| {a} | {k[:34]} | {(kc or kv)['inst_per_rep']:.0f} | {fmt(kc)} | {fmt(kv)} | {'yes' if same else '**NO**'} |")
print("\n".join(out)); print(f"\n{n} kernels audited, {len(bad)} with different launch shapes between omp-clang and omp-nvhpc: " + ", ".join(f"{a}:{k[:20]}" for a, k in bad))
