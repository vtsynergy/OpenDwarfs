#!/usr/bin/env python3
"""Compare two read-back streams ({uint64 bytes, data}...) written by SYCL (dwarf_sycl.hpp) and OpenCL (ocl_dump.so).

usage: compare_readback.py ocl.bin sycl.bin [--type exact|f32] [--rtol R] [--atol A]
exact : byte-identical records.
f32   : records read as float32; element passes if |a-b| <= atol + rtol*max(|a|,|b|) (NaN/inf must match exactly).
Record count and record sizes must always match. Exit status 0 = match."""
import sys, struct, argparse, array, math
def records(path):
    d = open(path, 'rb').read(); o = 0; out = []
    while o < len(d):
        n = struct.unpack_from('<Q', d, o)[0]; o += 8
        out.append(d[o:o+n]); o += n
    return out
ap = argparse.ArgumentParser(); ap.add_argument('ref'); ap.add_argument('new')
ap.add_argument('--type', default='exact', choices=['exact', 'f32']); ap.add_argument('--rtol', type=float, default=0.0); ap.add_argument('--atol', type=float, default=0.0)
a = ap.parse_args()
R, N = records(a.ref), records(a.new)
if len(R) != len(N): print("DIFF record count: ocl=%d sycl=%d" % (len(R), len(N))); sys.exit(1)
if not R: print("DIFF no records"); sys.exit(1)
total = sum(map(len, R)); identical = 0; bad = 0; nel = 0; maxabs = 0.0; maxrel = 0.0
for k, (x, y) in enumerate(zip(R, N)):
    if len(x) != len(y): print("DIFF record %d size: ocl=%d sycl=%d" % (k, len(x), len(y))); sys.exit(1)
    if x == y: identical += 1; nel += len(x) // 4; continue
    if a.type == 'exact':
        bad += sum(1 for i in range(len(x)) if x[i] != y[i]); continue
    u = array.array('f'); u.frombytes(x[:len(x)//4*4]); v = array.array('f'); v.frombytes(y[:len(y)//4*4]); nel += len(u)
    for p, q in zip(u, v):
        if p == q: continue
        if not (math.isfinite(p) and math.isfinite(q)):
            if not (math.isnan(p) and math.isnan(q)): bad += 1
            continue
        d = abs(p - q); m = max(abs(p), abs(q)); maxabs = max(maxabs, d); maxrel = max(maxrel, d / m)
        if d > a.atol + a.rtol * m: bad += 1
unit = "bytes" if a.type == 'exact' else "floats"
print("%s records=%d bytes=%d identical_records=%d differing_%s_beyond_tol=%d max_abs=%.3g max_rel=%.3g (type=%s rtol=%g atol=%g)" %
      ("MATCH" if not bad else "DIFF", len(R), total, identical, unit, bad, maxabs, maxrel, a.type, a.rtol, a.atol))
sys.exit(1 if bad else 0)
