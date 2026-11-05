#!/bin/bash
# Generates the deterministic (-R) CRC and CSR inputs shared by the OpenCL and OpenMP trees into
# <repo>/data/inputs/ (gitignored), using the OpenCL generators createcrc / createcsr.
set -e
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
D="$ROOT/data/inputs"; mkdir -p "$D"
make -C "$ROOT/opencl/combinational-logic/crc" build >/dev/null
make -C "$ROOT/opencl/sparse-linear-algebra/SPMV" build >/dev/null
G1="$ROOT/opencl/combinational-logic/crc/build/createcrc"; G2="$ROOT/opencl/sparse-linear-algebra/SPMV/build/createcsr"
[ -f "$D/crc_small.dat" ]  || "$G1" -R -n 1024  -s 1024 -f "$D/crc_small.dat"  >/dev/null
[ -f "$D/crc_medium.dat" ] || "$G1" -R -n 16384 -s 1024 -f "$D/crc_medium.dat" >/dev/null
[ -f "$D/csr_small.dat" ]  || "$G2" -R -n 2048 -d 5000 -f "$D/csr_small.dat"  >/dev/null
[ -f "$D/csr_medium.dat" ] || "$G2" -R -n 8192 -d 5000 -f "$D/csr_medium.dat" >/dev/null
ls -la "$D"
