# OpenACC tree

Two unrelated parts live here.

## 1. OpenACC translations of the 14 OpenCL dwarfs

`lud spmv bfs crc nw srad cfd kmeans nqueens hmm swat tdm gem fft`, translated from `../opencl/` (ground truth), with the
corrections of `../opencl/FIXES.md`. One source per application for both targets. Same command line as the OpenCL
programs (without the `-p 0 -d 0 --` prefix) and the shared timer `../common/dwarf_timer.h` (`model=openacc`).
Deliberate differences: `FIXES.md`. `gem` is in `n-body-methods/gem-dwarf` (`n-body-methods/gem` is legacy code).

Compiler: NVIDIA HPC SDK `nvc`/`nvc++` (found under `/opt/nvidia/hpc_sdk`; override `NVHPC_ROOT`, `ACC_CC`, `ACC_CXX`).
Rules: `dwarf.mk`; launch shapes: `acc-common/dwarf_acc.h` (from `../common/dwarf_launch.h`).

```bash
make build                          # all 14, GPU (-acc=gpu -gpu=cc<detected>)
make build ACC_TARGET=multicore     # all 14, host threads (-acc=multicore)
make -C structured-grids/srad run ACC_TARGET=gpu GPU_INDEX=0 RUN_ARGS="512 512 0 127 0 127 0.5 2"
./run_all.sh [--target gpu|multicore|both] [--size small|medium] [--only a,b]
make clean
```

Binaries are `build/<app>_<target>`. `run_all.sh` reports each application's own validator; comparisons against the OpenCL
outputs were done separately (see `FIXES.md`, "Validation").

## 2. Legacy 2014 PGI/CoreTSAR code (as imported from `openacc-beta`)

PolyBench-GPU ports, `helmholtz` (`structured-grid/`) and `n-body-methods/gem`. The sources are **unmodified**.
Build files: `common.mk`, per-benchmark `Makefile`s, and the `libctsar` rule that builds the CoreTSAR-lite runtime with
`nvc++` and adds stand-ins for the removed PGI OpenCL accessors. They have not been ported.

```bash
source ~/.local/apps/openacc-<host>/env.sh
make legacy-build                   # prints BUILD OK/FAIL per benchmark   (formerly `make build`)
make -C polybench-gpu/dense-linear-algebra/2MM run GPU_INDEX=0
make legacy-clean
```
Benchmarks that use OpenACC 1.0 syntax (`acc data region`, `acc for`) fail with a current `nvc` until ported; see
`build.log` in each directory.
