# OpenDwarfs — refactored layout (`refactor` branch)

One repository, three implementations of the Berkeley-dwarf benchmarks, each in its own
directory with the same conventions: a `Makefile` in every benchmark directory, a
`run_all.sh` for the whole model, and toolchains installed by scripts in `toolchains/`.

| Directory | Model | Origin | Toolchain (V100 + x86) |
|---|---|---|---|
| `opencl/`  | OpenCL 1.2 | `master` (autotools replaced by plain Makefiles) | NVIDIA driver ICD (GPU), PoCL (x86) |
| `omp/`     | OpenMP target offload | `omp-offload` | Clang/LLVM 21 built with NVPTX offload |
| `openacc/` | OpenACC | `openacc-beta` (2014, kept unmodified) | NVIDIA HPC SDK `nvc` (site module `nvhpc/25.5`) |
| `toolchains/` | — | new | download/install scripts + `env.sh` (all knobs) |
| `scripts/` | — | new | `run_all.sh <model|all>` dispatcher, comparison tools |

## Quick start (milan2: x86 EPYC + Tesla V100)
```bash
# 1. toolchains (once). Prefixes are ~/.local/apps/{llvm-offload,openacc,opencl}-<host>
toolchains/install-omp-llvm-offload.sh        # Clang 21 w/ OpenMP offload (~1 h build)
toolchains/install-openacc-nvhpc-module.sh    # records the site NVHPC module (no download)
toolchains/install-opencl-pocl.sh             # PoCL for x86 OpenCL (NVIDIA OpenCL = driver)
# 2. run
scripts/run_all.sh opencl --vendor nvidia --size small
scripts/run_all.sh opencl --vendor pocl   --size small
scripts/run_all.sh omp    --device v100   --size small
scripts/run_all.sh omp    --device cpu    --size small
scripts/run_all.sh openacc --run
# 3. single benchmark
make -C opencl/dense-linear-algebra/lud run OCL_VENDOR=nvidia RUN_ARGS="-i ../../test/dense-linear-algebra/lud/256.dat -v"
make -C omp/dense-linear-algebra/gemm run DEVICE=v100 RUN_SIZE=2048
```
All settings (target GPU arch, GPU index, versions, prefixes) are in `toolchains/env.sh`
and are documented, with instructions for retargeting, in `toolchains/README.md`.

## Notes
* Input data under `opencl/test` and `omp/test` are near-duplicates (the omp copy has three extra files).
* The `.aocx` FPGA binaries and the autotools build were dropped on this branch; FPGA `.cl`
  sources are kept in `opencl/`.
* Original OpenCL documentation: `opencl/README.legacy`; OpenMP documentation: `omp/README.md`.
