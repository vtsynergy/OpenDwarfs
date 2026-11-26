# OpenDwarfs, SYCL tree

SYCL translations of the 14 applications of `opencl/` (the ground truth), built with Intel oneAPI DPC++ (`icpx -fsycl`).
Same directory layout, same program arguments, same validators, same timer (`common/dwarf_timer.h`, `model=sycl`).
Deliberate differences from the OpenCL tree are listed in [NOTES.md](NOTES.md).

## Toolchain

See `toolchains/install-sycl-oneapi.sh` (checks the install and documents how it is laid out; it installs nothing).
Defaults, overridable on the `make` command line or in the environment:

| Variable | Default | Meaning |
|---|---|---|
| `ONEAPI_ROOT` | `$HOME/.local/apps/oneapi-<host>` | oneAPI install (2025.2.1 tested) |
| `CUDA_PATH` | `/opt/nvidia/hpc_sdk/Linux_x86_64/25.7/cuda/12.9` | CUDA toolkit for the nvptx64 target |
| `GPU_ARCH_SM` | `90` | NVIDIA compute capability (H100) |
| `SYCL_GPU` | auto | `1` if the CUDA adapter and toolkit exist; `0` builds for the CPU only |

Only H100 (`sm_90`) and the CPU of one AMD EPYC 9454 host were tested.

## Build and run

```bash
make -C sycl/graph-traversal/bfs                               # one binary: spir64 (CPU) + sm_90 (GPU)
make -C sycl/graph-traversal/bfs run SYCL_DEVICE=cpu RUN_ARGS="$PWD/opencl/test/graph-traversal/bfs/graph4096.txt"
make -C sycl/graph-traversal/bfs run SYCL_DEVICE=gpu GPU_INDEX=0 RUN_ARGS=...

sycl/run_all.sh --device cpu --size small [--cpuset LIST] [--repeats N] [--only a,b]
sycl/tools/validate.sh --devices "cpu gpu" [--cpuset LIST] [--only a,b]
```

* The leading `-p P -d D -t T --` arguments are accepted as in OpenCL. The device type comes from `DWARF_SYCL_DEVICE`
  (set by `make run` from `SYCL_DEVICE`), otherwise from `-t` (1 = GPU, else CPU). `-d` indexes devices of that type.
  `-p` is ignored. The program prints the chosen device and exits with status 97 if that type is not available.
* `DWARF_WARMUP` / `DWARF_REPEATS` work as in the other trees.
* `crc` and `spmv` read `data/inputs/{crc,csr}_{small,medium}.dat` (override with `DATA=`).

## Validation

`tools/validate.sh` runs, per application and device: the application's own validator (`valid=1`), then the OpenCL
reference on the same device (CPU: Intel CPU OpenCL runtime; GPU: NVIDIA OpenCL) and a comparison. The reference is a
**copy** of the prebuilt `opencl/<app>/build/` run in the output directory (build `opencl/` first); nothing is written
under `opencl/`. `tools/ocl_dump.c` is preloaded into it and records every device-to-host transfer; the SYCL side
records the same stream when `DWARF_DUMP_DIR` is set (only used for validation, never for timing runs).

| App | Comparison with OpenCL | Why this criterion |
|---|---|---|
| bfs crc nw kmeans tdm | read-backs byte-identical | integer outputs (kmeans: memberships of every iteration, which determine the centres) |
| spmv | f32, rtol 1e-6, atol 1e-6 | one sum per row; only contraction could differ |
| srad | f32, rtol 1e-5, atol 1e-6 | GPU backends differ in the last ulps |
| lud | f32, rtol 1e-4, atol 1e-4 | rounding accumulates through the elimination |
| cfd | f32, rtol 1e-4, atol 1e-5 | atol covers momentum components near zero (their relative error is meaningless) |
| hmm | f32, rtol 1e-4, atol 1e-6 | all intermediate dot-product and model read-backs |
| gem | f32, rtol 1e-4, atol 1e-5 | `native_exp` has implementation-defined accuracy |
| fft | f32, rtol 1e-4, atol 1e-3 | `native_sin/cos`; error scales with the spectrum magnitude |
| nqueens | solution and unique counts printed | the raw result buffer contains unwritten (uninitialised) slots |
| swat | printed alignment block | output buffers are only defined up to the alignment length |

The tolerances are upper limits; the measured maximum differences are printed in `summary.csv`.
Results cover the small sizes of `opencl/run_all.sh`. The medium sizes have not been validated.
Additionally checked once: fft 1-D at all seven sizes (128 to 8192) and srad on non-square grids.
