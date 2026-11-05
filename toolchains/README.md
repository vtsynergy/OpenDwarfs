# Toolchains

Scripts that download/build/register the compilers used by the three implementations, plus
`env.sh` (single source of configuration). Every script is idempotent-ish, prints what it
does, and is controlled only by environment variables, so retargeting never needs edits.

| Script | Installs | Where | Notes |
|---|---|---|---|
| `install-omp-llvm-offload.sh` | LLVM/Clang 21.1.8 with OpenMP target offload (host x86 + NVPTX) | `$OMP_PREFIX` (`~/.local/apps/llvm-offload-<host>`) | needs CUDA <13 for Volta; uses the CUDA bundled in NVHPC. ~1 h at `JOBS=32` |
| `install-openacc-nvhpc-module.sh` | nothing (verifies `module load nvhpc/25.5`, writes `$ACC_PREFIX/env.sh`) | `$ACC_PREFIX` (`~/.local/apps/openacc-<host>`) | preferred where a module exists |
| `install-openacc-nvhpc-download.sh` | NVIDIA HPC SDK from developer.download.nvidia.com (~15 GB) | `$ACC_PREFIX` | fallback for machines without the module |
| `install-opencl-pocl.sh` | PoCL (x86 CPU OpenCL); links NVIDIA's ICD | `$OCL_PREFIX` (`~/.local/apps/opencl-<host>`) | NVIDIA OpenCL itself is part of the driver |

## Retargeting (what to change)
* **Different NVIDIA GPU:** `GPU_ARCH_SM=80` (A100) / `90` (H100), `GPU_NAME=a100`. The omp build also knows
  `DEVICE=a100|h100|h200|v100` (see `omp/common/benchmark.mk`, variables `V100_ARCH`, `A100_ARCH`...).
  OpenACC uses `-gpu=cc$GPU_ARCH_SM`. **Volta needs CUDA < 13.**
* **Different GPU index:** `GPU_INDEX=3` (exported as `CUDA_VISIBLE_DEVICES` by all run scripts).
* **AMD GPU:** `ENABLE_AMDGPU=1 ROCM_ROOT=/opt/rocm` for the OpenMP script (the `omp/` Makefiles already have mi100/mi210/mi300a);
  OpenACC (nvc) and the NVIDIA ICD do not apply; PoCL can target AMD via its HIP backend (not set up here).
* **Different LLVM/PoCL/NVHPC versions:** `LLVM_VER`, `POCL_VER`, `NVHPC_VER`, `NVHPC_MODULE`.
* **Different install location / machine name:** `APPS_DIR`, `OMP_PREFIX`, `ACC_PREFIX`, `OCL_PREFIX`, `SYSTEM_NAME`.
* **Fewer/more cores for builds:** `JOBS`.
* **GCC as an alternative OpenACC/OpenMP offload compiler:** build GCC with
  `--enable-offload-targets=nvptx-none` (not scripted here; `nvc` is the supported OpenACC compiler).
