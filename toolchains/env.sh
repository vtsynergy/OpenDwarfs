#!/bin/bash
# =============================================================================
# Shared configuration for every toolchain/run script in this repository.
# Source it:   source toolchains/env.sh
# Every value is overridable from the calling environment, e.g.
#   SYSTEM_NAME=a100node GPU_ARCH_SM=80 source toolchains/env.sh
#
# TARGET HARDWARE (edit/override these to retarget the whole repo):
#   x86 host CPU            -> nothing to set (host arch is auto-detected)
#   NVIDIA V100 (Volta)     -> GPU_ARCH_SM=70  (cc70 / sm_70)
#   NVIDIA A100 / H100      -> GPU_ARCH_SM=80 / 90
# =============================================================================
export SYSTEM_NAME="${SYSTEM_NAME:-$(uname -n | cut -d. -f1)}"
export APPS_DIR="${APPS_DIR:-$HOME/.local/apps}"
export SOURCES_DIR="${SOURCES_DIR:-$HOME/SOURCES}"
export DOWNLOAD_DIR="${DOWNLOAD_DIR:-$SOURCES_DIR/toolchain-downloads}"

# --- GPU target ---------------------------------------------------------------
export GPU_ARCH_SM="${GPU_ARCH_SM:-70}"          # V100 = 70
export GPU_NAME="${GPU_NAME:-v100}"               # label used by omp/ DEVICE=
export GPU_INDEX="${GPU_INDEX:-0}"                # physical GPU used for runs

# --- install prefixes (one per programming model) ----------------------------
export OMP_PREFIX="${OMP_PREFIX:-$APPS_DIR/llvm-offload-$SYSTEM_NAME}"
export ACC_PREFIX="${ACC_PREFIX:-$APPS_DIR/openacc-$SYSTEM_NAME}"
export OCL_PREFIX="${OCL_PREFIX:-$APPS_DIR/opencl-$SYSTEM_NAME}"

# --- versions (bump here; see toolchains/README.md) ---------------------------
export LLVM_VER="${LLVM_VER:-21.1.8}"             # OpenMP offload (Clang)
export NVHPC_VER="${NVHPC_VER:-25.5}"            # OpenACC (nvc/nvc++; milan2 has it as module nvhpc/25.5)
export NVHPC_CUDA_VER="${NVHPC_CUDA_VER:-12.9}"   # must be <13: CUDA 13 dropped sm_70
export POCL_VER="${POCL_VER:-7.1}"                # OpenCL on x86 CPU

# Repo root (this file lives in <root>/toolchains)
export REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

# Where NVHPC lands once installed.
# Site module (milan2): `module load nvhpc/$NVHPC_VER` -> /opt/nvidia/hpc_sdk
export NVHPC_MODULE="${NVHPC_MODULE:-nvhpc/$NVHPC_VER}"
export NVHPC_SYSTEM_ROOT="${NVHPC_SYSTEM_ROOT:-/opt/nvidia/hpc_sdk/Linux_x86_64/$NVHPC_VER}"
# Download-installed copy lives under $ACC_PREFIX; the module copy is preferred.
if [ -d "$NVHPC_SYSTEM_ROOT" ]; then _d="$NVHPC_SYSTEM_ROOT"; else _d="$ACC_PREFIX/Linux_x86_64/$NVHPC_VER"; fi
export NVHPC_ROOT="${NVHPC_ROOT:-$_d}"
export NVHPC_CUDA_ROOT="${NVHPC_CUDA_ROOT:-$NVHPC_ROOT/cuda/$NVHPC_CUDA_VER}"
