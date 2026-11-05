#!/bin/bash
set -euo pipefail

# -------------------------------------------------
# Generic LLVM/Clang OpenMP-target-offload installer.
#
# Builds a single Clang toolchain capable of compiling
# OpenMP `#pragma omp target` code for the host CPU and
# for whichever GPU vendor(s) are available on this
# machine (NVIDIA via CUDA, AMD via ROCm, or both). GPU
# vendor support is auto-detected by default.
#
# This is the toolchain used to build and run the
# openmp-offload/ dwarfs in this repository; see its
# README.md for how to point `make`/`run.sh` at it via
# CXX=<install>/bin/clang++.
#
# If you use this project or this toolchain in your own
# work, please cite:
#
#   N. Chaudhury and W. Feng, "OpenDwarfs 2.0: Evolving
#   the Berkeley Dwarfs in OpenCL and OpenMP," in Proc.
#   IEEE High Performance Extreme Computing Conference
#   (HPEC), 2026. Best Student Paper Award.
#   [TBD: DOI / IEEE Xplore link / page numbers once the
#   proceedings are indexed]
#
# -------------------------------------------------
# Configuration
#
# Everything below is overridable via environment
# variables; the defaults auto-detect what this specific
# machine has available. Examples:
#
#   # Force an NVIDIA-only build, skip ROCm entirely:
#   ENABLE_AMDGPU=0 ./install-llvm-offload.sh
#
#   # Point at specific vendor toolkit locations:
#   CUDA_ROOT=/opt/cuda-12.9 ROCM_ROOT=/opt/rocm-6.4 \
#     ./install-llvm-offload.sh
#
#   # Use environment-modules to load vendor toolkits
#   # first (only attempted if these are set):
#   CUDA_MODULE=cuda/12.6.3 ROCM_MODULE=rocm/6.4.4 \
#     ./install-llvm-offload.sh
#
#   # Install somewhere specific instead of the
#   # hostname-scoped default:
#   INSTALL_DIR=/opt/llvm-offload ./install-llvm-offload.sh
# -------------------------------------------------
# --- repo integration: shared defaults come from toolchains/env.sh ----------
# CUDA_ROOT defaults to the CUDA bundled with the NVHPC module (needs bin/ptxas
# and nvvm/libdevice; Ubuntu's split /usr/lib/cuda layout does not work).
# Volta (sm_70) needs CUDA < 13.
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"
[ -d "${CUDA_ROOT:-}" ] || { [ -d "$NVHPC_CUDA_ROOT" ] && export CUDA_ROOT="$NVHPC_CUDA_ROOT"; }
export ENABLE_AMDGPU="${ENABLE_AMDGPU:-0}"      # x86 + NVIDIA only on this repo's targets
export DEFAULT_NVIDIA_ARCH="${DEFAULT_NVIDIA_ARCH:-sm_$GPU_ARCH_SM}"
export JOBS="${JOBS:-32}"                        # shared machine: do not take all cores
export INSTALL_DIR="${INSTALL_DIR:-$OMP_PREFIX}"

SOURCES_DIR="${SOURCES_DIR:-$HOME/SOURCES}"
LLVM_VER="${LLVM_VER:-21.1.8}"

LLVM_SRC_DIR="$SOURCES_DIR/llvm-$LLVM_VER"
LLVM_BUILD_DIR="$SOURCES_DIR/build-llvm-$LLVM_VER-$SYSTEM_NAME"
INSTALL_DIR="${INSTALL_DIR:-$HOME/.local/apps/llvm-offload-$SYSTEM_NAME}"

JOBS="${JOBS:-$(command -v nproc >/dev/null 2>&1 && nproc || echo 8)}"
LINK_JOBS="${LINK_JOBS:-8}"

# GPU architecture used only for the optional post-install
# smoke test and the printed usage examples -- an
# application's own build always selects its target
# architecture explicitly via --offload-arch, so these are
# not baked into the toolchain itself. Left blank (best
# effort auto-detect, otherwise unset) unless overridden.
DEFAULT_NVIDIA_ARCH="${DEFAULT_NVIDIA_ARCH:-}"
DEFAULT_AMDGPU_ARCH="${DEFAULT_AMDGPU_ARCH:-}"

# Set to 1/0 to force-enable/disable a vendor; default
# "auto" detects each via CUDA_ROOT/ROCM_ROOT (if set) or
# by finding nvcc/hipcc on PATH.
ENABLE_CUDA="${ENABLE_CUDA:-auto}"
ENABLE_AMDGPU="${ENABLE_AMDGPU:-auto}"

# Only attempted when explicitly set -- systems without
# environment-modules (or that already have the vendor
# toolkits on PATH/CUDA_ROOT/ROCM_ROOT) don't need these.
ROCM_MODULE="${ROCM_MODULE:-}"
CUDA_MODULE="${CUDA_MODULE:-}"

# Set to 1 to compile and run a multi-device OpenMP probe
# after installation, on whichever vendor(s) got enabled.
RUN_SMOKE_TEST="${RUN_SMOKE_TEST:-0}"

# -------------------------------------------------
# Host LLVM target backend, derived from `uname -m` so
# this isn't tied to x86_64-only clusters.
# -------------------------------------------------
case "$(uname -m)" in
  x86_64|amd64)   HOST_LLVM_TARGET="X86" ;;
  aarch64|arm64)  HOST_LLVM_TARGET="AArch64" ;;
  ppc64le|ppc64)  HOST_LLVM_TARGET="PowerPC" ;;
  riscv64)        HOST_LLVM_TARGET="RISCV" ;;
  *)
    echo "WARNING: unrecognized host architecture '$(uname -m)'," \
         "guessing LLVM target name 'X86'; override with HOST_LLVM_TARGET=..." >&2
    HOST_LLVM_TARGET="X86"
    ;;
esac
HOST_LLVM_TARGET="${HOST_LLVM_TARGET_OVERRIDE:-$HOST_LLVM_TARGET}"

# -------------------------------------------------
# Initialize the environment-modules command, only if
# module loading was actually requested.
# -------------------------------------------------
if [ -n "$ROCM_MODULE" ] || [ -n "$CUDA_MODULE" ]; then
  if ! type module >/dev/null 2>&1; then
    for init_script in \
      /etc/profile.d/modules.sh \
      /usr/share/Modules/init/bash \
      /usr/share/lmod/lmod/init/bash
    do
      if [ -f "$init_script" ]; then
        # shellcheck source=/dev/null
        source "$init_script"
        break
      fi
    done
  fi

  type module >/dev/null 2>&1 || {
    echo "ERROR: ROCM_MODULE/CUDA_MODULE were set, but the environment-modules" \
         "'module' command is unavailable"
    exit 1
  }

  if [ -n "$ROCM_MODULE" ]; then
    echo "Loading ROCm module: $ROCM_MODULE"
    module load "$ROCM_MODULE"
  fi

  if [ -n "$CUDA_MODULE" ]; then
    echo "Loading CUDA module: $CUDA_MODULE"
    module load "$CUDA_MODULE"
  fi

  echo
  echo "== Loaded modules =="
  module list 2>&1 || true
  echo
fi

# -------------------------------------------------
# Resolve ROCm installation and decide whether AMDGPU
# offload support gets built.
# -------------------------------------------------
ROCM_ROOT="${ROCM_ROOT:-${ROCM_PATH:-${HIP_PATH:-}}}"

if [ -z "$ROCM_ROOT" ] && command -v hipcc >/dev/null 2>&1; then
  HIPCC_PATH="$(readlink -f "$(command -v hipcc)")"
  ROCM_ROOT="$(dirname "$(dirname "$HIPCC_PATH")")"
fi
[ -n "$ROCM_ROOT" ] && ROCM_ROOT="$(readlink -f "$ROCM_ROOT")"

if [ "$ENABLE_AMDGPU" = "auto" ]; then
  if [ -n "$ROCM_ROOT" ] && [ -x "$ROCM_ROOT/bin/hipcc" ]; then
    ENABLE_AMDGPU=1
  else
    ENABLE_AMDGPU=0
  fi
fi

if [ "$ENABLE_AMDGPU" = "1" ]; then
  [ -n "$ROCM_ROOT" ] || {
    echo "ERROR: ENABLE_AMDGPU=1 but no ROCm installation was found."
    echo "       Set ROCM_ROOT (or ROCM_PATH/HIP_PATH), or ROCM_MODULE to load one."
    exit 1
  }

  HSA_INCLUDE="$ROCM_ROOT/include"
  HSA_LIB=""
  for candidate in \
    "$ROCM_ROOT/lib/libhsa-runtime64.so" \
    "$ROCM_ROOT/lib64/libhsa-runtime64.so"
  do
    if [ -f "$candidate" ]; then
      HSA_LIB="$candidate"
      break
    fi
  done

  [ -d "$ROCM_ROOT" ] || { echo "ERROR: ROCm root does not exist: $ROCM_ROOT"; exit 1; }
  [ -d "$HSA_INCLUDE" ] || { echo "ERROR: ROCm include directory not found: $HSA_INCLUDE"; exit 1; }
  [ -n "$HSA_LIB" ] || {
    echo "ERROR: libhsa-runtime64.so was not found under:"
    echo "       $ROCM_ROOT/lib"
    echo "       $ROCM_ROOT/lib64"
    exit 1
  }
  [ -x "$ROCM_ROOT/bin/hipcc" ] || { echo "ERROR: hipcc not found at $ROCM_ROOT/bin/hipcc"; exit 1; }
fi

# -------------------------------------------------
# Resolve CUDA installation and decide whether NVPTX
# offload support gets built.
# -------------------------------------------------
CUDA_ROOT="${CUDA_ROOT:-${CUDA_PATH:-${CUDA_HOME:-}}}"

if [ -z "$CUDA_ROOT" ] && command -v nvcc >/dev/null 2>&1; then
  NVCC_PATH="$(readlink -f "$(command -v nvcc)")"
  CUDA_ROOT="$(dirname "$(dirname "$NVCC_PATH")")"
fi
[ -n "$CUDA_ROOT" ] && CUDA_ROOT="$(readlink -f "$CUDA_ROOT")"

if [ "$ENABLE_CUDA" = "auto" ]; then
  if [ -n "$CUDA_ROOT" ] && [ -x "$CUDA_ROOT/bin/nvcc" ]; then
    ENABLE_CUDA=1
  else
    ENABLE_CUDA=0
  fi
fi

if [ "$ENABLE_CUDA" = "1" ]; then
  [ -n "$CUDA_ROOT" ] || {
    echo "ERROR: ENABLE_CUDA=1 but no CUDA installation was found."
    echo "       Set CUDA_ROOT (or CUDA_PATH/CUDA_HOME), or CUDA_MODULE to load one."
    exit 1
  }

  CUDA_LIBDEVICE="$(
    find "$CUDA_ROOT/nvvm/libdevice" \
      -maxdepth 1 -type f -name 'libdevice*.bc' \
      -print -quit 2>/dev/null
  )"

  [ -d "$CUDA_ROOT" ] || { echo "ERROR: CUDA root does not exist: $CUDA_ROOT"; exit 1; }
  [ -x "$CUDA_ROOT/bin/nvcc" ] || { echo "ERROR: nvcc not found at $CUDA_ROOT/bin/nvcc"; exit 1; }
  [ -n "$CUDA_LIBDEVICE" ] || {
    echo "ERROR: CUDA libdevice bitcode not found under:"
    echo "       $CUDA_ROOT/nvvm/libdevice"
    exit 1
  }
fi

if [ "$ENABLE_CUDA" != "1" ] && [ "$ENABLE_AMDGPU" != "1" ]; then
  echo "No usable CUDA or ROCm installation was found/enabled --" \
       "building a host-only (CPU target) OpenMP toolchain."
  echo "Set CUDA_ROOT/ROCM_ROOT (or *_MODULE, or ENABLE_CUDA=1/ENABLE_AMDGPU=1) to add GPU offload."
  echo
fi

# -------------------------------------------------
# Resolve host compilers.
#
# HOST_CC and HOST_CXX can be supplied explicitly, which
# is useful when this machine offers a specific GCC module
# or a non-default toolchain.
# -------------------------------------------------
HOST_CC="${HOST_CC:-$(command -v gcc || true)}"
HOST_CXX="${HOST_CXX:-$(command -v g++ || true)}"

[ -n "$HOST_CC" ] && [ -x "$HOST_CC" ] || { echo "ERROR: no usable host C compiler was found"; exit 1; }
[ -n "$HOST_CXX" ] && [ -x "$HOST_CXX" ] || { echo "ERROR: no usable host C++ compiler was found"; exit 1; }

# -------------------------------------------------
# Locate the host GCC installation for the Clang runtime
# sub-build.
# -------------------------------------------------
GCC_LIBGCC="$("$HOST_CXX" -print-libgcc-file-name)"
GCC_LIBGCC="$(readlink -f "$GCC_LIBGCC")"
GCC_INSTALL_DIR="${GCC_INSTALL_DIR:-$(dirname "$GCC_LIBGCC")}"

[ -f "$GCC_LIBGCC" ] || { echo "ERROR: GCC runtime library not found: $GCC_LIBGCC"; exit 1; }
[ -d "$GCC_INSTALL_DIR" ] || { echo "ERROR: GCC installation directory not found: $GCC_INSTALL_DIR"; exit 1; }

# -------------------------------------------------
# Configuration summary.
# -------------------------------------------------
echo "== LLVM OpenMP-offload build configuration =="
echo "SYSTEM_NAME     = $SYSTEM_NAME"
echo "HOST_LLVM_TARGET = $HOST_LLVM_TARGET"
echo "LLVM_VER        = $LLVM_VER"
echo "LLVM_SRC_DIR    = $LLVM_SRC_DIR"
echo "LLVM_BUILD_DIR  = $LLVM_BUILD_DIR"
echo "INSTALL_DIR     = $INSTALL_DIR"
echo
echo "ENABLE_AMDGPU   = $ENABLE_AMDGPU"
[ "$ENABLE_AMDGPU" = "1" ] && {
  echo "ROCM_ROOT           = $ROCM_ROOT"
  echo "HSA_INCLUDE         = $HSA_INCLUDE"
  echo "HSA_LIB             = ${HSA_LIB:-not found}"
  echo "DEFAULT_AMDGPU_ARCH = ${DEFAULT_AMDGPU_ARCH:-(unset -- pass --offload-arch explicitly)}"
}
echo
echo "ENABLE_CUDA     = $ENABLE_CUDA"
[ "$ENABLE_CUDA" = "1" ] && {
  echo "CUDA_ROOT           = $CUDA_ROOT"
  echo "CUDA_LIBDEVICE      = ${CUDA_LIBDEVICE:-not found}"
  echo "DEFAULT_NVIDIA_ARCH = ${DEFAULT_NVIDIA_ARCH:-(unset -- pass --offload-arch explicitly)}"
}
echo
echo "HOST_CC         = $HOST_CC"
echo "HOST_CXX        = $HOST_CXX"
echo "GCC_INSTALL_DIR = $GCC_INSTALL_DIR"
echo "JOBS            = $JOBS"
echo "LINK_JOBS       = $LINK_JOBS"
echo

# -------------------------------------------------
# Basic prerequisite checks.
# -------------------------------------------------
[ -d "$SOURCES_DIR" ] || { echo "ERROR: source directory does not exist: $SOURCES_DIR"; exit 1; }

for tool in cmake ninja tar; do
  command -v "$tool" >/dev/null 2>&1 || { echo "ERROR: missing required tool: $tool"; exit 1; }
done

# -------------------------------------------------
# Fetch LLVM sources when absent.
# -------------------------------------------------
if [ ! -d "$LLVM_SRC_DIR" ]; then
  command -v wget >/dev/null 2>&1 || { echo "ERROR: wget is required to fetch LLVM sources"; exit 1; }

  cd "$SOURCES_DIR"

  TARBALL="llvm-project-$LLVM_VER.src.tar.xz"
  URL="https://github.com/llvm/llvm-project/releases/download/llvmorg-$LLVM_VER/$TARBALL"

  echo "Fetching LLVM $LLVM_VER"
  echo "URL: $URL"

  rm -f "$TARBALL"
  wget "$URL"

  echo "Extracting $TARBALL"
  tar -xf "$TARBALL"

  EXTRACTED_DIR="llvm-project-$LLVM_VER.src"
  [ -d "$EXTRACTED_DIR" ] || {
    echo "ERROR: expected extracted directory was not created: $EXTRACTED_DIR"
    exit 1
  }

  mv "$EXTRACTED_DIR" "$LLVM_SRC_DIR"
  rm -f "$TARBALL"
fi

[ -d "$LLVM_SRC_DIR/llvm" ] || {
  echo "ERROR: LLVM source tree is incomplete: $LLVM_SRC_DIR/llvm does not exist"
  exit 1
}

# -------------------------------------------------
# Make enabled vendor toolchains discoverable.
#
# LLVM's CUDA and AMDGPU plugins can dynamically load
# libcuda and libhsa at execution time. The vendor
# installations are still needed for headers, device
# libraries, compilation, and runtime use.
# -------------------------------------------------
EXISTING_CMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH:-}"
EXISTING_LD_LIBRARY_PATH="${LD_LIBRARY_PATH:-}"

CMAKE_PREFIX_LIST=""
RUNTIME_LIBRARY_PATHS=()

if [ "$ENABLE_AMDGPU" = "1" ]; then
  CMAKE_PREFIX_LIST+="${CMAKE_PREFIX_LIST:+;}$ROCM_ROOT"
  RUNTIME_LIBRARY_PATHS+=("$ROCM_ROOT/lib" "$ROCM_ROOT/lib64")
  export ROCM_PATH="$ROCM_ROOT"
  export HIP_PATH="$ROCM_ROOT"
fi

if [ "$ENABLE_CUDA" = "1" ]; then
  CMAKE_PREFIX_LIST+="${CMAKE_PREFIX_LIST:+;}$CUDA_ROOT"
  RUNTIME_LIBRARY_PATHS+=("$CUDA_ROOT/lib64")
  export CUDA_HOME="$CUDA_ROOT"
  export CUDA_PATH="$CUDA_ROOT"
fi

if [ -n "$EXISTING_CMAKE_PREFIX_PATH" ]; then
  CMAKE_PREFIX_LIST+="${CMAKE_PREFIX_LIST:+;}${EXISTING_CMAKE_PREFIX_PATH//:/;}"
fi

NEW_LD_LIBRARY_PATH=""
for libdir in "${RUNTIME_LIBRARY_PATHS[@]:-}"; do
  if [ -n "$libdir" ] && [ -d "$libdir" ]; then
    NEW_LD_LIBRARY_PATH+="${NEW_LD_LIBRARY_PATH:+:}$libdir"
  fi
done
if [ -n "$EXISTING_LD_LIBRARY_PATH" ]; then
  NEW_LD_LIBRARY_PATH+="${NEW_LD_LIBRARY_PATH:+:}$EXISTING_LD_LIBRARY_PATH"
fi

export CMAKE_PREFIX_PATH="$CMAKE_PREFIX_LIST"
export LD_LIBRARY_PATH="$NEW_LD_LIBRARY_PATH"

# -------------------------------------------------
# Assemble the enabled-vendor plugin/target lists.
# -------------------------------------------------
LLVM_TARGETS="$HOST_LLVM_TARGET"
PLUGINS=""

if [ "$ENABLE_AMDGPU" = "1" ]; then
  LLVM_TARGETS+=";AMDGPU"
  PLUGINS+="${PLUGINS:+;}amdgpu"
fi
if [ "$ENABLE_CUDA" = "1" ]; then
  LLVM_TARGETS+=";NVPTX"
  PLUGINS+="${PLUGINS:+;}cuda"
fi

# -------------------------------------------------
# Configure LLVM.
#
#   LLVM_TARGETS_TO_BUILD:
#       Host backend plus whichever GPU backends are enabled.
#
#   LIBOMPTARGET_PLUGINS_TO_BUILD / _DLOPEN_PLUGINS:
#       Only the enabled vendor runtime plugin(s); dlopen'd
#       at runtime so the toolchain isn't linked against one
#       specific libcuda/libhsa version.
#
# GPU architecture selection happens when an application is
# compiled, via --offload-arch.
# -------------------------------------------------
mkdir -p "$LLVM_BUILD_DIR"

CMAKE_ARGS=(
  -S "$LLVM_SRC_DIR/llvm"
  -B "$LLVM_BUILD_DIR"
  -G Ninja

  -D CMAKE_BUILD_TYPE=Release
  -D CMAKE_INSTALL_PREFIX="$INSTALL_DIR"
  -D CMAKE_C_COMPILER="$HOST_CC"
  -D CMAKE_CXX_COMPILER="$HOST_CXX"

  -D LLVM_ENABLE_PROJECTS="clang;clang-tools-extra;lld"
  -D LLVM_ENABLE_RUNTIMES="openmp;offload"
  -D LLVM_TARGETS_TO_BUILD="$LLVM_TARGETS"
  -D CLANG_DEFAULT_CXX_STDLIB="libstdc++"
  -D "RUNTIMES_CMAKE_ARGS=-DCMAKE_C_FLAGS=--gcc-install-dir=${GCC_INSTALL_DIR};-DCMAKE_CXX_FLAGS=--gcc-install-dir=${GCC_INSTALL_DIR}"

  -D LLVM_BUILD_LLVM_DYLIB=ON
  -D LLVM_LINK_LLVM_DYLIB=ON
  -D BUILD_SHARED_LIBS=OFF
  -D LLVM_ENABLE_RTTI=ON
  -D LLVM_PARALLEL_LINK_JOBS="$LINK_JOBS"
  -D LLVM_INCLUDE_TESTS=OFF
  -D COMPILER_RT_INCLUDE_TESTS=OFF
  -D OFFLOAD_INCLUDE_TESTS=OFF
  -D LLVM_INCLUDE_EXAMPLES=OFF
  -D LLVM_INCLUDE_BENCHMARKS=OFF

  -D LIBOMP_OMPD_GDB_SUPPORT=OFF
  -D OPENMP_ENABLE_LIBOMPTARGET=ON
  -D LIBOMPTARGET_BUILD_SHARED_LIBS=ON
)

if [ -n "$CMAKE_PREFIX_LIST" ]; then
  CMAKE_ARGS+=(-D "CMAKE_PREFIX_PATH=$CMAKE_PREFIX_LIST")
fi
if [ -n "$PLUGINS" ]; then
  CMAKE_ARGS+=(
    -D "LIBOMPTARGET_PLUGINS_TO_BUILD=$PLUGINS"
    -D "LIBOMPTARGET_DLOPEN_PLUGINS=$PLUGINS"
  )
fi
[ "$ENABLE_CUDA" = "1" ] && CMAKE_ARGS+=(-D "CUDAToolkit_ROOT=$CUDA_ROOT")

cmake "${CMAKE_ARGS[@]}"

# -------------------------------------------------
# Build and install.
# -------------------------------------------------
cmake --build "$LLVM_BUILD_DIR" --parallel "$JOBS"
cmake --install "$LLVM_BUILD_DIR"

echo
echo "LLVM installed to: $INSTALL_DIR"
echo

# -------------------------------------------------
# Post-install sanity checks.
# -------------------------------------------------
echo "================================================="
echo "Post-install sanity checks"
echo "================================================="
echo

echo "-- Clang version --"
"$INSTALL_DIR/bin/clang" --version
echo

echo "-- Registered LLVM targets --"
"$INSTALL_DIR/bin/llc" --version | sed -n '/Registered Targets:/,$p'
echo

echo "-- Detected offload architectures --"
if [ -x "$INSTALL_DIR/bin/offload-arch" ]; then
  "$INSTALL_DIR/bin/offload-arch" || true
else
  echo "offload-arch was not installed"
fi
echo

echo "-- Host OpenMP/offload runtimes --"
find "$INSTALL_DIR" -maxdepth 8 -type f \
  \( -name 'libomp.so*' -o -name 'libomptarget.so*' -o -name 'liboffload.so*' \) \
  -print | sort
echo

if [ "$ENABLE_AMDGPU" = "1" ]; then
  echo "-- AMDGPU OpenMP plugin --"
  find "$INSTALL_DIR" -maxdepth 8 -type f \
    \( -name 'libomptarget.rtl.amdgpu.so*' -o -name 'libomptarget-amdgpu.so*' -o -name '*omptarget*amdgpu*.so*' \) \
    -print | sort
  echo
fi

if [ "$ENABLE_CUDA" = "1" ]; then
  echo "-- CUDA OpenMP plugin --"
  find "$INSTALL_DIR" -maxdepth 8 -type f \
    \( -name 'libomptarget.rtl.cuda.so*' -o -name 'libomptarget-cuda.so*' -o -name '*omptarget*cuda*.so*' \) \
    -print | sort
  echo
fi

echo "-- Device RTL bitcode --"
find "$INSTALL_DIR" -maxdepth 10 -type f \
  \( -name 'libomptarget-amdgpu*.bc' -o -name 'libomptarget-nvptx*.bc' -o -name 'libompdevice*.bc' -o -name 'libompdevice*.a' \) \
  -print | sort
echo

if command -v nvidia-smi >/dev/null 2>&1; then
  echo "-- NVIDIA visibility --"
  nvidia-smi --query-gpu=index,name,compute_cap --format=csv,noheader 2>/dev/null || nvidia-smi || true
  echo
fi

if command -v rocminfo >/dev/null 2>&1; then
  echo "-- AMD visibility --"
  rocminfo 2>/dev/null | grep -E '^[[:space:]]*Name:[[:space:]]+gfx' | sort -u || true
  echo
fi

# -------------------------------------------------
# Best-effort default arch detection for the smoke test
# and usage examples, if not already supplied.
# -------------------------------------------------
if [ "$ENABLE_CUDA" = "1" ] && [ -z "$DEFAULT_NVIDIA_ARCH" ] && command -v nvidia-smi >/dev/null 2>&1; then
  CC="$(nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -n1 | tr -d '.')"
  [ -n "$CC" ] && DEFAULT_NVIDIA_ARCH="sm_$CC"
fi
if [ "$ENABLE_AMDGPU" = "1" ] && [ -z "$DEFAULT_AMDGPU_ARCH" ] && command -v rocminfo >/dev/null 2>&1; then
  DEFAULT_AMDGPU_ARCH="$(rocminfo 2>/dev/null | grep -oE 'gfx[0-9a-fA-F]+' | head -n1 || true)"
fi

# -------------------------------------------------
# Optional runtime probe: enumerate every OpenMP target
# device exposed by libomptarget and launch a target
# region on each.
# -------------------------------------------------
if [ "$RUN_SMOKE_TEST" = "1" ]; then
  echo "================================================="
  echo "Building and running offload smoke test"
  echo "================================================="

  SMOKE_DIR="$LLVM_BUILD_DIR/offload-smoke-test"
  SMOKE_SRC="$SMOKE_DIR/offload_probe.c"
  SMOKE_EXE="$SMOKE_DIR/offload_probe"
  mkdir -p "$SMOKE_DIR"

  cat > "$SMOKE_SRC" <<'EOF'
#include <omp.h>
#include <stdio.h>

int main(void) {
  const int num_devices = omp_get_num_devices();
  printf("OpenMP target devices: %d\n", num_devices);
  if (num_devices <= 0) {
    fprintf(stderr, "ERROR: libomptarget found no accelerator devices\n");
    return 1;
  }
  for (int device = 0; device < num_devices; ++device) {
    int ran_on_initial_device = -1;
    int observed_device = -1;
#pragma omp target device(device) map(from : ran_on_initial_device, observed_device)
    {
      ran_on_initial_device = omp_is_initial_device();
      observed_device = omp_get_device_num();
    }
    printf("device %d: omp_get_device_num()=%d, initial_device=%d\n",
           device, observed_device, ran_on_initial_device);
    if (ran_on_initial_device != 0) {
      fprintf(stderr, "ERROR: target region for device %d fell back to the host\n", device);
      return 2;
    }
  }
  return 0;
}
EOF

  SMOKE_TARGETS=""
  SMOKE_FLAGS=()
  if [ "$ENABLE_AMDGPU" = "1" ]; then
    SMOKE_TARGETS+="${SMOKE_TARGETS:+,}amdgcn-amd-amdhsa"
    [ -n "$DEFAULT_AMDGPU_ARCH" ] && SMOKE_FLAGS+=(--offload-arch="$DEFAULT_AMDGPU_ARCH")
    SMOKE_FLAGS+=(--rocm-path="$ROCM_ROOT")
  fi
  if [ "$ENABLE_CUDA" = "1" ]; then
    SMOKE_TARGETS+="${SMOKE_TARGETS:+,}nvptx64-nvidia-cuda"
    [ -n "$DEFAULT_NVIDIA_ARCH" ] && SMOKE_FLAGS+=(--offload-arch="$DEFAULT_NVIDIA_ARCH")
    SMOKE_FLAGS+=(--cuda-path="$CUDA_ROOT")
  fi

  if [ -z "$SMOKE_TARGETS" ]; then
    echo "No GPU vendor enabled -- skipping smoke test (host-only toolchain has nothing to probe)."
  else
    "$INSTALL_DIR/bin/clang" -O2 -fopenmp \
      -fopenmp-targets="$SMOKE_TARGETS" \
      "${SMOKE_FLAGS[@]}" \
      "$SMOKE_SRC" -o "$SMOKE_EXE"

    echo "-- Embedded offload images --"
    "$INSTALL_DIR/bin/llvm-objdump" --offloading "$SMOKE_EXE"
    echo

    echo "-- Runtime execution --"
    OMP_TARGET_OFFLOAD=MANDATORY \
    LIBOMPTARGET_INFO="${LIBOMPTARGET_INFO:-4}" \
      "$SMOKE_EXE"
    echo
  fi
fi

# -------------------------------------------------
# Usage instructions.
# -------------------------------------------------
cat <<EOF

=================================================
LLVM OpenMP offload toolchain installed
=================================================

Installation:
  $INSTALL_DIR

Environment:
  export PATH="$INSTALL_DIR/bin:\$PATH"
  export LD_LIBRARY_PATH="$INSTALL_DIR/lib:${ROCM_ROOT:+$ROCM_ROOT/lib:$ROCM_ROOT/lib64:}${CUDA_ROOT:+$CUDA_ROOT/lib64:}\${LD_LIBRARY_PATH:-}"

CPU-only compilation (always available):
  clang++ app.cpp -O3 -fopenmp -o app-cpu
EOF

if [ "$ENABLE_AMDGPU" = "1" ]; then
cat <<EOF

AMD GPU compilation:
  clang++ app.cpp -O3 -fopenmp \\
    -fopenmp-targets=amdgcn-amd-amdhsa \\
    --offload-arch=${DEFAULT_AMDGPU_ARCH:-<your-gfxNNN>} \\
    --rocm-path="$ROCM_ROOT" \\
    -o app-amdgpu
EOF
fi

if [ "$ENABLE_CUDA" = "1" ]; then
cat <<EOF

NVIDIA GPU compilation:
  clang++ app.cpp -O3 -fopenmp \\
    -fopenmp-targets=nvptx64-nvidia-cuda \\
    --offload-arch=${DEFAULT_NVIDIA_ARCH:-<your-smNN>} \\
    --cuda-path="$CUDA_ROOT" \\
    -o app-nvidia
EOF
fi

if [ "$ENABLE_AMDGPU" = "1" ] && [ "$ENABLE_CUDA" = "1" ]; then
cat <<EOF

Dual-vendor fat binary:
  clang++ app.cpp -O3 -fopenmp \\
    -fopenmp-targets=amdgcn-amd-amdhsa,nvptx64-nvidia-cuda \\
    --offload-arch=${DEFAULT_AMDGPU_ARCH:-<your-gfxNNN>} \\
    --offload-arch=${DEFAULT_NVIDIA_ARCH:-<your-smNN>} \\
    --rocm-path="$ROCM_ROOT" \\
    --cuda-path="$CUDA_ROOT" \\
    -o app-dual
EOF
fi

cat <<EOF

Inspect embedded GPU images:
  llvm-objdump --offloading app-<...>

Runtime diagnostics:
  OMP_TARGET_OFFLOAD=MANDATORY LIBOMPTARGET_INFO=4 ./app-<...>

Build and run the optional device probe:
  RUN_SMOKE_TEST=1 ./install-llvm-offload.sh

To use this toolchain with the dwarfs in this repository:
  make build DEVICE=<device> CXX="$INSTALL_DIR/bin/clang++"
See README.md for the full list of known DEVICE names and
how to add a new GPU architecture.

EOF
