#!/bin/bash
# =============================================================================
# OpenCL toolchain.
#
#  * NVIDIA GPU (V100): no compiler install is needed -- the NVIDIA driver ships
#    the OpenCL runtime + JIT (libnvidia-opencl.so, registered through
#    /etc/OpenCL/vendors/nvidia.icd). This script only *links* that ICD into the
#    private vendor directory below so the vendor can be selected per run.
#  * x86 CPU: PoCL (Portable Computing Language), the maintained open-source CPU
#    OpenCL implementation, is built from source and installed under $OCL_PREFIX.
#
# Result in $OCL_PREFIX:
#   lib/libpocl.so, include/ (PoCL), etc/OpenCL/vendors/{pocl,nvidia}.icd
# The benchmark Makefiles (opencl/common.mk) select the implementation with a private OCL_ICD_VENDORS directory holding exactly one
# .icd file (OCL_ICD_FILENAMES would only ADD to the system ICDs, so a CPU run could silently land on the GPU platform).
#
# Requirements (checked): cmake, ninja, an LLVM/Clang dev install, libhwloc-dev,
# ocl-icd-opencl-dev, opencl-headers, pkg-config.
#
# HOW TO MODIFY
#   PoCL version   : POCL_VER=7.1
#   LLVM used      : LLVM_CONFIG=/path/to/llvm-config  (PoCL supports a window of
#                    LLVM majors; PoCL 7.x -> LLVM 14..20. Default below is the
#                    existing ~/.local/apps/llvm-<host> if present.)
#   PoCL CUDA backend (PoCL driving the V100 through CUDA instead of the NVIDIA
#     OpenCL runtime -- useful for a compiler-vs-compiler comparison):
#                    ENABLE_POCL_CUDA=1 CUDA_ROOT=<toolkit with bin/ptxas>
#                    (needs an LLVM built with the NVPTX target).
#   other prefix   : OCL_PREFIX=/opt/opencl
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"

POCL_VER="${POCL_VER:-7.1}"
LLVM_CONFIG="${LLVM_CONFIG:-$APPS_DIR/llvm-$SYSTEM_NAME/bin/llvm-config}"
ENABLE_POCL_CUDA="${ENABLE_POCL_CUDA:-0}"
JOBS="${JOBS:-32}"

for t in cmake ninja pkg-config wget; do command -v $t >/dev/null || { echo "missing $t"; exit 1; }; done
[ -x "$LLVM_CONFIG" ] || { echo "ERROR: LLVM_CONFIG not found: $LLVM_CONFIG (set LLVM_CONFIG=...)"; exit 1; }
pkg-config --exists hwloc || { echo "ERROR: libhwloc-dev missing"; exit 1; }
LLVM_PREFIX="$("$LLVM_CONFIG" --prefix)"
echo "PoCL $POCL_VER  <- LLVM $("$LLVM_CONFIG" --version) at $LLVM_PREFIX"

mkdir -p "$DOWNLOAD_DIR" "$OCL_PREFIX"
cd "$DOWNLOAD_DIR"
SRC="pocl-$POCL_VER"
if [ ! -d "$SRC" ]; then
  wget -c -O "pocl-v$POCL_VER.tar.gz" "https://github.com/pocl/pocl/archive/refs/tags/v$POCL_VER.tar.gz"
  tar xzf "pocl-v$POCL_VER.tar.gz"
fi
BUILD="$DOWNLOAD_DIR/build-pocl-$POCL_VER-$SYSTEM_NAME"
EXTRA=()
if [ "$ENABLE_POCL_CUDA" = 1 ]; then EXTRA+=(-DENABLE_CUDA=ON -DCUDAToolkit_ROOT="${CUDA_ROOT:?set CUDA_ROOT}"); fi
cmake -S "$SRC" -B "$BUILD" -G Ninja \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="$OCL_PREFIX" \
  -DWITH_LLVM_CONFIG="$LLVM_CONFIG" \
  -DENABLE_ICD=ON -DENABLE_LOADABLE_DRIVERS=ON \
  -DKERNELLIB_HOST_CPU_VARIANTS=distro \
  -DENABLE_TESTS=OFF -DENABLE_EXAMPLES=OFF \
  -DCMAKE_INSTALL_RPATH="$OCL_PREFIX/lib;$LLVM_PREFIX/lib" "${EXTRA[@]}"
cmake --build "$BUILD" -j "$JOBS"
cmake --install "$BUILD"

# --- private ICD vendor directory ------------------------------------------
V="$OCL_PREFIX/etc/OpenCL/vendors"; mkdir -p "$V"
echo "$OCL_PREFIX/lib/libpocl.so" > "$V/pocl.icd"
if [ -r /etc/OpenCL/vendors/nvidia.icd ]; then
  cp /etc/OpenCL/vendors/nvidia.icd "$V/nvidia.icd"
else
  echo "WARNING: no system nvidia.icd found; NVIDIA OpenCL unavailable"
fi
cat > "$OCL_PREFIX/env.sh" <<EOT
export OCL_ICD_VENDORS="$V"
export LD_LIBRARY_PATH="$OCL_PREFIX/lib:$LLVM_PREFIX/lib:\${LD_LIBRARY_PATH:-}"
export OPENCL_INC=/usr/include
export OPENCL_LIB=/usr/lib/x86_64-linux-gnu
EOT
echo "Installed. Check with:  source $OCL_PREFIX/env.sh && clinfo -l"
