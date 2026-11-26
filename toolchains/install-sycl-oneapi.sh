#!/bin/bash
# =============================================================================
# SYCL toolchain: Intel oneAPI DPC++ (icpx -fsycl) as used by sycl/.
#
# This script DOCUMENTS the install and CHECKS it. It downloads nothing, installs nothing and accepts no licence.
#   toolchains/install-sycl-oneapi.sh            check the install at $INTEL_ONEAPI_ROOT and print what is there
#
# What the install on host "hudson" is (read from its logs and file dates):
#   * Intel oneAPI Base Toolkit 2025.2.1+41, installed 2025-11-22 19:52 (-0500) with the offline/online installer in
#     silent mode ("-s --eula accept --install-dir $HOME/intel/oneapi"; logs/installer.install.*.log).
#     Compiler: DPC++ 2025.2.1 (2025.2.0.20250806). It was later moved to $HOME/.local/apps/oneapi-hudson
#     (a symlink remains at $HOME/intel/oneapi). The layout has no top-level setvars.sh; 2025.2/oneapi-vars.sh exists.
#   * NVIDIA GPU support: three files added 2025-11-22 19:59, after the base install and not listed in its log:
#       compiler/2025.2/lib/libur_adapter_cuda.so.0.12.0   (+ two symlinks)
#       licensing/2025.2/LICENSE_oneAPI_for_nvidia_GPUs.md (Apache-2.0 with LLVM exceptions, "Unified Runtime Project")
#       compiler/2025.2/share/doc/compiler/oneAPI_for_nvidia_GPUs/index.html (redirects to developer.codeplay.com,
#                                                           "oneAPI for NVIDIA GPUs" 2025.2.0)
#     This matches Codeplay's separately distributed "oneAPI for NVIDIA GPUs" plugin. The files on disk do not record
#     which download was used or which terms were shown when it was obtained. Anyone reproducing the install must get
#     the plugin from developer.codeplay.com themselves and read the terms shown there.
#
# Reproducing on another host (manual steps, each a licensing decision for the person doing it):
#   1. Install the Intel oneAPI Base Toolkit (or the DPC++/C++ compiler) of the SAME version into
#      $APPS_DIR/oneapi-$SYSTEM_NAME. Do not move an install made elsewhere if you intend to use Intel's updater.
#   2. CPU target: nothing else is needed (the Intel CPU OpenCL runtime libintelocl.so ships with the compiler and
#      also works on AMD x86).
#   3. GPU target (optional): install the "oneAPI for NVIDIA GPUs" plugin matching the compiler version exactly, and
#      provide a CUDA toolkit for the device libraries (CUDA_PATH; here NVHPC 25.7's cuda/12.9).
#      sycl/common.mk enables the nvptx64 target only when both are present (override with SYCL_GPU=0|1).
#
# Modifying the install: treat it as read-only. The plugin must match the compiler version; replacing the compiler
# means reinstalling the plugin. The install is shared with the OpenCL tree's OCL_VENDOR=intel runs (libintelocl.so,
# TBB), so upgrading it changes that OpenCL CPU configuration too.
# =============================================================================
set -u
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$HERE/env.sh" >/dev/null
R="$INTEL_ONEAPI_ROOT"; L="$R/compiler/latest/lib"; rc=0
echo "oneAPI root: $R"
if [ -x "$R/compiler/latest/bin/icpx" ]; then
  LD_LIBRARY_PATH="$L:${INTEL_TBB_LIB:-}" "$R/compiler/latest/bin/icpx" --version | head -1
else echo "MISSING: $R/compiler/latest/bin/icpx"; rc=1; fi
if [ -e "$L/libintelocl.so" ]; then echo "CPU runtime: $L/libintelocl.so"; else echo "MISSING: Intel CPU OpenCL runtime"; rc=1; fi
if [ -e "$L/libur_adapter_cuda.so" ]; then echo "CUDA adapter: $(readlink -f "$L/libur_adapter_cuda.so")"
else echo "CUDA adapter not present: sycl/ builds for the CPU target only"; fi
if [ -x "$R/compiler/latest/bin/sycl-ls" ]; then
  echo "devices:"; env LD_LIBRARY_PATH="$L:${INTEL_TBB_LIB:-}" OCL_ICD_FILENAMES="$L/libintelocl.so" "$R/compiler/latest/bin/sycl-ls" 2>/dev/null
fi
exit $rc
