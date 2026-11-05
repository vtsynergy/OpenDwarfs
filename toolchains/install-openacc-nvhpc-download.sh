#!/bin/bash
# =============================================================================
# FALLBACK (use install-openacc-nvhpc-module.sh when the site provides a module).
# Download + install the NVIDIA HPC SDK (nvc / nvc++ / nvfortran) -- the reference,
# fully supported OpenACC compiler. Targets x86 (-acc=multicore) and NVIDIA
# GPUs (-acc=gpu -gpu=ccXX). No root needed.
#
# Default install: $ACC_PREFIX = ~/.local/apps/openacc-<hostname>
#
# HOW TO MODIFY
#   other version : NVHPC_VER=26.3  (check the URL exists: see $URL below)
#   other CUDA    : NVHPC_CUDA_VER=12.9. Keep <13 for Volta (sm_70) -- CUDA 13
#                   removed Volta. For A100/H100 any version works.
#   other GPU     : nothing here; pass -gpu=ccXX when compiling (GPU_ARCH_SM).
#   other prefix  : ACC_PREFIX=/opt/nvhpc
#   smaller dl    : NVHPC_FLAVOR=cuda_12.9 (single-CUDA tarball, if published)
# The tarball is ~15 GB; it is kept in $DOWNLOAD_DIR and resumed if partial.
# =============================================================================
set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/env.sh"

NVHPC_FLAVOR="${NVHPC_FLAVOR:-cuda_multi}"
VER_NODOT="${NVHPC_VER//./}"
YEAR="20${NVHPC_VER%%.*}"
TARBALL="nvhpc_${YEAR}_${VER_NODOT}_Linux_x86_64_${NVHPC_FLAVOR}.tar.gz"
URL="${NVHPC_URL:-https://developer.download.nvidia.com/hpc-sdk/${NVHPC_VER}/${TARBALL}}"

if [ -x "$NVHPC_ROOT/compilers/bin/nvc" ]; then
  echo "NVHPC $NVHPC_VER already installed at $NVHPC_ROOT"; exit 0
fi
mkdir -p "$DOWNLOAD_DIR" "$ACC_PREFIX"
cd "$DOWNLOAD_DIR"
echo "Downloading $URL"
wget -c "$URL"
EXTRACT="$DOWNLOAD_DIR/${TARBALL%.tar.gz}"
[ -d "$EXTRACT" ] || tar xpzf "$TARBALL"
export NVHPC_SILENT=true
export NVHPC_INSTALL_DIR="$ACC_PREFIX"
export NVHPC_INSTALL_TYPE=single
"$EXTRACT/install"
"$NVHPC_ROOT/compilers/bin/nvc" --version | head -2
echo "Done. Delete $DOWNLOAD_DIR/$TARBALL and $EXTRACT to reclaim ~30 GB."
