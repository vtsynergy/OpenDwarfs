# OpenCL build for finite-state-machine/tdm (replaces the autotools Makefile.mk).
# Usage: make [OCL_VENDOR=nvidia|pocl] [RUN_ARGS=...] {build|run|clean}   (see ../../common.mk)
NAME    := tdm
SRCS    := tdm_ocl.cpp dataio.cpp
KERNELS := tdm_ocl_kernel.cl types.h

include ../../common.mk
