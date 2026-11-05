# OpenCL build for dense-linear-algebra/lud (replaces the autotools Makefile.mk).
# Usage: make [OCL_VENDOR=nvidia|pocl] [RUN_ARGS=...] {build|run|clean}   (see ../../common.mk)
NAME    := lud
SRCS    := lud.c common.c
KERNELS := lud_kernel.cl

include ../../common.mk
