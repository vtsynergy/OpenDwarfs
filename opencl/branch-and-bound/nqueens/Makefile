# OpenCL build for branch-and-bound/nqueens (replaces the autotools Makefile.mk).
# Usage: make [OCL_VENDOR=nvidia|pocl] [RUN_ARGS=...] {build|run|clean}   (see ../../common.mk)
NAME    := nqueens
SRCS    := main.cpp nqueen_cl.cpp nqueen_cpu.cpp
KERNELS := kernels_nqueens.cl

include ../../common.mk
