# OpenCL build for graphical-models/hmm (replaces the autotools Makefile.mk).
# Usage: make [OCL_VENDOR=nvidia|pocl] [RUN_ARGS=...] {build|run|clean}   (see ../../common.mk)
NAME    := bwa_hmm
SRCS    := main_bwa_hmm.c
KERNELS := bwa_hmm_opencl.cl

include ../../common.mk
