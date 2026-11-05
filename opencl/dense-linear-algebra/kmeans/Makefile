# OpenCL build for dense-linear-algebra/kmeans (replaces the autotools Makefile.mk).
# Usage: make [OCL_VENDOR=nvidia|pocl] [RUN_ARGS=...] {build|run|clean}   (see ../../common.mk)
NAME    := kmeans
SRCS    := kmeans.c cluster.c getopt.c kmeans_clustering.c kmeans_opencl.cpp rmse.c
KERNELS := kmeans_opencl_kernel.cl

include ../../common.mk
