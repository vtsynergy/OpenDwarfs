# OpenMP-offload translation of opencl/sparse-linear-algebra/SPMV (CSR sparse matrix-vector multiply)
BENCHMARK := spmv
ROOT_DIR := ../..
SRCS := src/csr.c src-common/sparse_formats.c src-common/ziggurat.c src-common/common.c
RUN_INPUT ?= $(ROOT_DIR)/../data/inputs/csr_small.dat
APP_ARGS = -i $(RUN_INPUT) -a
include $(ROOT_DIR)/common/benchmark.mk
