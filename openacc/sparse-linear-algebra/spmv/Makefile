# OpenACC translation of opencl/sparse-linear-algebra/SPMV (CSR sparse matrix-vector multiply)
NAME     := spmv
SRCS     := src/csr.c src-common/sparse_formats.c src-common/ziggurat.c src-common/common.c
RUN_ARGS ?= -i $(DATA_DIR)/csr_small.dat -a
EXTRA_LIBS := -lm
include ../../dwarf.mk
