# OpenMP-offload translation of opencl/dynamic-programming/swat (affine-gap Smith-Waterman database search)
BENCHMARK := swat
ROOT_DIR := ../..
SRCS := swat.cpp alignments.cpp param.cpp prints.cpp sequences.cpp timeRec.cpp
RUN_INPUT ?= $(ROOT_DIR)/test/dynamic-programming/swat/query1K1
RUN_INPUT2 ?= $(ROOT_DIR)/test/dynamic-programming/swat/sampledb1K1
APP_ARGS = $(RUN_INPUT) $(RUN_INPUT2)
include $(ROOT_DIR)/common/benchmark.mk
