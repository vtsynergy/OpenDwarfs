# OpenMP-offload translation of opencl/unstructured-grids/cfd (Rodinia CFD Euler solver)
BENCHMARK := cfd
ROOT_DIR := ../..
SRCS := cfd.cpp
RUN_INPUT ?= $(ROOT_DIR)/test/unstructured-grids/cfd/fvcorr.domn.097K
APP_ARGS = $(RUN_INPUT)
include $(ROOT_DIR)/common/benchmark.mk
