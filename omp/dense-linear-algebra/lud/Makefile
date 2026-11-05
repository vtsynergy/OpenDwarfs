# OpenMP-offload translation of opencl/dense-linear-algebra/lud
BENCHMARK := lud
ROOT_DIR := ../..
SRCS := lud.c common.c
RUN_INPUT ?= $(ROOT_DIR)/test/dense-linear-algebra/lud/256.dat
APP_ARGS = -i $(RUN_INPUT)
include $(ROOT_DIR)/common/benchmark.mk
