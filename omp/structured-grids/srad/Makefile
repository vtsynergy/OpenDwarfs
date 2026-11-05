# OpenMP-offload translation of opencl/structured-grids/srad (Rodinia SRAD)
BENCHMARK := srad
ROOT_DIR := ../..
SRCS := srad.c
RUN_SIZE ?= 512
APP_ARGS = $(RUN_SIZE) $(RUN_SIZE) 0 127 0 127 0.5 2
include $(ROOT_DIR)/common/benchmark.mk
