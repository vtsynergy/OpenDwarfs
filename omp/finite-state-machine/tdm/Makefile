# OpenMP-offload translation of opencl/finite-state-machine/tdm (GPU temporal data mining)
BENCHMARK := tdm
ROOT_DIR := ../..
SRCS := tdm_ocl.cpp dataio.cpp
T := $(ROOT_DIR)/test/finite-state-machine/tdm
RUN_INPUT ?= $(T)/sim-64-size50.csv
RUN_INPUT2 ?= $(T)/ivl.txt
RUN_INPUT3 ?= $(T)/30-episodes.txt
APP_ARGS = $(RUN_INPUT) $(RUN_INPUT2) $(RUN_INPUT3) 128
include $(ROOT_DIR)/common/benchmark.mk
