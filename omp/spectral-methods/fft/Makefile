# OpenMP-offload translation of opencl/spectral-methods/fft (SHOC FFT, 1-D path)
BENCHMARK := fft
ROOT_DIR := ../..
SRCS := fft.cpp
RUN_SIZE ?= 4
APP_ARGS = --pts $(RUN_SIZE) --passes 5
include $(ROOT_DIR)/common/benchmark.mk
