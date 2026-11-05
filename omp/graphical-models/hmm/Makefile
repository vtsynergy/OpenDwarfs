# OpenMP-offload translation of opencl/graphical-models/hmm (Baum-Welch HMM training)
BENCHMARK := hmm
ROOT_DIR := ../..
SRCS := main_bwa_hmm.c
RUN_SIZE ?= 60
APP_ARGS = -n $(RUN_SIZE) -v n
EXTRA_CFLAGS := -std=gnu99
EXTRA_LIBS := -lm
include $(ROOT_DIR)/common/benchmark.mk
