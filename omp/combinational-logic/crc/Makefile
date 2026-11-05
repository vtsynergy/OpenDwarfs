# OpenMP-offload translation of opencl/combinational-logic/crc (Slice-by-8 Ethernet CRC32)
BENCHMARK := crc
ROOT_DIR := ../..
SRCS := src/crc_algo.c src-common/crc_formats.c
RUN_INPUT ?= $(ROOT_DIR)/../data/inputs/crc_small.dat
APP_ARGS = -i $(RUN_INPUT) -a
include $(ROOT_DIR)/common/benchmark.mk
