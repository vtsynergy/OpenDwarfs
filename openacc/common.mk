# =============================================================================
# Shared rules for the OpenACC benchmarks (NVIDIA HPC SDK `nvc`).
# Toolchain comes from the module recorded by toolchains/install-openacc-nvhpc-module.sh
# (source $ACC_PREFIX/env.sh first, or let run_all.sh do it).
#
# Knobs:  ACC_TARGET=gpu|multicore   GPU_CC=cc70   CC=nvc
#   gpu       : -acc=gpu -gpu=$(GPU_CC)           (V100: cc70)
#   multicore : -acc=multicore                     (x86 host threads)
# A benchmark Makefile sets NAME and SRCS (relative to its dir) and includes this file.
# CAVEAT: the sources are the 2014 PGI/CoreTSAR tree, kept unmodified; many use OpenACC
# 1.0-era syntax (`acc data region`) that current compilers reject. `make` reports it.
# =============================================================================
ACC_ROOT   := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
BUILD_DIR  ?= build
CC         ?= nvc
CXX        ?= nvc++
ACC_TARGET ?= gpu
GPU_CC     ?= cc70
OPT        ?= -O3
ifeq ($(ACC_TARGET),gpu)
ACCFLAGS   := -acc=gpu -gpu=$(GPU_CC)
else
ACCFLAGS   := -acc=multicore
endif
CTSAR_LIB  := $(ACC_ROOT)/build/libctsar.a
INCS       := -I$(ACC_ROOT)/include -I$(ACC_ROOT)/build/shim -I$(ACC_ROOT)/polybench-gpu/common $(EXTRA_INCS)
CFLAGS_ALL := $(OPT) $(ACCFLAGS) -mp -Dacc_device_radeon=acc_device_not_host $(INCS) $(CFLAGS)
CXXFLAGS_ALL := $(filter-out -Dacc_device_radeon%,$(CFLAGS_ALL)) -std=c++11 $(CXXFLAGS)
LDLIBS_ALL := $(CTSAR_LIB) -lstdc++ -lpthread -ldl -lm $(EXTRA_LIBS)
BIN        := $(BUILD_DIR)/$(NAME)

.PHONY: build run clean
build: $(BIN)
$(CTSAR_LIB):
	$(MAKE) -C $(ACC_ROOT) libctsar
OBJS := $(addprefix $(BUILD_DIR)/obj/,$(addsuffix .o,$(basename $(SRCS))))
$(BUILD_DIR)/obj/%.o: %.c
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS_ALL) -c $< -o $@
$(BUILD_DIR)/obj/%.o: %.cpp
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS_ALL) -c $< -o $@
$(BIN): $(OBJS) $(CTSAR_LIB)
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(ACCFLAGS) -mp $(OBJS) -o $@ $(LDLIBS_ALL)
RUN_ARGS ?=
GPU_INDEX ?= 0
run: build
	CUDA_VISIBLE_DEVICES=$(GPU_INDEX) $(BIN) $(RUN_ARGS)
clean:
	rm -rf $(BUILD_DIR)
