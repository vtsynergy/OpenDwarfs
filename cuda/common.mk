# =============================================================================
# Shared rules for every CUDA benchmark Makefile (counterpart of ../opencl/common.mk).
#
# A benchmark Makefile sets (before `include`):
#   NAME      binary name                                      (required)
#   SRCS      sources relative to the bench dir: .cu (CUDACC), .c (CC), .cpp (CXX)
#   EXTRA_CPPFLAGS / EXTRA_CXXFLAGS / EXTRA_CUDAFLAGS / EXTRA_LDLIBS   optional
# Then:   make            build into $(BUILD_DIR) (default ./build)
#         make run        run with RUN_ARGS on GPU_INDEX
#         make clean      prints what to delete (this file removes nothing itself)
#
# Compiler knobs:
#   CUDACC      CUDA compiler driver; any nvcc-compatible driver (default: NVIDIA nvcc from the HPC SDK).
#               SCALE:  make CUDACC=<scale-prefix>/targets/<target>/bin/nvcc BUILD_DIR=build-scale
#               (use a separate BUILD_DIR per compiler)
#   GPU_ARCH_SM target architecture, 90 = H100 (sm_90)
#   CC, CXX     host compilers for the plain C/C++ sources (default gcc/g++)
# =============================================================================
CUDA_TREE   := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
BUILD_DIR   ?= build
CC          ?= gcc
CXX         ?= g++
OPT         ?= -O3
NVHPC_CUDA_ROOT ?= /opt/nvidia/hpc_sdk/Linux_x86_64/25.7/cuda/12.9
CUDACC      ?= $(NVHPC_CUDA_ROOT)/bin/nvcc
GPU_ARCH_SM ?= 90
GPU_INDEX   ?= 0
CUDA_ARCHFLAGS ?= -gencode arch=compute_$(GPU_ARCH_SM),code=sm_$(GPU_ARCH_SM)

# Floating point: OpenCL (no build option) allows 2.5/3 ulp fp32 divide/sqrt; nvcc defaults to correctly rounded. OpenCL is the
# reference, so match it (with these flags the lud result is bitwise identical to NVIDIA OpenCL). `make FPFLAGS=` restores nvcc defaults.
FPFLAGS     ?= -prec-div=false -prec-sqrt=false

CPPFLAGS_ALL := -I$(CUDA_TREE)/include -I$(CUDA_TREE)/../common $(EXTRA_CPPFLAGS)
WARN         ?= -w
CFLAGS_ALL   := $(OPT) $(WARN) -std=gnu89 -fcommon $(CFLAGS)
CXXFLAGS_ALL := $(OPT) $(WARN) -fpermissive -std=gnu++14 $(EXTRA_CXXFLAGS) $(CXXFLAGS)
CUDAFLAGS_ALL := $(OPT) $(CUDA_ARCHFLAGS) -std=c++14 -w -Xcompiler -fpermissive $(FPFLAGS) $(EXTRA_CUDAFLAGS) $(CUDAFLAGS)
LDLIBS_ALL   := -lm -ldl $(EXTRA_LDLIBS)

OBJ  := $(BUILD_DIR)/obj
define objname
$(OBJ)/$(subst /,_,$(subst ..,UP,$(1))).o
endef
OBJS := $(foreach s,$(SRCS),$(call objname,$(s)))
BIN  := $(BUILD_DIR)/$(NAME)
HDRS := $(wildcard *.h *.cuh inc/*.h src/*.cuh $(CUDA_TREE)/include/*.h $(CUDA_TREE)/../common/*.h)

.PHONY: build run clean
build: $(BIN)

$(BIN): $(OBJS)
	@mkdir -p $(BUILD_DIR)
	$(CUDACC) $(CUDA_ARCHFLAGS) $(OBJS) -o $@ $(LDLIBS_ALL)

define mkrule
$(call objname,$(1)): $(1) $$(HDRS) $$(wildcard $$(dir $(1))*.cuh $$(dir $(1))*.h)
	@mkdir -p $(OBJ)
	$$(if $$(filter %.cu,$(1)),$$(CUDACC) $$(CUDAFLAGS_ALL),$$(if $$(filter %.c,$(1)),$$(CC) $$(CFLAGS_ALL),$$(CXX) $$(CXXFLAGS_ALL))) $$(CPPFLAGS_ALL) -c $(1) -o $$@
endef
$(foreach s,$(SRCS),$(eval $(call mkrule,$(s))))

RUN_ENV = CUDA_VISIBLE_DEVICES=$(GPU_INDEX) DWARF_TOOLCHAIN=$(DWARF_TOOLCHAIN_LABEL) DWARF_DEVICE=gpu
DWARF_TOOLCHAIN_LABEL ?= nvcc
PLATFORM_ARGS ?= -p 0 -d 0 --
RUN_ARGS ?=
run: build
	cd $(BUILD_DIR) && env $(RUN_ENV) ./$(NAME) $(PLATFORM_ARGS) $(RUN_ARGS)

clean:
	@echo "To clean, delete the directory $(abspath $(BUILD_DIR)) by hand."
