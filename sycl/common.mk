# =============================================================================
# Shared rules for every SYCL benchmark Makefile (counterpart of opencl/common.mk).
#
# A benchmark Makefile sets (before `include`):
#   NAME       binary name                                              (required)
#   SRCS       SYCL translation units, compiled with icpx -fsycl        (required)
#   HOST_SRCS  host-only sources (.c: gcc, .cpp: g++) with the OpenCL tree's host flags
#   EXTRA_CPPFLAGS / EXTRA_CXXFLAGS / EXTRA_HOSTCXXFLAGS / EXTRA_CFLAGS / EXTRA_LDLIBS   optional
# Then:  make            build into ./build
#        make run SYCL_DEVICE=cpu|gpu RUN_ARGS="..."
#        make clean
# =============================================================================
SYCL_ROOT   := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
REPO_ROOT   := $(abspath $(SYCL_ROOT)/..)
BUILD_DIR   ?= build
CC          := gcc
OPT         ?= -O3
ONEAPI_ROOT ?= $(HOME)/.local/apps/oneapi-$(shell uname -n | cut -d. -f1)
ICPX        ?= $(ONEAPI_ROOT)/compiler/latest/bin/icpx
ONEAPI_LIBS := $(ONEAPI_ROOT)/compiler/latest/lib:$(ONEAPI_ROOT)/2025.2/lib

# NVIDIA target: built when the CUDA adapter is in the oneAPI install and a CUDA toolkit is found.
CUDA_PATH   ?= /opt/nvidia/hpc_sdk/Linux_x86_64/25.7/cuda/12.9
GPU_ARCH_SM ?= 90
SYCL_GPU    ?= $(if $(and $(wildcard $(ONEAPI_ROOT)/compiler/latest/lib/libur_adapter_cuda.so),$(wildcard $(CUDA_PATH)/nvvm)),1,0)
ifeq ($(SYCL_GPU),1)
SYCL_FLAGS  := -fsycl -fsycl-targets=spir64,nvptx64-nvidia-cuda --cuda-path=$(CUDA_PATH) \
               -Xsycl-target-backend=nvptx64-nvidia-cuda --cuda-gpu-arch=sm_$(GPU_ARCH_SM)
else
SYCL_FLAGS  := -fsycl
endif
# icpx defaults to fast-math; the OpenCL kernels are built with no relaxed-math option. These flags give the
# OpenCL default: no value-unsafe optimisation, FP contraction on, 2.5/3 ulp device divide/sqrt allowed.
FPFLAGS     := -fp-model=precise -fno-offload-fp32-prec-div -fno-offload-fp32-prec-sqrt

CPPFLAGS_ALL := -I$(SYCL_ROOT)/include -I$(REPO_ROOT)/common $(EXTRA_CPPFLAGS)
WARN        ?= -w
CFLAGS_ALL  := $(OPT) $(WARN) -std=gnu89 -fcommon $(EXTRA_CFLAGS)
HOSTCXX     := g++
HOSTCXXFLAGS = $(OPT) $(WARN) -fpermissive -std=gnu++14 $(EXTRA_HOSTCXXFLAGS)
CXXFLAGS_ALL =$(OPT) $(WARN) -std=c++17 $(FPFLAGS) $(SYCL_FLAGS) $(EXTRA_CXXFLAGS)
LDLIBS_ALL   = -lm -ldl $(EXTRA_LDLIBS)

OBJ  := $(BUILD_DIR)/obj
objname = $(OBJ)/$(subst /,_,$(basename $(1))).o
OBJS := $(foreach s,$(SRCS) $(HOST_SRCS),$(call objname,$(s)))
BIN  := $(BUILD_DIR)/$(NAME)
HDRS := $(SYCL_ROOT)/include/dwarf_sycl.hpp $(REPO_ROOT)/common/dwarf_timer.h $(REPO_ROOT)/common/dwarf_launch.h $(wildcard *.h *.hpp)

.PHONY: build run clean
build: $(BIN)

$(BIN): $(OBJS)
	env LD_LIBRARY_PATH=$(ONEAPI_LIBS):$$LD_LIBRARY_PATH $(ICPX) $(FPFLAGS) $(SYCL_FLAGS) $(OBJS) -o $@ $(LDLIBS_ALL)

define sycl_rule
$(call objname,$(1)): $(1) $(HDRS)
	@mkdir -p $(OBJ)
	env LD_LIBRARY_PATH=$(ONEAPI_LIBS):$$$$LD_LIBRARY_PATH $(ICPX) $$(CXXFLAGS_ALL) $$(CPPFLAGS_ALL) -c $(1) -o $$@
endef
define host_rule
$(call objname,$(1)): $(1) $(HDRS)
	@mkdir -p $(OBJ)
	$$(if $$(filter %.c,$(1)),$$(CC) $$(CFLAGS_ALL),$$(HOSTCXX) $$(HOSTCXXFLAGS)) $$(CPPFLAGS_ALL) -c $(1) -o $$@
endef
$(foreach s,$(SRCS),$(eval $(call sycl_rule,$(s))))
$(foreach s,$(HOST_SRCS),$(eval $(call host_rule,$(s))))

# ---- run ---------------------------------------------------------------------
# SYCL_DEVICE=cpu : Intel CPU OpenCL runtime of the oneAPI install (the same libintelocl.so as OCL_VENDOR=intel);
#                   a private ICD directory keeps other OpenCL CPU platforms (PoCL) invisible.
# SYCL_DEVICE=gpu : CUDA backend, physical GPU $(GPU_INDEX).
SYCL_DEVICE ?= cpu
GPU_INDEX   ?= 0
ICD_DIR     := $(abspath $(BUILD_DIR))/icd-intel
ifeq ($(SYCL_DEVICE),gpu)
RUN_ENV := ONEAPI_DEVICE_SELECTOR=cuda:gpu CUDA_VISIBLE_DEVICES=$(GPU_INDEX)
else
RUN_ENV := ONEAPI_DEVICE_SELECTOR=opencl:cpu OCL_ICD_VENDORS=$(ICD_DIR)
endif
RUN_ENV += LD_LIBRARY_PATH=$(ONEAPI_LIBS):$$LD_LIBRARY_PATH DWARF_SYCL_DEVICE=$(SYCL_DEVICE) DWARF_TOOLCHAIN=dpcpp
PLATFORM_ARGS ?= -p 0 -d 0 --
RUN_ARGS ?=
run: build
	@mkdir -p $(ICD_DIR) && echo $(ONEAPI_ROOT)/compiler/latest/lib/libintelocl.so > $(ICD_DIR)/intel.icd
	cd $(BUILD_DIR) && env $(RUN_ENV) ./$(NAME) $(PLATFORM_ARGS) $(RUN_ARGS)

clean:
	rm -rf build
