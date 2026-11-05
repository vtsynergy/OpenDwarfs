# =============================================================================
# Shared rules for every OpenCL benchmark Makefile (replaces autotools).
#
# A benchmark Makefile sets (before `include`):
#   NAME      binary name                        (required)
#   SRCS      sources relative to the bench dir  (required, .c and/or .cpp)
#   KERNELS   .cl/.h files the program fopen()s at run time; copied next to the
#             binary because every app loads them from its working directory
#   EXTRA_CPPFLAGS / EXTRA_LDLIBS   optional
#   (apps with extra programs, e.g. input generators, add their own rules)
# Then:   make            build into ./build
#         make run        run with RUN_ARGS on the selected OpenCL platform
#         make clean
#
# Toolchain knobs (all overridable on the command line / environment):
#   CC, CXX           host compilers (default gcc/g++; clang works)
#   OCL_PREFIX        installed OpenCL toolchain (toolchains/install-opencl-pocl.sh)
#   OCL_VENDOR        nvidia | pocl : which ICD to load at run time
#   OPENCL_INC/LIB    headers / libOpenCL.so (defaults: system ocl-icd)
# =============================================================================
OCL_ROOT    := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
BUILD_DIR   ?= build
CC          ?= gcc
CXX         ?= g++
OPT         ?= -O3 -U_FORTIFY_SOURCE -D_FORTIFY_SOURCE=0  # legacy Rodinia code overflows benign buffers; Ubuntu gcc enables fortify by default
OPENCL_INC  ?= /usr/include
OPENCL_LIB  ?= /usr/lib/x86_64-linux-gnu
OCL_PREFIX  ?= $(HOME)/.local/apps/opencl-$(shell uname -n | cut -d. -f1)
OCL_VENDOR  ?= nvidia
GPU_INDEX   ?= 0

CPPFLAGS_ALL := -I$(OCL_ROOT) -include CL/cl.h -I$(OCL_ROOT)/include -I$(OCL_ROOT)/opts -I$(OCL_ROOT)/../common -I$(OPENCL_INC) \
  -DCL_TARGET_OPENCL_VERSION=120 -DCL_USE_DEPRECATED_OPENCL_1_1_APIS \
  -DCL_USE_DEPRECATED_OPENCL_1_2_APIS $(EXTRA_CPPFLAGS)
WARN        ?= -w
CFLAGS_ALL  := $(OPT) $(WARN) -std=gnu89 -fcommon $(CFLAGS)
CXXFLAGS_ALL:= $(OPT) $(WARN) -fpermissive -std=gnu++14 $(EXTRA_CXXFLAGS) $(CXXFLAGS)
LDLIBS_ALL  := -L$(OPENCL_LIB) -Wl,-rpath,$(OPENCL_LIB) -lOpenCL -lm $(EXTRA_LDLIBS)

COMMON_SRCS := $(OCL_ROOT)/include/common_args.c $(OCL_ROOT)/include/rdtsc.c $(OCL_ROOT)/opts/opts.c
OBJ         := $(BUILD_DIR)/obj
define objname
$(OBJ)/$(subst /,_,$(subst ..,UP,$(basename $(1)))).o
endef
OBJS        := $(foreach s,$(SRCS) $(COMMON_SRCS),$(call objname,$(s)))
BIN         := $(BUILD_DIR)/$(NAME)

.PHONY: build run clean stage
build: $(BIN) stage

$(BIN): $(OBJS)
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(OBJS) -o $@ $(LDLIBS_ALL)

# One pattern rule per source language; the object-name mangling above is
# reversed by explicit rules generated here.
define mkrule
$(call objname,$(1)): $(1)
	@mkdir -p $(OBJ)
	$$(if $$(filter %.c,$(1)),$$(CC) $$(CFLAGS_ALL),$$(CXX) $$(CXXFLAGS_ALL)) $$(CPPFLAGS_ALL) -c $(1) -o $$@
endef
$(foreach s,$(SRCS) $(COMMON_SRCS),$(eval $(call mkrule,$(s))))

stage: $(BIN)
	@for k in $(KERNELS); do cp -f $$k $(BUILD_DIR)/; done; true

# ---- run ---------------------------------------------------------------------
# Selects the OpenCL implementation through ocl-icd (no recompilation needed):
#   OCL_VENDOR=nvidia -> NVIDIA driver (V100)   OCL_VENDOR=pocl -> PoCL (x86)
ifeq ($(OCL_VENDOR),pocl)
RUN_ENV := OCL_ICD_FILENAMES=$(OCL_PREFIX)/lib/libpocl.so LD_LIBRARY_PATH=$(OCL_PREFIX)/lib:$$LD_LIBRARY_PATH DWARF_TOOLCHAIN=pocl DWARF_DEVICE=cpu
else
RUN_ENV := OCL_ICD_FILENAMES=libnvidia-opencl.so.1 CUDA_VISIBLE_DEVICES=$(GPU_INDEX) DWARF_TOOLCHAIN=nvidia-ocl DWARF_DEVICE=gpu
endif
PLATFORM_ARGS ?= -p 0 -d 0 --
RUN_ARGS ?=
run: build
	cd $(BUILD_DIR) && env $(RUN_ENV) ./$(NAME) $(PLATFORM_ARGS) $(RUN_ARGS)

clean:
	rm -rf $(BUILD_DIR)
