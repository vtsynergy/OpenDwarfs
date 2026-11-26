ROOT_DIR ?= ../..
COMMON_DIR := $(ROOT_DIR)/common
BUILD_DIR ?= build

# Toolchain: the Clang installed by toolchains/install-omp-llvm-offload.sh
# (default prefix ~/.local/apps/llvm-offload-<host>); override with CXX=... or OMP_PREFIX=...
OMP_PREFIX ?= $(HOME)/.local/apps/llvm-offload-$(shell uname -n | cut -d. -f1)
TOOLCHAIN ?= clang
# The toolchain choice is authoritative: CXX/CC inherited from the environment (e.g. set by
# `module load nvhpc`) are ignored; only an explicit `make CXX=...` overrides the Clang path.
ifeq ($(TOOLCHAIN),clang)
ifneq ($(origin CXX),command line)
CXX := $(if $(wildcard $(OMP_PREFIX)/bin/clang++),$(OMP_PREFIX)/bin/clang++,clang++)
endif
endif
CXXFLAGS ?= -O3 -std=c++17 -fopenmp -I$(COMMON_DIR)
LDFLAGS ?= -fopenmp -Wl,--allow-shlib-undefined

# Embed an rpath to this CXX's own lib dirs so built binaries find their
# OpenMP-target-offload runtime (libomptarget.so etc.) regardless of the
# caller's LD_LIBRARY_PATH -- without this, a binary built against a
# non-default clang (the common case here) silently fails at run time with
# "error while loading shared libraries: libomptarget.so.* ... No such file".
CXX_RESOLVED := $(shell command -v $(CXX) 2>/dev/null)
CXX_ROOT := $(patsubst %/bin/$(notdir $(CXX)),%,$(CXX_RESOLVED))
ifneq ($(CXX_ROOT),)
LDFLAGS += -Wl,-rpath,$(CXX_ROOT)/lib -Wl,-rpath,$(CXX_ROOT)/lib/x86_64-unknown-linux-gnu

# Some LLVM installs (this offload toolchain included) ship libomptarget.so
# with its own DT_RUNPATH=$ORIGIN, which -- per ELF/ld.so semantics --
# shadows our executable's rpath when resolving *its* dependencies, so our
# rpath above alone doesn't reach libomptarget's own libLLVM.so dependency
# if that lib lives outside libomptarget.so's directory. Work around it by
# linking libLLVM directly into our binary too: as one of the executable's
# own direct NEEDED entries it resolves via our rpath, and the dynamic
# linker then reuses that already-loaded copy to satisfy libomptarget's
# request for the same soname instead of re-searching for it.
CXX_LIBLLVM := $(firstword $(wildcard $(CXX_ROOT)/lib/libLLVM.so.*) $(wildcard $(CXX_ROOT)/lib/libLLVM.so))
ifneq ($(CXX_LIBLLVM),)
LDFLAGS += -Wl,--no-as-needed -L$(CXX_ROOT)/lib -l:$(notdir $(CXX_LIBLLVM))
endif
endif

NVPTX_TARGET ?= nvptx64-nvidia-cuda
NVPTX_ARCH ?= sm_90
AMDGPU_TARGET ?= amdgcn-amd-amdhsa
AMDGPU_ARCH ?= gfx90a
MI300A_ARCH ?= gfx942
MI100_ARCH ?= gfx908
A100_ARCH ?= sm_80
V100_ARCH ?= sm_70

# Set CUDA_PATH if clang can't find libdevice for the NVPTX targets
# (e.g. "cannot find libdevice for sm_80"); needed on systems where the
# CUDA toolkit isn't clang's default search path (only the driver is
# installed, or CUDA ships inside something like the NVIDIA HPC SDK).
# Default: the CUDA 12.x bundled with the site NVHPC module (CUDA 13 dropped sm_70).
# Version: NVHPC_VER (set by toolchains/env.sh: the newest of 25.7 25.5 26.1 25.9 present on the host); else the highest installed.
NVHPC_ROOT_DIR ?= $(if $(NVHPC_VER),$(wildcard /opt/nvidia/hpc_sdk/Linux_x86_64/$(NVHPC_VER)),$(lastword $(sort $(wildcard /opt/nvidia/hpc_sdk/Linux_x86_64/2[0-9].[0-9]*))))
NVHPC_CUDA_DEFAULT := $(firstword $(wildcard $(NVHPC_ROOT_DIR)/cuda/12.9) $(wildcard $(NVHPC_ROOT_DIR)/cuda/12.*))
CUDA_PATH ?= $(NVHPC_CUDA_DEFAULT)
CUDA_PATH_FLAG := $(if $(CUDA_PATH),--cuda-path=$(CUDA_PATH),)

# Every build/run rule below takes the target device as an argument:
#   make build DEVICE=a100
#   make run   DEVICE=a100
# rather than a separate phony target per device.
DEVICE ?= cpu
KNOWN_DEVICES := cpu h100 h200 a100 v100 mi210 mi300a mi100

ifeq ($(filter $(DEVICE),$(KNOWN_DEVICES)),)
$(error Unknown DEVICE '$(DEVICE)'. Known devices: $(KNOWN_DEVICES). \
  Add a new one by extending KNOWN_DEVICES and the arch/flags rules in common/benchmark.mk)
endif

ifeq ($(DEVICE),cpu)
  DEVICE_USE_TARGET := 0
  DEVICE_FLAGS :=
else ifneq ($(filter $(DEVICE),h100 h200 a100 v100),)
  DEVICE_USE_TARGET := 1
  ifeq ($(DEVICE),a100)
    DEVICE_ARCH := $(A100_ARCH)
  else ifeq ($(DEVICE),v100)
    DEVICE_ARCH := $(V100_ARCH)
  else
    DEVICE_ARCH := $(NVPTX_ARCH)
  endif
  DEVICE_FLAGS := -fopenmp-targets=$(NVPTX_TARGET) -Xopenmp-target=$(NVPTX_TARGET) -march=$(DEVICE_ARCH) $(CUDA_PATH_FLAG)
else
  DEVICE_USE_TARGET := 1
  ifeq ($(DEVICE),mi100)
    DEVICE_ARCH := $(MI100_ARCH)
  else ifeq ($(DEVICE),mi300a)
    DEVICE_ARCH := $(MI300A_ARCH)
  else
    DEVICE_ARCH := $(AMDGPU_ARCH)
  endif
  DEVICE_FLAGS := -fopenmp-targets=$(AMDGPU_TARGET) -Xopenmp-target=$(AMDGPU_TARGET) -march=$(DEVICE_ARCH)
endif

# ---- compiler toolchain selection (toolchain-only comparisons) ---------------------
#   TOOLCHAIN=clang  LLVM/Clang 21 built by toolchains/install-omp-llvm-offload.sh (default)
#   TOOLCHAIN=nvhpc  NVIDIA HPC SDK nvc/nvc++ (OpenMP offload: -mp=gpu, host: -mp=multicore)
# The SAME source is compiled by both; apps that one toolchain cannot build are reported as
# blocked by the run scripts (no source changes per toolchain).
comma := ,
# NVHPC without `module load` (a loaded module would also leak its libomp into Clang runs).
NVHPC_BIN := $(NVHPC_ROOT_DIR)/compilers/bin
GPU_ARCH_SM_NUM ?= $(if $(GPU_ARCH_SM),$(GPU_ARCH_SM),70)
NVHPC_CUDA_SHORT ?= 12.9
REPO_COMMON := $(abspath $(ROOT_DIR)/../common)
ifeq ($(TOOLCHAIN),nvhpc)
  CXX := $(NVHPC_BIN)/nvc++
  CC  := $(NVHPC_BIN)/nvc
  CXXFLAGS := -O3 -std=c++17 -mp$(if $(filter-out cpu,$(DEVICE)),=gpu -gpu=cc$(GPU_ARCH_SM_NUM)$(if $(NVHPC_CUDA_SHORT),$(comma)cuda$(NVHPC_CUDA_SHORT),),=multicore) -I$(COMMON_DIR) -I$(REPO_COMMON)
  CFLAGS_APP := -O3 -mp$(if $(filter-out cpu,$(DEVICE)),=gpu -gpu=cc$(GPU_ARCH_SM_NUM)$(if $(NVHPC_CUDA_SHORT),$(comma)cuda$(NVHPC_CUDA_SHORT),),=multicore) -I$(COMMON_DIR) -I$(REPO_COMMON)
  LDFLAGS := -mp$(if $(filter-out cpu,$(DEVICE)),=gpu -gpu=cc$(GPU_ARCH_SM_NUM),=multicore) -Wl,-rpath,$(NVHPC_ROOT_DIR)/compilers/lib
  RUN_LIBPATH := $(NVHPC_ROOT_DIR)/compilers/lib:$(NVHPC_ROOT_DIR)/cuda/$(NVHPC_CUDA_SHORT)/lib64
  DEVICE_FLAGS :=
  TC_SUFFIX := _nvhpc
else
  CC := $(dir $(CXX))clang
  CFLAGS_APP := -O3 -fopenmp -I$(COMMON_DIR) -I$(REPO_COMMON)
  CXXFLAGS += -I$(REPO_COMMON)
  # [FIX-OMP-15] clang CPU builds targeted baseline x86-64 (no AVX/AVX-512/FMA) while nvc++ (-tp host) and the OpenCL CPU JITs target the host CPU
  ifeq ($(DEVICE),cpu)
    CXXFLAGS += -march=native
    CFLAGS_APP += -march=native
  endif
  TC_SUFFIX :=
  RUN_LIBPATH := $(CXX_ROOT)/lib:$(CXX_ROOT)/lib/x86_64-unknown-linux-gnu
endif
export DWARF_TOOLCHAIN := $(TOOLCHAIN)
export DWARF_DEVICE := $(DEVICE)

BIN := $(BUILD_DIR)/$(BENCHMARK)_$(DEVICE)$(TC_SUFFIX)

# DEVICES is a space-separated list used only by the `all` target, e.g.:
#   make all DEVICES="cpu a100 mi100"
DEVICES ?= cpu

RUN_SIZE ?= 1024
RUN_SECONDARY_SIZE ?= $(RUN_SIZE)
RUN_ITERATIONS ?= 1
RUN_REPEATS ?= 3
RUN_WARMUP ?= 1
DWARF_THREAD_LIMIT ?= 256
RUN_INPUT ?=
RUN_INPUT2 ?=
RUN_INPUT3 ?=
RUN_INPUT_ARG := $(if $(RUN_INPUT),--input $(RUN_INPUT),)
RUN_INPUT2_ARG := $(if $(RUN_INPUT2),--input2 $(RUN_INPUT2),)
RUN_INPUT3_ARG := $(if $(RUN_INPUT3),--input3 $(RUN_INPUT3),)

.PHONY: build run all clean

SRCS ?= main.cpp
# Sources without target regions may be compiled without the offload flags (needed when app headers define
# min/max macros that clash with clang's OpenMP-offload math wrappers). Link step still uses the offload flags.
HOST_ONLY_SRCS ?=
# Translated (OpenCL-derived) apps set APP_ARGS and take raw command-line arguments;
# the older proxy apps use the RUN_* variables below.
APP_ARGS ?=
OBJDIR := $(BUILD_DIR)/obj_$(DEVICE)$(TC_SUFFIX)
OBJS := $(addprefix $(OBJDIR)/,$(addsuffix .o,$(SRCS)))
DEFS := -DDWARF_USE_TARGET=$(DEVICE_USE_TARGET) -DDWARF_TARGET_NAME='"$(DEVICE)"' -DDWARF_THREAD_LIMIT=$(DWARF_THREAD_LIMIT)

build: $(BIN)

$(OBJDIR)/%.cpp.o: %.cpp $(wildcard $(COMMON_DIR)/*.hpp) $(wildcard $(REPO_COMMON)/*.h) Makefile $(COMMON_DIR)/benchmark.mk
	@mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) $(EXTRA_CXXFLAGS) $(DEFS) $(if $(filter $<,$(HOST_ONLY_SRCS)),,$(DEVICE_FLAGS)) -c $< -o $@
$(OBJDIR)/%.c.o: %.c $(wildcard $(REPO_COMMON)/*.h) Makefile $(COMMON_DIR)/benchmark.mk
	@mkdir -p $(dir $@)
	$(CC) $(CFLAGS_APP) $(DEFS) $(if $(filter $<,$(HOST_ONLY_SRCS)),,$(DEVICE_FLAGS)) $(EXTRA_CFLAGS) -c $< -o $@
$(BIN): $(OBJS)
	@mkdir -p $(BUILD_DIR)
	$(CXX) $(OBJS) -o $@ $(LDFLAGS) $(DEVICE_FLAGS) $(EXTRA_LIBS)

run: build
	LD_LIBRARY_PATH=$(RUN_LIBPATH) DWARF_WARMUP=$(RUN_WARMUP) DWARF_REPEATS=$(RUN_REPEATS) $(BIN) $(if $(APP_ARGS),$(APP_ARGS),--size $(RUN_SIZE) --secondary-size $(RUN_SECONDARY_SIZE) --iterations $(RUN_ITERATIONS) --warmup $(RUN_WARMUP) --repeats $(RUN_REPEATS) $(RUN_INPUT_ARG) $(RUN_INPUT2_ARG) $(RUN_INPUT3_ARG))

all:
	for d in $(DEVICES); do $(MAKE) build DEVICE=$$d || exit $$?; done

clean:
	rm -rf $(BUILD_DIR)
