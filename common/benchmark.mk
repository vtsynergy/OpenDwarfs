ROOT_DIR ?= ../..
COMMON_DIR := $(ROOT_DIR)/common
BUILD_DIR ?= build

ifeq ($(origin CXX),default)
CXX := clang++
endif
CXX ?= clang++
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

# Set CUDA_PATH if clang can't find libdevice for the NVPTX targets
# (e.g. "cannot find libdevice for sm_80"); needed on systems where the
# CUDA toolkit isn't clang's default search path (only the driver is
# installed, or CUDA ships inside something like the NVIDIA HPC SDK).
CUDA_PATH ?=
CUDA_PATH_FLAG := $(if $(CUDA_PATH),--cuda-path=$(CUDA_PATH),)

# Every build/run rule below takes the target device as an argument:
#   make build DEVICE=a100
#   make run   DEVICE=a100
# rather than a separate phony target per device.
DEVICE ?= cpu
KNOWN_DEVICES := cpu h100 h200 a100 mi210 mi300a mi100

ifeq ($(filter $(DEVICE),$(KNOWN_DEVICES)),)
$(error Unknown DEVICE '$(DEVICE)'. Known devices: $(KNOWN_DEVICES). \
  Add a new one by extending KNOWN_DEVICES and the arch/flags rules in common/benchmark.mk)
endif

ifeq ($(DEVICE),cpu)
  DEVICE_USE_TARGET := 0
  DEVICE_FLAGS :=
else ifneq ($(filter $(DEVICE),h100 h200 a100),)
  DEVICE_USE_TARGET := 1
  ifeq ($(DEVICE),a100)
    DEVICE_ARCH := $(A100_ARCH)
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

BIN := $(BUILD_DIR)/$(BENCHMARK)_$(DEVICE)

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

build: $(BIN)

$(BIN): main.cpp $(wildcard $(COMMON_DIR)/*.hpp) Makefile $(COMMON_DIR)/benchmark.mk
	mkdir -p $(BUILD_DIR)
	$(CXX) $(CXXFLAGS) -DDWARF_USE_TARGET=$(DEVICE_USE_TARGET) -DDWARF_TARGET_NAME='"$(DEVICE)"' \
	  $(DEVICE_FLAGS) -DDWARF_THREAD_LIMIT=$(DWARF_THREAD_LIMIT) \
	  $< -o $@ $(LDFLAGS)

run: build
	$(BIN) --size $(RUN_SIZE) --secondary-size $(RUN_SECONDARY_SIZE) --iterations $(RUN_ITERATIONS) --warmup $(RUN_WARMUP) --repeats $(RUN_REPEATS) $(RUN_INPUT_ARG) $(RUN_INPUT2_ARG) $(RUN_INPUT3_ARG)

all:
	for d in $(DEVICES); do $(MAKE) build DEVICE=$$d || exit $$?; done

clean:
	rm -rf $(BUILD_DIR)
