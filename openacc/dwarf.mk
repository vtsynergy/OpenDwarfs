# =============================================================================
# Shared rules for the 14 OpenACC translations of the OpenCL dwarfs (NVIDIA HPC SDK).
# (common.mk is the separate rule file of the legacy PGI/CoreTSAR sources.)
#   make build ACC_TARGET=gpu|multicore            (default gpu)
#   make run   ACC_TARGET=... [RUN_ARGS="..."]     same command line as the OpenCL binary
#                                                  (without the OpenCL "-p 0 -d 0 --" prefix)
# A benchmark Makefile sets NAME, SRCS, RUN_ARGS (default) and includes this file.
# =============================================================================
ACC_ROOT   := $(abspath $(dir $(lastword $(MAKEFILE_LIST))))
REPO_ROOT  := $(abspath $(ACC_ROOT)/..)
TEST_DIR   := $(REPO_ROOT)/opencl/test
DATA_DIR   := $(REPO_ROOT)/data/inputs
NVHPC_VER  ?= $(firstword $(foreach v,25.7 25.5 26.1 25.9,$(if $(wildcard /opt/nvidia/hpc_sdk/Linux_x86_64/$(v)),$(v))))
NVHPC_ROOT ?= /opt/nvidia/hpc_sdk/Linux_x86_64/$(NVHPC_VER)
# a CC/CXX inherited from the environment is ignored; override with `make ACC_CC=...`
ACC_CC     ?= $(NVHPC_ROOT)/compilers/bin/nvc
ACC_CXX    ?= $(NVHPC_ROOT)/compilers/bin/nvc++
ACC_TARGET ?= gpu
GPU_ARCH_SM ?= $(shell nvidia-smi --query-gpu=compute_cap --format=csv,noheader 2>/dev/null | head -1 | tr -d '. ')
GPU_CC     ?= cc$(if $(GPU_ARCH_SM),$(GPU_ARCH_SM),90)
OPT        ?= -O3
ifeq ($(ACC_TARGET),gpu)
ACCFLAGS   := -acc=gpu -gpu=$(GPU_CC)
USE_TARGET := 1
else ifeq ($(ACC_TARGET),multicore)
ACCFLAGS   := -acc=multicore
USE_TARGET := 0
else
$(error ACC_TARGET must be gpu or multicore)
endif
BUILD_DIR  ?= build
OBJDIR     := $(BUILD_DIR)/obj_$(ACC_TARGET)
BIN        := $(BUILD_DIR)/$(NAME)_$(ACC_TARGET)
CPPFLAGS_ALL := -I$(ACC_ROOT)/acc-common -I$(REPO_ROOT)/common -DDWARF_USE_TARGET=$(USE_TARGET) $(EXTRA_CPPFLAGS)
# -Mvect=noassoc: nvc -O3 otherwise vectorises floating-point reductions by reassociating them (host loops
# included); the OpenCL host (gcc -O3, no -ffast-math) does not. No -tp is given: nvc targets the build host CPU.
VECTFLAGS    ?= -Mvect=noassoc
CFLAGS_ALL   := $(OPT) $(VECTFLAGS) $(ACCFLAGS) -Minfo=accel $(EXTRA_CFLAGS)
CXXFLAGS_ALL := $(OPT) $(VECTFLAGS) -std=c++17 $(ACCFLAGS) -Minfo=accel $(EXTRA_CXXFLAGS)
LDFLAGS_ALL  := $(ACCFLAGS) $(EXTRA_LDFLAGS)

.PHONY: build run clean
build: $(BIN)
HDRS := $(wildcard *.h */*.h $(ACC_ROOT)/acc-common/*.h $(REPO_ROOT)/common/*.h)
$(OBJDIR)/%.c.o: %.c $(HDRS) Makefile $(ACC_ROOT)/dwarf.mk
	@mkdir -p $(dir $@)
	$(ACC_CC) $(CFLAGS_ALL) $(CPPFLAGS_ALL) -c $< -o $@
$(OBJDIR)/%.cpp.o: %.cpp $(HDRS) Makefile $(ACC_ROOT)/dwarf.mk
	@mkdir -p $(dir $@)
	$(ACC_CXX) $(CXXFLAGS_ALL) $(CPPFLAGS_ALL) -c $< -o $@
OBJS := $(addprefix $(OBJDIR)/,$(addsuffix .o,$(SRCS)))
$(BIN): $(OBJS)
	@mkdir -p $(BUILD_DIR)
	$(ACC_CXX) $(OBJS) -o $@ $(LDFLAGS_ALL) $(EXTRA_LIBS)

GPU_INDEX ?= 0
export DWARF_TOOLCHAIN := nvhpc-$(NVHPC_VER)
export DWARF_DEVICE := $(ACC_TARGET)
ifeq ($(ACC_TARGET),gpu)
RUN_ENV := CUDA_VISIBLE_DEVICES=$(GPU_INDEX)
else
RUN_ENV :=
endif
# the programs run inside the build directory (several write result files into the working directory)
run: build
	cd $(BUILD_DIR) && env $(RUN_ENV) ./$(NAME)_$(ACC_TARGET) $(RUN_ARGS)

# guarded: only the relative directory `build` or a subdirectory of it, without "..", is ever removed
clean:
	@case '$(BUILD_DIR)' in \
	  *..*|*[[:space:]]*) echo "clean refused: unsafe BUILD_DIR='$(BUILD_DIR)'" >&2; exit 1 ;; \
	  build|build/?*) rm -rf -- './$(BUILD_DIR)' ;; \
	  *) echo "clean refused: BUILD_DIR='$(BUILD_DIR)' is not build or build/<sub>" >&2; exit 1 ;; \
	esac
