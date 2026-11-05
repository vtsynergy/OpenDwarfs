# OpenMP-offload translation of opencl/n-body-methods/gem (GEM electrostatics, no-UI variant)
BENCHMARK := gem
ROOT_DIR := ../..
SRCS := dump_vertices.c floating_centers.c radix_sort.c read_pqr.c vector_math.c write_grid.c extrapolate_bonds.c \
        populate_stats.c read_msms.c run_msms.c write_avs.c write_xyzr.c \
        calculate_potential.cpp check_cmdline.cpp estimate_a.cpp gem_no_ui.cpp open_pqr_run_msms.cpp
HOST_ONLY_SRCS := $(filter-out calculate_potential.cpp,$(SRCS))
G := $(ROOT_DIR)/test/n-body-methods/gem
RUN_INPUT ?= $(G)/1uwo_A
APP_ARGS = $(RUN_INPUT) 80 1 0
EXTRA_CFLAGS := -std=gnu99 -Iinclude -Iinclude/visualize/dialogs -DNO_UI -w
EXTRA_CXXFLAGS := -w -Iinclude -Iinclude/visualize/dialogs -DNO_UI -include exception -include iostream -include fstream -include sstream -include string -include vector -include algorithm -include cstdlib -include cmath
EXTRA_LIBS := -lm
include $(ROOT_DIR)/common/benchmark.mk
