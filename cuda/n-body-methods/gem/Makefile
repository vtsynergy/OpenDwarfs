# CUDA build for n-body-methods/gem.  Usage: make [CUDACC=...] [RUN_ARGS=...] {build|run}   (see ../../common.mk)
NAME    := gemnoui
SRCS    := dump_vertices.c floating_centers.c radix_sort.c read_pqr.c vector_math.c write_grid.c extrapolate_bonds.c populate_stats.c read_msms.c run_msms.c write_avs.c write_xyzr.c calculate_potential.cu check_cmdline.cpp estimate_a.cpp gem_no_ui.cu open_pqr_run_msms.cpp
EXTRA_CPPFLAGS := -I$(CURDIR)/include -I$(CURDIR)/include/visualize/dialogs -DNO_UI
# Pre-include libstdc++ headers: they must come before the min/max macros defined by the app headers
PRE_HDRS := exception iostream fstream sstream string vector algorithm cstdlib cmath
EXTRA_CXXFLAGS := $(addprefix -include ,$(PRE_HDRS))
EXTRA_CUDAFLAGS := $(addprefix -include ,$(PRE_HDRS))

include ../../common.mk
