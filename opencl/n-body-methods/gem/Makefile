# OpenCL build for n-body-methods/gem (replaces the autotools Makefile.mk).
# Usage: make [OCL_VENDOR=nvidia|pocl] [RUN_ARGS=...] {build|run|clean}   (see ../../common.mk)
NAME    := gemnoui
SRCS    := dump_vertices.c floating_centers.c radix_sort.c read_pqr.c vector_math.c write_grid.c extrapolate_bonds.c populate_stats.c read_msms.c run_msms.c write_avs.c write_xyzr.c calculate_potential.cpp check_cmdline.cpp estimate_a.cpp gem_no_ui.cpp open_pqr_run_msms.cpp
KERNELS := calculate_potential.cl
EXTRA_CPPFLAGS := -I$(CURDIR)/include -I$(CURDIR)/include/visualize/dialogs -DNO_UI
# Pre-include libstdc++ headers: they must come before the min/max macros defined by the app headers
EXTRA_CXXFLAGS := -include exception -include iostream -include fstream -include sstream -include string -include vector -include algorithm -include cstdlib -include cmath

include ../../common.mk
