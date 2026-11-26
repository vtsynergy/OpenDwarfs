# SYCL build for n-body-methods/gem (see ../../common.mk)
NAME      := gemnoui
SRCS      := calculate_potential.cpp gem_no_ui.cpp
HOST_SRCS := dump_vertices.c floating_centers.c radix_sort.c read_pqr.c vector_math.c write_grid.c extrapolate_bonds.c populate_stats.c read_msms.c run_msms.c write_avs.c write_xyzr.c check_cmdline.cpp estimate_a.cpp open_pqr_run_msms.cpp
EXTRA_CPPFLAGS := -I$(CURDIR)/include -I$(CURDIR)/include/visualize/dialogs -DNO_UI
# Pre-include libstdc++ headers: they must come before the min/max macros defined by the app headers
PREINC := -include exception -include iostream -include fstream -include sstream -include string -include vector -include algorithm -include cstdlib -include cmath
EXTRA_HOSTCXXFLAGS := $(PREINC)
EXTRA_CXXFLAGS := -include sycl/sycl.hpp $(PREINC) -Wno-c++11-narrowing
include ../../common.mk
