# OpenACC translation of opencl/n-body-methods/gem (GEM electrostatics, no-UI variant).
# Named gem-dwarf: ../gem is the legacy 2014 PGI/CoreTSAR code and is left untouched.
NAME     := gem
SRCS     := dump_vertices.c floating_centers.c radix_sort.c read_pqr.c vector_math.c write_grid.c extrapolate_bonds.c \
            populate_stats.c read_msms.c run_msms.c write_avs.c write_xyzr.c \
            calculate_potential.cpp check_cmdline.cpp estimate_a.cpp gem_no_ui.cpp open_pqr_run_msms.cpp
RUN_ARGS ?= $(TEST_DIR)/n-body-methods/gem/1uwo_A 80 1 0
EXTRA_CPPFLAGS := -Iinclude -Iinclude/visualize/dialogs -DNO_UI -w
# libstdc++ headers must come before the min/max macros of the application headers
EXTRA_CXXFLAGS := --preinclude=exception --preinclude=iostream --preinclude=fstream --preinclude=sstream --preinclude=string \
                  --preinclude=vector --preinclude=algorithm --preinclude=cstdlib --preinclude=cmath
EXTRA_LIBS := -lm
include ../../dwarf.mk
