# OpenACC build (gemopenacc); see ../../common.mk for knobs (ACC_TARGET=gpu|multicore, GPU_CC=cc70).
NAME := gemopenacc
SRCS := calculate_potential.c check_cmdline.cpp cmain.c dump_vertices.c estimate_a.cpp extrapolate_bonds.c floating_centers.c gem_no_ui.cpp open_pqr_run_msms.cpp populate_stats.c radix_sort.c read_msms.c read_pqr.c run_msms.c vector_math.c write_avs.c write_grid.c write_xyzr.c
EXTRA_INCS := -Iinclude -Iinclude/visualize/dialogs -DNO_UI

include ../../common.mk
