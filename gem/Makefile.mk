#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += gemnoui

gemnoui_SOURCES = $(top_srcdir)/gem/dump_vertices.c $(top_srcdir)/gem/floating_centers.c \
	$(top_srcdir)/gem/radix_sort.c $(top_srcdir)/gem/read_pqr.c \
	$(top_srcdir)/gem/vector_math.c $(top_srcdir)/gem/write_grid.c \
	$(top_srcdir)/gem/extrapolate_bonds.c $(top_srcdir)/gem/populate_stats.c \
	$(top_srcdir)/gem/read_msms.c $(top_srcdir)/gem/run_msms.c \
	$(top_srcdir)/gem/write_avs.c $(top_srcdir)/gem/write_xyzr.c \
	$(top_srcdir)/gem/calculate_potential.cpp $(top_srcdir)/gem/check_cmdline.cpp \
	$(top_srcdir)/gem/estimate_a.cpp $(top_srcdir)/gem/gem_no_ui.cpp \
	$(top_srcdir)/gem/open_pqr_run_msms.cpp

gemnoui_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/gem/include \
	-I$(top_srcdir)/gem/include/visualize/dialogs -DNO_UI

all_local += gem-all-local
exec_local += gem-exec-local

gem-all-local:
	cp $(top_srcdir)/gem/calculate_potential.cl .

gem-exec-local:
	cp $(top_srcdir)/gem/calculate_potential.cl ${DESTDIR}${bindir}
