#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += gemnoui

gemnoui_SOURCES = n-body-methods/gem/dump_vertices.c n-body-methods/gem/floating_centers.c \
	n-body-methods/gem/radix_sort.c n-body-methods/gem/read_pqr.c \
	n-body-methods/gem/vector_math.c n-body-methods/gem/write_grid.c \
	n-body-methods/gem/extrapolate_bonds.c n-body-methods/gem/populate_stats.c \
	n-body-methods/gem/read_msms.c n-body-methods/gem/run_msms.c \
	n-body-methods/gem/write_avs.c n-body-methods/gem/write_xyzr.c \
	n-body-methods/gem/calculate_potential.cpp n-body-methods/gem/check_cmdline.cpp \
	n-body-methods/gem/estimate_a.cpp n-body-methods/gem/gem_no_ui.cpp \
	n-body-methods/gem/open_pqr_run_msms.cpp

gemnoui_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/n-body-methods/gem/include \
	-I$(top_srcdir)/n-body-methods/gem/include/visualize/dialogs -DNO_UI

all_local += gem-all-local
exec_local += gem-exec-local

gem-all-local:
	cp $(top_srcdir)/n-body-methods/gem/calculate_potential.cl .

gem-exec-local:
	cp $(top_srcdir)/n-body-methods/gem/calculate_potential.cl ${DESTDIR}${bindir}
