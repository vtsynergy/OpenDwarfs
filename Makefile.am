#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

AUTOMAKE_OPTIONS = subdir-objects

ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -Werror

bin_PROGRAMS =
all_local =
exec_local =

if build_gem
include n-body-methods/gem/Makefile.mk
endif

if build_nw
include dynamic-programming/nw/Makefile.mk
endif

if build_swat
include dynamic-programming/swat/Makefile.mk
endif

if build_srad
include structured-grids/srad/Makefile.mk
endif

if build_lud
include dense-linear-algebra/lud/Makefile.mk
endif

if build_kmeans
include dense-linear-algebra/kmeans/Makefile.mk
endif

if build_samplecl
include samplecl/Makefile.mk
endif

if build_fft
include spectral-methods/fft/Makefile.mk
endif

if build_csr
include sparse-linear-algebra/SPMV/Makefile.mk
endif

if build_crc
include combinational-logic/crc/Makefile.mk
endif

if build_cfd
include unstructured-grids/cfd/Makefile.mk
endif

if build_bfs
include graph-traversal/bfs/Makefile.mk
endif

#if build_bsort
#include graph-traversal/bitonic-sort/Makefile.mk
#endif

#if build_oesort
#include graph-traversal/oddeven-sort/Makefile.mk
#endif

if build_tdm
include finite-state-machine/tdm/Makefile.mk
endif

if build_astar
include branch-and-bound/astar/Makefile.mk
endif

all-local: $(all_local)

install-exec-local: $(exec_local)

clean-local:
	rm -f *.cl
