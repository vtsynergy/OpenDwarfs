#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

#AUTOMAKE_OPTIONS = subdir-objects

ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -Werror

bin_PROGRAMS =
all_local =
exec_local =

if build_gem
include n-body-methods/gem/Makefile.mk
endif

if build_nw
include dynamic-programming/rodinia-nw/Makefile.mk
endif

if build_srad
include structured-grids/rodinia-srad/Makefile.mk
endif

if build_lud
include dense-linear-algebra/rodinia-lud/Makefile.mk
endif

if build_kmeans
include dense-linear-algebra/rodinia-kmeans/Makefile.mk
endif

if build_samplecl
include samplecl/Makefile.mk
endif

if build_crc
include combinational-logic/dwarf-crc/Makefile.mk
endif

if build_cfd
include unstructured-grids/dwarf-cfd/Makefile.mk
endif

if build_bfs
include graph-traversal/rodinia-bfs/Makefile.mk
endif

if build_bsort
include graph-traversal/bitonic-sort/Makefile.mk
endif

if build_oesort
include graph-traversal/oddeven-sort/Makefile.mk
endif

all-local: $(all_local)

install-exec-local: $(exec_local)

clean-local:
	rm -f *.cl
