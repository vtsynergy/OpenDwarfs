#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += bsort

bsort_SOURCES = $(top_srcdir)/graph-traversal/bitonic-sort/main.c

all_local += bitonic-sort-all-local
exec_local += bitonic-sort-exec-local

bitonic-sort-all-local:
	cp $(top_srcdir)/graph-traversal/bitonic-sort/bitonicSort.cl .

bitonic-sort-exec-local:
	cp $(top_srcdir)/graph-traversal/bitonic-sort/bitonicSort.cl ${DESTDIR}${bindir}
