#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += needle

needle_SOURCES = dynamic-programming/rodinia-nw/needle.c

all_local += rodinia-nw-all-local
exec_local += rodinia-nw-exec-local

rodinia-nw-all-local:
	cp $(top_srcdir)/dynamic-programming/rodinia-nw/needle_kernel.cl .

rodinia-nw-exec-local:
	cp $(top_srcdir)/dynamic-programming/rodinia-nw/needle_kernel.cl ${DESTDIR}${bindir}
