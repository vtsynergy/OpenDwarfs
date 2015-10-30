#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += needle

needle_SOURCES = dynamic-programming/nw/needle.c

all_local += nw-all-local
exec_local += nw-exec-local

nw-all-local:
	cp $(top_srcdir)/dynamic-programming/nw/needle_kernel.cl .
	cp $(top_srcdir)/dynamic-programming/nw/needle_kernel_opt_gpu.cl .

nw-exec-local:
	cp $(top_srcdir)/dynamic-programming/nw/needle_kernel.cl ${DESTDIR}${bindir}
	cp $(top_srcdir)/dynamic-programming/nw/needle_kernel_opt_gpu.cl ${DESTDIR}${bindir}
