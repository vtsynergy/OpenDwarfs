#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += spmv

lud_SOURCES = $(top_srcdir)/spmv/spmv.c $(top_srcdir)/spmv/common.c

all_local += spmv-all-local
exec_local += spmv-lud-exec-local

rodinia-lud-all-local:
	cp $(top_srcdir)/spmv/spmv_kernel.cl .

rodinia-lud-exec-local:
	cp $(top_srcdir)/spmv/spmv_kernel.cl ${DESTDIR}${bindir}
