#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += csr

csr_SOURCES = $(top_srcdir)/SPMV/csr.c

all_local += csr-all-local
exec_local += csr-exec-local

csr-all-local:
	cp $(top_srcdir)/SPMV/spmv_csr_kernel.cl .

csr-exec-local:
	cp $(top_srcdir)/SPMV/spmv_csr_kernel.cl ${DESTDIR}${bindir}
