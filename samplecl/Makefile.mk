#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += scl

scl_SOURCES = $(top_srcdir)/samplecl/samplecl.c

all_local += samplecl-all-local
exec_local += samplecl-exec-local

samplecl-all-local:
	cp $(top_srcdir)/samplecl/samplecl_kernel.cl .

samplecl-exec-local:
	cp $(top_srcdir)/samplecl/samplecl_kernel.cl ${DESTDIR}${bindir}
