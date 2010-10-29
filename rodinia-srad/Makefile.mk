#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += srad

srad_SOURCES = $(top_srcdir)/rodinia-srad/srad.c

all_local += rodinia-srad-all-local
exec_local += rodinia-srad-exec-local

rodinia-srad-all-local:
	cp $(top_srcdir)/rodinia-srad/srad_kernel.cl .

rodinia-srad-exec-local:
	cp $(top_srcdir)/rodinia-srad/srad_kernel.cl ${DESTDIR}${bindir}
