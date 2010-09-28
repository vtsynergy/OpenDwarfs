#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += sw dnagen

sw_SOURCES = $(top_srcdir)/swat/swat.cpp
dnagen_SOURCES = $(top_srcdir)/swat/dnagen.cpp

all_local += swat-all-local
exec_local += swat-exec-local

swat-all-local:
	cp $(top_srcdir)/swat/swat_kernels.cl .

swat-exec-local:
	cp $(top_srcdir)/swat/swat_kernels.cl ${DESTDIR}${bindir}
