#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += srad

srad_SOURCES = structured-grids/srad/srad.c

srad_LDFLAGS = -lm

all_local += srad-all-local
exec_local += srad-exec-local

srad-all-local:
	cp $(top_srcdir)/structured-grids/srad/srad_kernel.cl .
	cp $(top_srcdir)/structured-grids/srad/srad_kernel_opt_gpu.cl .

srad-exec-local:
	cp $(top_srcdir)/structured-grids/srad/srad_kernel.cl ${DESTDIR}${bindir}
	cp $(top_srcdir)/structured-grids/srad/srad_kernel_opt_gpu.cl ${DESTDIR}${bindir}
