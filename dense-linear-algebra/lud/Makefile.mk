#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += lud

lud_SOURCES = dense-linear-algebra/rodinia-lud/lud.c dense-linear-algebra/rodinia-lud/common.c

all_local += rodinia-lud-all-local
exec_local += rodinia-lud-exec-local

rodinia-lud-all-local:
	cp $(top_srcdir)/dense-linear-algebra/rodinia-lud/lud_kernel.cl .

rodinia-lud-exec-local:
	cp $(top_srcdir)/dense-linear-algebra/rodinia-lud/lud_kernel.cl ${DESTDIR}${bindir}
