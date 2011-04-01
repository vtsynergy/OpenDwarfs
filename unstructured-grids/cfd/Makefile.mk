#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += cfd 

cfd_SOURCES = unstructured-grids/cfd/cfd.cpp

all_local += dwarf-cfd-all-local
exec_local += dwarf-cfd-exec-local

dwarf-cfd-all-local:
	cp $(top_srcdir)/unstructured-grids/cfd/cfd_kernel.cl .

dwarf-cfd-exec-local:
	cp $(top_srcdir)/unstructured-grids/cfd/cfd_kernel.cl ${DESTDIR}${bindir}
