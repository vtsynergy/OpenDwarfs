#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += bfs

bfs_SOURCES = $(top_srcdir)/rodinia-bfs/bfs.cpp

all_local += bfs-all-local
exec_local += bfs-exec-local

bfs-all-local:
	cp $(top_srcdir)/rodinia-bfs/bfs_kernel.cl .

bfs-exec-local:
	cp $(top_srcdir)/rodinia-bfs/bfs_kernel.cl ${DESTDIR}${bindir}
