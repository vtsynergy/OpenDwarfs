#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += astar

astar_SOURCES = Branch-and-Bound/astar/astar.c

all_local += astar-all-local
exec_local += astar-exec-local

astar-all-local:
	cp $(top_srcdir)/Branch-and-Bound/astar/astar.cl .

astar-exec-local:
	cp $(top_srcdir)/Branch-and-Bound/astar/astar.cl ${DESTDIR}${bindir}
