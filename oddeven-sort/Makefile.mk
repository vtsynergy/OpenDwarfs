#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += oesort

oesort_SOURCES = $(top_srcdir)/oddeven-sort/main2.c

all_local += oesort-all-local
exec_local += oesort-exec-local

oesort-all-local:
	cp $(top_srcdir)/oddeven-sort/OddEvenSort.cl .

oesort-exec-local:
	cp $(top_srcdir)/oddeven-sort/OddEvenSort.cl ${DESTDIR}${bindir}
