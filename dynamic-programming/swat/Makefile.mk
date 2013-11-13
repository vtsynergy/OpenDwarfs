#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += swat

swat_SOURCES = dynamic-programming/swat/alignments.cpp dynamic-programming/swat/prints.cpp \
				 dynamic-programming/swat/sequences.cpp dynamic-programming/swat/swat.cpp \
				 dynamic-programming/swat/param.cpp dynamic-programming/swat/timeRec.cpp






all_local += swat-all-local
exec_local += swat-exec-local

swat-all-local:
	cp $(top_srcdir)/dynamic-programming/swat/kernels.cl .

swat-exec-local:
	cp $(top_srcdir)/dynamic-programming/swat/kernels.cl ${DESTDIR}${bindir}
