#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += bwa_hmm

bwa_hmm_LDFLAGS = -lm

bwa_hmm_SOURCES = graphical-models/hmm/main_bwa_hmm.c

all_local += bwa_hmm-all-local
exec_local += bwa_hmm-exec-local

bwa_hmm-all-local:
	cp $(top_srcdir)/graphical-models/hmm/bwa_hmm_opencl.cl .

bwa_hmm-exec-local:
	cp $(top_srcdir)/graphical-models/hmm/bwa_hmm_opencl.cl ${DESTDIR}${bindir}


	

