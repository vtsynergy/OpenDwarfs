#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += crc 
bin_PROGRAMS += createcrc

crc_LDFLAGS = -lm 
#@SEARCHFLAGS@ @LIBFLAGS@ @RPATHFLAGS@
crc_SOURCES = combinational-logic/crc/src/crc_algo.c combinational-logic/crc/src-common/crc_formats.c

createcrc_SOURCES = combinational-logic/crc/src-test/createcrc.c combinational-logic/crc/src-common/crc_formats.c

##createcrc does not need to be linked with any of the opencl common files
createcrc_LDADD = include/rdtsc.o include/common_args.o opts/opts.o
createcrc_LINK = $(CCLD) -lm -o $@

all_local += dwarf-crc-all-local
exec_local += dwarf-crc-exec-local

dwarf-crc-all-local:
	cp $(top_srcdir)/combinational-logic/crc/src/crc_kernel.cl .

dwarf-crc-exec-local:
	cp $(top_srcdir)/combinational-logic/crc/src/crc_algo_kernel_2.cl ${DESTDIR}${bindir}
