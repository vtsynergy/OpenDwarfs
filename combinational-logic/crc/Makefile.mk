#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += crc 

crc_SOURCES = combinational-logic/crc/crc_algo.c

all_local += dwarf-crc-all-local
exec_local += dwarf-crc-exec-local

dwarf-crc-all-local:
	cp $(top_srcdir)/combinational-logic/crc/crc_algo_kernel.cl .

dwarf-crc-exec-local:
	cp $(top_srcdir)/combinational-logic/crc/crc_algo_kernel_2.cl ${DESTDIR}${bindir}
