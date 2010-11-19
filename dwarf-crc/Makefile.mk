#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += crc 

crc_SOURCES = $(top_srcdir)/dwarf-crc/crc.c

all_local += dwarf-crc-all-local
exec_local += dwarf-crc-exec-local

dwarf-crc-all-local:
	cp $(top_srcdir)/dwarf-crc/crc_kernel.cl .

dwarf-crc-exec-local:
	cp $(top_srcdir)/dwarf-crc/crc_kernel.cl ${DESTDIR}${bindir}
