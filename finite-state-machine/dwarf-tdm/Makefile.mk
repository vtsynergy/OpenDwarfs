#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += tdm

tdm_SOURCES = finite-state-machine/dwarf-tdm/GpuTemporalDataMining.cpp

all_local += tdm-all-local
exec_local += tdm-exec-local

tdm-all-local:
	cp $(top_srcdir)/finite-state-machine/dwarf-tdm/GpuTemporalDataMining.cl .
	cp $(top_srcdir)/finite-state-machine/dwarf-tdm/strncpy.cl .
	cp $(top_srcdir)/finite-state-machine/dwarf-tdm/strncmp.cl .

tdm-exec-local:
	cp $(top_srcdir)/finite-state-machine/dwarf-tdm/GpuTemporalDataMining.cl ${DESTDIR}${bindir}
	cp $(top_srcdir)/finite-state-machine/dwarf-tdm/strncmp.cl ${DESTDIR}${bindir}
	cp $(top_srcdir)/finite-state-machine/dwarf-tdm/strncpy.cl ${DESTDIR}${bindir}
