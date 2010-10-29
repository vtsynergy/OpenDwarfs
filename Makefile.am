#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -Werror

bin_PROGRAMS =
all_local =
exec_local =

#if build_gem
#include gem/Makefile.mk
#endif

#if build_matmult
#include matmult/Makefile.mk
#endif

#if build_mattran
#include mattran/Makefile.mk
#endif

#if build_swat
#include swat/Makefile.mk
#endif

if build_nw
include rodinia-nw/Makefile.mk
endif

if build_srad
include rodinia-srad/Makefile.mk
endif

if build_samplecl
include samplecl/Makefile.mk
endif

all-local: $(all_local)

install-exec-local: $(exec_local)

clean-local:
	rm -f *.cl
