#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

#This build currently assumes access to the PGI compiler suite for OpenACC
#support.  Patches for other compilers would be welcome, as would access to
#said compilers for our testing

# PGCC=${PGI_PATH}/bin/pgcc
# PGCC=${PGI_PATH}/bin/pgcc

lib_LIBRARIES += libctsar.a

#base
libctsar_a_SOURCES = openacc-tests/src/ctsar.cpp openacc-tests/src/pgi_shim.c openacc-tests/include/ctsar.h \
	openacc-tests/include/ctsar_impl.hpp openacc-tests/include/ctsar_dbg.h openacc-tests/include/pgi_shim.h \
	openacc-tests/include/device.hpp

ACC_INCLUDES=-I${top_srcdir}/openacc-tests/include
ACC_LIBS=-L. -lpthread -ldl

libctsar_a_CPPFLAGS=${ACC_INCLUDES}


bin_PROGRAMS += gemopenacc


gemopenacc_SOURCES = openacc-tests/n-body-methods/gem/dump_vertices.c openacc-tests/n-body-methods/gem/floating_centers.c \
	openacc-tests/n-body-methods/gem/radix_sort.c openacc-tests/n-body-methods/gem/read_pqr.c \
	openacc-tests/n-body-methods/gem/vector_math.c openacc-tests/n-body-methods/gem/write_grid.c \
	openacc-tests/n-body-methods/gem/extrapolate_bonds.c openacc-tests/n-body-methods/gem/populate_stats.c \
	openacc-tests/n-body-methods/gem/read_msms.c openacc-tests/n-body-methods/gem/run_msms.c \
	openacc-tests/n-body-methods/gem/write_avs.c openacc-tests/n-body-methods/gem/write_xyzr.c \
	openacc-tests/n-body-methods/gem/calculate_potential.c openacc-tests/n-body-methods/gem/check_cmdline.cpp \
	openacc-tests/n-body-methods/gem/estimate_a.cpp openacc-tests/n-body-methods/gem/gem_no_ui.cpp \
	openacc-tests/n-body-methods/gem/open_pqr_run_msms.cpp openacc-tests/n-body-methods/gem/cmain.c

gemopenacc_CPPFLAGS = $(AM_CPPFLAGS) -I$(top_srcdir)/openacc-tests/n-body-methods/gem/include \
	-I$(top_srcdir)/openacc-tests/n-body-methods/gem/include/visualize/dialogs -DNO_UI  -I${top_srcdir}/openacc-tests/include 

gemopenacc_LDFLAGS=${ACC_LIBS} 
gemopenacc_LDADD   = libctsar.a

bin_PROGRAMS += helmholtz
helmholtz_SOURCES = openacc-tests/structured-grid/helmholtz.c
helmholtz_LDFLAGS=${ACC_LIBS} 
helmholtz_LDADD   = libctsar.a
helmholtz_CPPFLAGS = -I${top_srcdir}/openacc-tests/include
nodist_EXTRA_helmholtz_SOURCES   = dummy.cpp


#polybench derivatives
bin_PROGRAMS += 2dconv 2mm 3dconv 3mm atax bicg corr covar fdtd2d gemm gesummv gramschm mvt syrk syr2k

#structured grid
2dconv_SOURCES = openacc-tests/polybench-gpu/structured-grid/2DCONV/twodconv.c
3dconv_SOURCES = openacc-tests/polybench-gpu/structured-grid/3DCONV/threedconv.c
fdtd2d_SOURCES = openacc-tests/polybench-gpu/structured-grid/FDTD-2D/fdtd2d.c

#dense linear algebra
2mm_SOURCES      = openacc-tests/polybench-gpu/dense-linear-algebra/2MM/twomm.c
3mm_SOURCES      = openacc-tests/polybench-gpu/dense-linear-algebra/3MM/threemm.c
atax_SOURCES     = openacc-tests/polybench-gpu/dense-linear-algebra/ATAX/atax.c
bicg_SOURCES     = openacc-tests/polybench-gpu/dense-linear-algebra/BICG/bicg.c
gemm_SOURCES     = openacc-tests/polybench-gpu/dense-linear-algebra/GEMM/gemm.c
gesummv_SOURCES  = openacc-tests/polybench-gpu/dense-linear-algebra/GESUMMV/gesummv.c
gramschm_SOURCES = openacc-tests/polybench-gpu/dense-linear-algebra/GRAMSCHM/gramschmidt.c
mvt_SOURCES      = openacc-tests/polybench-gpu/dense-linear-algebra/MVT/mvt.c
syrk_SOURCES     = openacc-tests/polybench-gpu/dense-linear-algebra/SYRK/syrk.c
syr2k_SOURCES    = openacc-tests/polybench-gpu/dense-linear-algebra/SYR2K/syr2k.c

#graph problems, expressed as dense matrices, classified as dense linear algebra for now
corr_SOURCES  = openacc-tests/polybench-gpu/dense-linear-algebra/CORR/corr.c
covar_SOURCES = openacc-tests/polybench-gpu/dense-linear-algebra/COVAR/covar.c

nodist_EXTRA_2dconv_SOURCES   = dummy.cpp
nodist_EXTRA_2mm_SOURCES      = dummy.cpp
nodist_EXTRA_3dconv_SOURCES   = dummy.cpp
nodist_EXTRA_3mm_SOURCES      = dummy.cpp
nodist_EXTRA_atax_SOURCES     = dummy.cpp
nodist_EXTRA_bicg_SOURCES     = dummy.cpp
nodist_EXTRA_corr_SOURCES     = dummy.cpp
nodist_EXTRA_covar_SOURCES    = dummy.cpp
nodist_EXTRA_fdtd2d_SOURCES   = dummy.cpp
nodist_EXTRA_gemm_SOURCES     = dummy.cpp
nodist_EXTRA_gesummv_SOURCES  = dummy.cpp
nodist_EXTRA_gramschm_SOURCES = dummy.cpp
nodist_EXTRA_mvt_SOURCES      = dummy.cpp
nodist_EXTRA_syrk_SOURCES     = dummy.cpp
nodist_EXTRA_syr2k_SOURCES    = dummy.cpp

2dconv_LDADD   = libctsar.a
2mm_LDADD      = libctsar.a
3dconv_LDADD   = libctsar.a
3mm_LDADD      = libctsar.a
atax_LDADD     = libctsar.a
bicg_LDADD     = libctsar.a
corr_LDADD     = libctsar.a
covar_LDADD    = libctsar.a
fdtd2d_LDADD   = libctsar.a
gemm_LDADD     = libctsar.a
gesummv_LDADD  = libctsar.a
gramschm_LDADD = libctsar.a
mvt_LDADD      = libctsar.a
syrk_LDADD     = libctsar.a
syr2k_LDADD    = libctsar.a

2dconv_LDFLAGS   = ${ACC_LIBS}
2mm_LDFLAGS      = ${ACC_LIBS}
3dconv_LDFLAGS   = ${ACC_LIBS}
3mm_LDFLAGS      = ${ACC_LIBS}
atax_LDFLAGS     = ${ACC_LIBS}
bicg_LDFLAGS     = ${ACC_LIBS}
corr_LDFLAGS     = ${ACC_LIBS}
covar_LDFLAGS    = ${ACC_LIBS}
fdtd2d_LDFLAGS   = ${ACC_LIBS}
gemm_LDFLAGS     = ${ACC_LIBS}
gesummv_LDFLAGS  = ${ACC_LIBS}
gramschm_LDFLAGS = ${ACC_LIBS}
mvt_LDFLAGS      = ${ACC_LIBS}
syrk_LDFLAGS     = ${ACC_LIBS}
syr2k_LDFLAGS    = ${ACC_LIBS}


2dconv_CPPFLAGS = -I${top_srcdir}/openacc-tests/include  -I${top_srcdir}/openacc-tests/polybench-gpu/common
2mm_CPPFLAGS = -I${top_srcdir}/openacc-tests/include  -I${top_srcdir}/openacc-tests/polybench-gpu/common
3dconv_CPPFLAGS = -I${top_srcdir}/openacc-tests/include  -I${top_srcdir}/openacc-tests/polybench-gpu/common
3mm_CPPFLAGS = -I${top_srcdir}/openacc-tests/include  -I${top_srcdir}/openacc-tests/polybench-gpu/common
atax_CPPFLAGS = -I${top_srcdir}/openacc-tests/include  -I${top_srcdir}/openacc-tests/polybench-gpu/common
bicg_CPPFLAGS = -I${top_srcdir}/openacc-tests/include  -I${top_srcdir}/openacc-tests/polybench-gpu/common
corr_CPPFLAGS = -I${top_srcdir}/openacc-tests/include  -I${top_srcdir}/openacc-tests/polybench-gpu/common
covar_CPPFLAGS = -I${top_srcdir}/openacc-tests/include  -I${top_srcdir}/openacc-tests/polybench-gpu/common
fdtd2d_CPPFLAGS = -I${top_srcdir}/openacc-tests/include  -I${top_srcdir}/openacc-tests/polybench-gpu/common
gemm_CPPFLAGS = -I${top_srcdir}/openacc-tests/include  -I${top_srcdir}/openacc-tests/polybench-gpu/common
gesummv_CPPFLAGS = -I${top_srcdir}/openacc-tests/include  -I${top_srcdir}/openacc-tests/polybench-gpu/common
gramschm_CPPFLAGS = -I${top_srcdir}/openacc-tests/include  -I${top_srcdir}/openacc-tests/polybench-gpu/common
mvt_CPPFLAGS = -I${top_srcdir}/openacc-tests/include  -I${top_srcdir}/openacc-tests/polybench-gpu/common
syrk_CPPFLAGS = -I${top_srcdir}/openacc-tests/include  -I${top_srcdir}/openacc-tests/polybench-gpu/common
syr2k_CPPFLAGS = -I${top_srcdir}/openacc-tests/include  -I${top_srcdir}/openacc-tests/polybench-gpu/common


