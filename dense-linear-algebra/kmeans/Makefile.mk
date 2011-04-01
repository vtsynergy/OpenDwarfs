#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += kmeans

kmeans_SOURCES = dense-linear-algebra/kmeans/kmeans.c \
	dense-linear-algebra/kmeans/cluster.c \
	dense-linear-algebra/kmeans/getopt.c \
	dense-linear-algebra/kmeans/kmeans_clustering.c \
	dense-linear-algebra/kmeans/kmeans_opencl.cpp \
	dense-linear-algebra/kmeans/rmse.c

all_local += kmeans-all-local
exec_local += kmeans-exec-local

kmeans-all-local:
	cp $(top_srcdir)/dense-linear-algebra/kmeans/kmeans_opencl_kernel.cl .

kmeans-exec-local:
	cp $(top_srcdir)/dense-linear-algebra/kmeans/kmeans_opencl_kernel.cl ${DESTDIR}${bindir}
