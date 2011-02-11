#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += kmeans

kmeans_SOURCES = dense-linear-algebra/rodinia-kmeans/kmeans.c \
	dense-linear-algebra/rodinia-kmeans/cluster.c \
	dense-linear-algebra/rodinia-kmeans/getopt.c \
	dense-linear-algebra/rodinia-kmeans/kmeans_clustering.c \
	dense-linear-algebra/rodinia-kmeans/kmeans_opencl.cpp \
	dense-linear-algebra/rodinia-kmeans/rmse.c

all_local += rodinia-kmeans-all-local
exec_local += rodinia-kmeans-exec-local

rodinia-kmeans-all-local:
	cp $(top_srcdir)/dense-linear-algebra/rodinia-kmeans/kmeans_opencl_kernel.cl .

rodinia-kmeans-exec-local:
	cp $(top_srcdir)/dense-linear-algebra/rodinia-kmeans/kmeans_opencl_kernel.cl ${DESTDIR}${bindir}
