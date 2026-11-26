# OpenACC translation of opencl/dense-linear-algebra/kmeans (Rodinia k-means, OpenDwarfs variant)
NAME     := kmeans
SRCS     := kmeans.c cluster.c kmeans_clustering.c rmse.c kmeans_acc.c
RUN_ARGS ?= -i $(TEST_DIR)/dense-linear-algebra/kmeans/204800.txt -r -o
EXTRA_CFLAGS := -D_GNU_SOURCE
EXTRA_LIBS := -lm
include ../../dwarf.mk
