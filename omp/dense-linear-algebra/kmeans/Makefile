# OpenMP-offload translation of opencl/dense-linear-algebra/kmeans (Rodinia k-means, OpenDwarfs variant)
BENCHMARK := kmeans
ROOT_DIR := ../..
SRCS := kmeans.c cluster.c kmeans_clustering.c rmse.c kmeans_omp.c
RUN_INPUT ?= $(ROOT_DIR)/test/dense-linear-algebra/kmeans/204800.txt
APP_ARGS = -i $(RUN_INPUT) -r -o
EXTRA_CFLAGS := -std=gnu99 -D_GNU_SOURCE
EXTRA_LIBS := -lm
include $(ROOT_DIR)/common/benchmark.mk
