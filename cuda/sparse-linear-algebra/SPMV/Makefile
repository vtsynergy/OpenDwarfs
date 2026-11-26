# CUDA build for sparse-linear-algebra/SPMV.  Usage: make [CUDACC=...] [RUN_ARGS=...] {build|run}   (see ../../common.mk)
NAME    := csr
SRCS    := src/csr.cu src-common/sparse_formats.c src-common/ziggurat.c src-common/common.c ../../include/common_args.c
EXTRA_CPPFLAGS := -I$(CURDIR)/inc
EXTRA_CUDAFLAGS := -I$(CURDIR)/src

include ../../common.mk

# Input generator (same program and options as the OpenCL tree)
$(BUILD_DIR)/createcsr: src-test/createcsr.c src-common/sparse_formats.c src-common/ziggurat.c src-common/common.c ../../include/common_args.c
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS_ALL) $(CPPFLAGS_ALL) $^ -o $@ -lm
build: $(BUILD_DIR)/createcsr
