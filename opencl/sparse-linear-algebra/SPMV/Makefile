# OpenCL build for sparse-linear-algebra/SPMV (replaces the autotools Makefile.mk).
# Usage: make [OCL_VENDOR=nvidia|pocl] [RUN_ARGS=...] {build|run|clean}   (see ../../common.mk)
NAME    := csr
SRCS    := src/csr.c src-common/sparse_formats.c src-common/ziggurat.c src-common/common.c
KERNELS := src/spmv_kernel.cl

include ../../common.mk

$(BUILD_DIR)/createcsr: src-test/createcsr.c src-common/sparse_formats.c src-common/ziggurat.c src-common/common.c $(COMMON_SRCS)
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS_ALL) $(CPPFLAGS_ALL) $^ -o $@ $(LDLIBS_ALL)
build: $(BUILD_DIR)/createcsr
