# OpenCL build for combinational-logic/crc (replaces the autotools Makefile.mk).
# Usage: make [OCL_VENDOR=nvidia|pocl] [RUN_ARGS=...] {build|run|clean}   (see ../../common.mk)
NAME    := crc
SRCS    := src/crc_algo.c src-common/crc_formats.c
KERNELS := src/crc_kernel.cl inc/eth_crc32_lut.h

include ../../common.mk

# Input generator (needed to make CRC inputs)
$(BUILD_DIR)/createcrc: src-test/createcrc.c src-common/crc_formats.c $(COMMON_SRCS)
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS_ALL) $(CPPFLAGS_ALL) $^ -o $@ $(LDLIBS_ALL)
build: $(BUILD_DIR)/createcrc
