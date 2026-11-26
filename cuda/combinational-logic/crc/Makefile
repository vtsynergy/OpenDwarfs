# CUDA build for combinational-logic/crc.  Usage: make [CUDACC=...] [RUN_ARGS=...] {build|run}   (see ../../common.mk)
NAME    := crc
SRCS    := src/crc_algo.cu src-common/crc_formats.c ../../include/common_args.c
EXTRA_CPPFLAGS := -I$(CURDIR)/inc
EXTRA_CUDAFLAGS := -I$(CURDIR)/src

include ../../common.mk

# Input generator (same program and options as the OpenCL tree)
$(BUILD_DIR)/createcrc: src-test/createcrc.c src-common/crc_formats.c ../../include/common_args.c
	@mkdir -p $(BUILD_DIR)
	$(CC) $(CFLAGS_ALL) $(CPPFLAGS_ALL) $^ -o $@ -lm
build: $(BUILD_DIR)/createcrc
