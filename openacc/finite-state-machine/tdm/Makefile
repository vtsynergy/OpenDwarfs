# OpenACC translation of opencl/finite-state-machine/tdm (GPU temporal data mining)
NAME     := tdm
SRCS     := tdm_ocl.cpp dataio.cpp
T        := $(TEST_DIR)/finite-state-machine/tdm
RUN_ARGS ?= $(T)/sim-64-size50.csv $(T)/ivl.txt $(T)/30-episodes.txt 128
include ../../dwarf.mk
