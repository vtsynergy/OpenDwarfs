# OpenACC translation of opencl/spectral-methods/fft (SHOC FFT, 1-D path)
# make EXTRA_CPPFLAGS=-DPRINT_RESULT prints the transformed data (as the OpenCL host with -DPRINT_RESULT)
NAME     := fft
SRCS     := fft.cpp
RUN_ARGS ?= --pts 1 --passes 5
EXTRA_LIBS := -lm
include ../../dwarf.mk
