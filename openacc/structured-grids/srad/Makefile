# OpenACC translation of opencl/structured-grids/srad (Rodinia SRAD)
NAME     := srad
SRCS     := srad.c
RUN_ARGS ?= 512 512 0 127 0 127 0.5 2
# make EXTRA_CPPFLAGS=-DOUTPUT prints the final image (as the OpenCL host does with -DOUTPUT)
EXTRA_LIBS := -lm
include ../../dwarf.mk
