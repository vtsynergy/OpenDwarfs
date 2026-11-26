# SYCL build for spectral-methods/fft (1-D only; see ../../common.mk)
NAME      := clfft
SRCS      := src/sycl/fft_sycl.cpp src/sycl/main.cpp
HOST_SRCS := src/common/OptionParser.cpp src/common/Option.cpp src/common/InvalidArgValue.cpp
EXTRA_CPPFLAGS := -I$(CURDIR)/src/common -I$(CURDIR)/src/sycl
include ../../common.mk
