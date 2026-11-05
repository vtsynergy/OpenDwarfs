# OpenCL build for spectral-methods/fft (replaces the autotools Makefile.mk).
# Usage: make [OCL_VENDOR=nvidia|pocl] [RUN_ARGS=...] {build|run|clean}   (see ../../common.mk)
NAME    := clfft
SRCS    := src/opencl/fft/fft.cpp src/opencl/fft/fftlib.cpp src/opencl/common/main.cpp src/opencl/common/Event.cpp src/opencl/common/OpenCLDeviceInfo.cpp src/opencl/common/OpenCLNodePlatformContainer.cpp src/opencl/common/OpenCLPlatform.cpp src/common/OptionParser.cpp src/common/ResultDatabase.cpp src/common/Timer.cpp src/common/Option.cpp src/common/InvalidArgValue.cpp
KERNELS := src/opencl/fft/fft.cl
EXTRA_CPPFLAGS := -I$(CURDIR)/src/common -I$(CURDIR)/src/opencl/common
# Pre-include libstdc++ headers: they must come before the min/max macros defined by the app headers
EXTRA_CXXFLAGS := -include exception -include iostream -include fstream -include sstream -include string -include vector -include algorithm -include cstdlib -include cmath

include ../../common.mk
