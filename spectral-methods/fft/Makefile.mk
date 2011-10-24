#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

bin_PROGRAMS += clfft

clfft_SOURCES = spectral-methods/fft/src/opencl/fft/fft.cpp spectral-methods/fft/src/opencl/fft/fftlib.cpp spectral-methods/fft/src/opencl/common/main.cpp spectral-methods/fft/src/opencl/common/Event.cpp spectral-methods/fft/src/opencl/common/OpenCLDeviceInfo.cpp spectral-methods/fft/src/opencl/common/OpenCLNodePlatformContainer.cpp spectral-methods/fft/src/opencl/common/OpenCLPlatform.cpp spectral-methods/fft/src/common/OptionParser.cpp spectral-methods/fft/src/common/ResultDatabase.cpp spectral-methods/fft/src/common/Timer.cpp spectral-methods/fft/src/common/Option.cpp spectral-methods/fft/src/common/InvalidArgValue.cpp 

clfft_CPPFLAGS = -I$(top_srcdir)/spectral-methods/fft/src/common -I$(top_srcdir)/spectral-methods/fft/src/opencl/common

all_local += clfft-all-local
exec_local += clfft-exec-local

clfft-all-local:
	cp $(top_srcdir)/spectral-methods/fft/src/opencl/fft/fft.cl .

clfft-exec-local:
	cp $(top_srcdir)/spectral-methods/fft/src/opencl/fft/fft.cl ${DESTDIR}${bindir}
