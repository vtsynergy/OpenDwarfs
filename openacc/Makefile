# OpenACC top level: builds the CoreTSAR-lite support library and then each benchmark.
#   make libctsar | build | clean     (requires nvc: source $ACC_PREFIX/env.sh)
include benchmarks.mk
CC  ?= nvc
CXX ?= nvc++
B := build
.PHONY: libctsar build clean
libctsar: $(B)/libctsar.a
$(B)/libctsar.a: src/ctsar.cpp src/pgi_shim.c
	@mkdir -p $(B)/shim
	@echo '/* empty: stands in for the autoconf config.h */' > $(B)/shim/config.h
	$(CXX) -O2 -std=c++11 -Iinclude -I$(B)/shim -c src/ctsar.cpp -o $(B)/ctsar.o
	$(CC)  -O2 -acc -Dacc_device_radeon=acc_device_not_host -Iinclude -I$(B)/shim -c src/pgi_shim.c -o $(B)/shim.o
	@printf '%s\n' '/* nvc has no OpenCL OpenACC backend: stand-ins for the removed PGI accessors */' \
	  'void* acc_get_opencl_queue(long a){(void)a;return 0;}' \
	  'void* acc_get_current_opencl_context(void){return 0;}' \
	  'void* acc_get_current_opencl_device(void){return 0;}' > $(B)/shim/opencl_stubs.c
	$(CC) -O2 -c $(B)/shim/opencl_stubs.c -o $(B)/opencl_stubs.o
	ar rcs $@ $(B)/ctsar.o $(B)/shim.o $(B)/opencl_stubs.o
build: libctsar
	@for d in $(OPENACC_BENCHMARK_DIRS); do $(MAKE) -C $$d build >$$d/build.log 2>&1 && echo "BUILD OK    $$d" || echo "BUILD FAIL  $$d (see $$d/build.log)"; done
clean:
	rm -rf $(B)
	@for d in $(OPENACC_BENCHMARK_DIRS); do $(MAKE) -C $$d clean; rm -f $$d/build.log; done
