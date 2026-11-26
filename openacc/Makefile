# OpenACC top level.
#   make build [ACC_TARGET=gpu|multicore]   the 14 OpenCL translations
#   make clean                              their build directories
#   make libctsar | legacy-build | legacy-clean   the 2014 PGI/CoreTSAR code (source $ACC_PREFIX/env.sh first)
include benchmarks.mk
CC  ?= nvc
CXX ?= nvc++
ACC_TARGET ?= gpu
.PHONY: build clean libctsar legacy-build legacy-clean

build:
	@fail=0; for d in $(DWARF_BENCHMARK_DIRS); do \
	  if $(MAKE) -s -C $$d build ACC_TARGET=$(ACC_TARGET) >$$d/build.log 2>&1; then echo "BUILD OK    $$d"; \
	  else echo "BUILD FAIL  $$d (see $$d/build.log)"; fail=1; fi; done; exit $$fail
clean:
	@for d in $(DWARF_BENCHMARK_DIRS); do $(MAKE) -s -C $$d clean && rm -f ./$$d/build.log; done

libctsar: build/libctsar.a
build/libctsar.a: src/ctsar.cpp src/pgi_shim.c
	@mkdir -p build/shim
	@echo '/* empty: stands in for the autoconf config.h */' > build/shim/config.h
	$(CXX) -O2 -std=c++11 -Iinclude -Ibuild/shim -c src/ctsar.cpp -o build/ctsar.o
	$(CC)  -O2 -acc -Dacc_device_radeon=acc_device_not_host -Iinclude -Ibuild/shim -c src/pgi_shim.c -o build/shim.o
	@printf '%s\n' '/* nvc has no OpenCL OpenACC backend: stand-ins for the removed PGI accessors */' \
	  'void* acc_get_opencl_queue(long a){(void)a;return 0;}' \
	  'void* acc_get_current_opencl_context(void){return 0;}' \
	  'void* acc_get_current_opencl_device(void){return 0;}' > build/shim/opencl_stubs.c
	$(CC) -O2 -c build/shim/opencl_stubs.c -o build/opencl_stubs.o
	ar rcs $@ build/ctsar.o build/shim.o build/opencl_stubs.o
legacy-build: libctsar
	@for d in $(OPENACC_BENCHMARK_DIRS); do $(MAKE) -C $$d build >$$d/build.log 2>&1 && echo "BUILD OK    $$d" || echo "BUILD FAIL  $$d (see $$d/build.log)"; done
legacy-clean:
	rm -rf ./build
	@for d in $(OPENACC_BENCHMARK_DIRS); do $(MAKE) -C $$d clean; rm -f ./$$d/build.log; done
