# Top-level OpenCL driver: build/run/clean every benchmark.
#   make build                     # all
#   make build-<dwarf-dir>         # e.g. make -C dense-linear-algebra/lud
#   make run OCL_VENDOR=pocl
include benchmarks.mk
.PHONY: build run clean
build:
	for d in $(OPENCL_BENCHMARK_DIRS); do $(MAKE) -C $$d build || echo "BUILD FAILED: $$d"; done
clean:
	for d in $(OPENCL_BENCHMARK_DIRS); do $(MAKE) -C $$d clean; done
