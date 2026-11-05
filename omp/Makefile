include common/benchmarks.mk

# Device to build/run against; passed straight through to each dwarf's
# Makefile (see common/benchmark.mk). Known: cpu h100 h200 a100 mi210 mi300a mi100.
DEVICE ?= cpu

# Space-separated device list used only by `all`, e.g.:
#   make all DEVICES="cpu a100 mi100"
DEVICES ?= cpu

.PHONY: build run all test clean

build:
	for d in $(BENCHMARK_DIRS); do $(MAKE) -C $$d build DEVICE=$(DEVICE) || exit $$?; done

run:
	for d in $(BENCHMARK_DIRS); do $(MAKE) -C $$d run DEVICE=$(DEVICE) || exit $$?; done

all:
	for dev in $(DEVICES); do \
	  for d in $(BENCHMARK_DIRS); do $(MAKE) -C $$d build DEVICE=$$dev || exit $$?; done; \
	done

test:
	$(MAKE) -C test run DEVICE=$(DEVICE)

clean:
	for d in $(BENCHMARK_DIRS); do $(MAKE) -C $$d clean || exit $$?; done
	$(MAKE) -C test clean
