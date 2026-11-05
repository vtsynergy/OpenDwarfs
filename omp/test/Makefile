include ../common/benchmarks.mk

# Smoke tests for all 19 apps, plus real-input cases (see ../README.md).
# Require a successful exit and a DWARF_RESULT line reporting valid=1.
#
# Usage:
#   make -C test run DEVICE=cpu
#   make -C test run DEVICE=a100

DEVICE ?= cpu
BUILD_DIR ?= build

.PHONY: run clean

run:
	python3 -m unittest discover -s . -p 'test_*.py'
	@root_dir="$$(cd .. && pwd)"; \
	fail=0; \
	run_case() { \
	  path="$$1"; shift; \
	  echo "=== $$path (DEVICE=$(DEVICE)) ==="; \
	  if out="$$($(MAKE) -C "$${root_dir}/$${path}" run DEVICE=$(DEVICE) RUN_WARMUP=0 RUN_REPEATS=1 "$$@" 2>&1)" && \
	     printf '%s\n' "$$out" | grep -Eq '^DWARF_RESULT,.*valid=1(,|$$)'; then \
	    echo "PASS: $$path"; \
	  else \
	    echo "FAIL: $$path"; echo "$$out" | tail -15; fail=1; \
	  fi; \
	}; \
	for path in $(BENCHMARK_DIRS); do \
	  run_case "$$path" RUN_SIZE=8 RUN_SECONDARY_SIZE=8 RUN_ITERATIONS=2; \
	done; \
	run_case dense-linear-algebra/lud RUN_INPUT=$$(pwd)/dense-linear-algebra/lud/64.dat; \
	run_case dense-linear-algebra/kmeans RUN_INPUT=$$(pwd)/dense-linear-algebra/kmeans/16384.txt; \
	run_case dynamic-programming/swat RUN_INPUT=$$(pwd)/dynamic-programming/swat/query1K1 RUN_INPUT2=$$(pwd)/dynamic-programming/swat/sampledb1K1; \
	run_case finite-state-machine/tdm RUN_INPUT=$$(pwd)/finite-state-machine/tdm/sim-64-size20.csv RUN_INPUT2=$$(pwd)/finite-state-machine/tdm/ivl.txt RUN_INPUT3=$$(pwd)/finite-state-machine/tdm/30-episodes.txt; \
	run_case graph-traversal/bfs RUN_INPUT=$$(pwd)/graph-traversal/bfs/graph4096.txt; \
	run_case n-body-methods/gem RUN_INPUT=$$(pwd)/n-body-methods/gem/1uwo_A.pqr RUN_INPUT2=$$(pwd)/n-body-methods/gem/1uwo_A.vert; \
	run_case unstructured-grids/cfd RUN_INPUT=$$(pwd)/unstructured-grids/cfd/fvcorr.domn.097K; \
	exit $$fail
	python3 check_binaries.py --device $(DEVICE) --build-dir $(BUILD_DIR)

clean:
	@:
