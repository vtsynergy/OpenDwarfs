#
# Copyright 2010 by Virginia Polytechnic Institute and State
# University. All rights reserved. Virginia Polytechnic Institute and
# State University (Virginia Tech) owns the software and its
# associated documentation.
#

AUTOMAKE_OPTIONS = subdir-objects

ACLOCAL_AMFLAGS = -I m4
AM_CPPFLAGS = -Werror

bin_PROGRAMS =
all_local =
exec_local =
testruns =

platform_n = 0
device_n = 0
device_t = 0

bin_PROGRAMS += ocd
ocd_SOURCES = include/common_args.c include/rdtsc.c
ocd_LINK = true

noinst_LIBRARIES = libopts.a
libopts_a_SOURCES = opts/opts.c 

LDADD = include/common_args.o include/rdtsc.o libopts.a

if build_gem
include n-body-methods/gem/Makefile.mk
testruns += gem
endif

if build_nw
include dynamic-programming/nw/Makefile.mk
testruns += nw
endif

if build_swat
include dynamic-programming/swat/Makefile.mk
testruns += swat
endif

if build_srad
include structured-grids/srad/Makefile.mk
testruns += srad
endif

if build_lud
include dense-linear-algebra/lud/Makefile.mk
testruns += lud
endif

if build_kmeans
include dense-linear-algebra/kmeans/Makefile.mk
testruns += kmeans
endif

if build_samplecl
include samplecl/Makefile.mk
testruns += samplecl
endif

if build_fft
include spectral-methods/fft/Makefile.mk
testruns += fft
endif

if build_csr
include sparse-linear-algebra/SPMV/Makefile.mk
testruns += csr
endif

if build_crc
include combinational-logic/crc/Makefile.mk
testruns += crc
endif

if build_cfd
include unstructured-grids/cfd/Makefile.mk
testruns += cfd
endif

if build_bfs
include graph-traversal/bfs/Makefile.mk
testruns += bfs
endif

#if build_bsort
#include graph-traversal/bitonic-sort/Makefile.mk
#endif

#if build_oesort
#include graph-traversal/oddeven-sort/Makefile.mk
#endif

if build_tdm
include finite-state-machine/tdm/Makefile.mk
testruns += tdm
endif

if build_astar
include branch-and-bound/astar/Makefile.mk
testruns += astar
endif

if build_astar
include branch-and-bound/nqueens/Makefile.mk
testruns += nqueens
endif

if build_bwa_hmm
include graphical-models/hmm/Makefile.mk
testruns += bwa_hmm
endif

all-local: $(all_local)

install-exec-local: $(exec_local)

test-run:
	@if [ $(device_t) -gt 3 ]; then \
	##NO CHECKING FOR PLATFORM_ID and DEVICE_ID.
	echo "Wrong device type number. Enter make run-test platform_n=<PLATFORM_ID> device_n=<DEVICE_ID> device_t=<DEVICE_TYPE>, where:"; \
	echo "PLATFORM_ID and DEVICE_ID must be valid OpenCL platform and device IDs and DEVICE_TYPE can be:"; \
	echo "0 => CPU"; \
	echo "1 => GPU"; \
	echo "2 => MIC"; \
	echo "3 => FPGA"; \
	echo "(Values for omitted options are set to zero.)"; \
	exit 9; \
	fi
	@echo "Running tests for the following dwarfs: $(testruns)";
	@$(foreach dwarf_test,$(testruns),$(call run_dwarf_test, $(dwarf_test)))

define run_dwarf_test
	if [ $(1) = "gem" ]; then \
	echo "Running gem..."; \
	./gemnoui -p $(platform_n) -d $(device_n) -t $(device_t) -- ../test/n-body-methods/gem/nucleosome 80 1 0; \
	elif [ $(1) = "nw" ]; then \
	echo "Running nw..."; \
	./needle -p $(platform_n) -d $(device_n) -t $(device_t) -- 2048 10; \
	elif [ $(1) = "swat" ]; then \
	echo "Running swat..."; \
	./swat -p $(platform_n) -d $(device_n) -t $(device_t) -- ../test/dynamic-programming/swat/query1K1 ../test/dynamic-programming/swat/sampledb1K1; \
	elif [ $(1) = "srad" ]; then \
	echo "Running srad..."; \
	./srad -p $(platform_n) -d $(device_n) -t $(device_t) -- 256 256 0 127 0 127 0.5 2; \
	elif [ $(1) = "lud" ]; then \
	echo "Running lud..."; \
	./lud -p $(platform_n) -d $(device_n) -t $(device_t) -- -i ../test/dense-linear-algebra/lud/512.dat; \
	elif [ $(1) = "kmeans" ]; then \
	echo "Running kmeans..."; \
	./kmeans -p $(platform_n) -d $(device_n) -t $(device_t) -- -i ../test/dense-linear-algebra/kmeans/204800.txt; \
	elif [ $(1) = "fft" ]; then \
	echo "Running fft..."; \
	./clfft -p $(platform_n) -d $(device_n) -t $(device_t) -- --pts 1; \
	##elif [ $(1) = "csr" ]; then \
	##echo "Running csr..."; \
	##./csr -p $(platform_n) -d $(device_n) -t $(device_t) -- -i ../test/sparse-linear-algebra/SPMV/csr_65536.txt; \
	##elif [ $(1) = "crc" ]; then \
	##echo "Running crc..."; \
	##./crc -p $(platform_n) -d $(device_n) -t $(device_t) --; \
	elif [ $(1) = "cfd" ]; then \
	echo "Running cfd..."; \
	./cfd -p $(platform_n) -d $(device_n) -t $(device_t) -- ../test/unstructured-grids/cfd/fvcorr.domn.097K; \
	elif [ $(1) = "bfs" ]; then \
	echo "Running bfs..."; \
	./bfs -p $(platform_n) -d $(device_n) -t $(device_t) -- ../test/graph-traversal/bfs/graph65536.txt; \
	elif [ $(1) = "tdm" ]; then \
	echo "Running tdm..."; \
	./tdm -p $(platform_n) -d $(device_n) -t $(device_t) -- ../test/finite-state-machine/tdm/stream-new-1.csv ../test/finite-state-machine/tdm/streamIntervals.txt 100 0.5 r s n; \
	##elif [ $(1) = "astar" ]; then \
	##echo "Running astar..."; \
	##./astar -p $(platform_n) -d $(device_n) -t $(device_t) --; \
	elif [ $(1) = "nqueens" ]; then \
	echo "Running nqueens..."; \
	./nqueens 16; \
	elif [ $(1) = "bwa_hmm" ]; then \
	echo "Running bwa_hmm..."; \
	./bwa_hmm -n 30 -v n; \
	fi

endef

clean-local:
	rm -f *.cl
