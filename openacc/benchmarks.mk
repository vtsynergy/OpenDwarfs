# The 14 OpenACC translations of the OpenCL dwarfs (rules: dwarf.mk)
DWARF_BENCHMARK_DIRS := \
	dense-linear-algebra/lud \
	sparse-linear-algebra/spmv \
	graph-traversal/bfs \
	combinational-logic/crc \
	dynamic-programming/nw \
	structured-grids/srad \
	unstructured-grids/cfd \
	dense-linear-algebra/kmeans \
	branch-and-bound/nqueens \
	graphical-models/hmm \
	dynamic-programming/swat \
	finite-state-machine/tdm \
	n-body-methods/gem-dwarf \
	spectral-methods/fft

# Legacy 2014 PGI/CoreTSAR sources (rules: common.mk)
OPENACC_BENCHMARK_DIRS := \
	polybench-gpu/dense-linear-algebra/2MM \
	polybench-gpu/dense-linear-algebra/3MM \
	polybench-gpu/dense-linear-algebra/ATAX \
	polybench-gpu/dense-linear-algebra/BICG \
	polybench-gpu/dense-linear-algebra/CORR \
	polybench-gpu/dense-linear-algebra/COVAR \
	polybench-gpu/dense-linear-algebra/GEMM \
	polybench-gpu/dense-linear-algebra/GESUMMV \
	polybench-gpu/dense-linear-algebra/GRAMSCHM \
	polybench-gpu/dense-linear-algebra/MVT \
	polybench-gpu/dense-linear-algebra/SYR2K \
	polybench-gpu/dense-linear-algebra/SYRK \
	polybench-gpu/structured-grid/2DCONV \
	polybench-gpu/structured-grid/3DCONV \
	polybench-gpu/structured-grid/FDTD-2D \
	structured-grid \
	n-body-methods/gem
