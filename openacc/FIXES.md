# OpenACC translations: deliberate differences, validation, audit

Applies to the 14 translations of `../opencl/` (not to the legacy PGI/CoreTSAR code). The corrected behaviour of
`../opencl/FIXES.md` is included (bound checks, validators, `srand(1)`, fft jj==0 stores, log2 `pow1`, pinned launches).

## Deliberate differences from OpenCL
| App | Difference |
|---|---|
| all (GPU) | Work-groups are gangs, work-items vector lanes (`acc-common/dwarf_acc.h`). Free-size kernels: vector length 128, ceil(N/128) gangs; the OpenCL `gid < n` guard of 1-D kernels is the loop bound. |
| all (multicore) | Launch clauses expand to nothing; thread count is `ACC_NUM_CORES`. Data directives are no-ops (host memory is device memory): the timed region has **no H2D/D2H copies**, whereas OpenCL on a CPU device still copies buffers. Host mirrors that the code fills by `memcpy` (kmeans clusters, lud/nw/srad restore, hmm) remain. This matters for any CPU timing comparison. |
| nw, lud, fft | OpenACC has no barrier inside a gang: a barrier-delimited phase is one vector loop. fft kernels are fft.cl split at its barriers; `__local` memory and per-work-item registers are gang-private arrays. fft 128/256: rows jj!=0 start from zero instead of uninitialised registers. `--2D` not translated. |
| swat | No cross-gang barrier: one launch per anti-diagonal (14 gangs x 64) instead of one resident kernel with a lock-based barrier. The trace-back is an `acc routine seq` in a 1x1 parallel region. The `../kernelTime.txt` append is not present. |
| nqueens | Only the OpenCL GPU variant is translated (scalar, `USE_ATOMICS` persistent workers, `WORK_ITEMS=256`, 61440 threads) and **multicore runs it too**. OpenCL on a CPU device builds a different variant (`-D FORCE_CPU`, `kernels_nqueens.cl` lines 82-289: no atomics, one task per work-item). CPU results agree; CPU work distribution is not the same. |
| crc, spmv | Non-blocking enqueues are one `async` queue + `wait` (OpenCL uses three queues with events). |
| cfd | `clEnqueueCopyBuffer` is `acc_memcpy_device`. |
| gem | `native_exp` is `expf`. The OpenCL program also uploads its buffers once before the timed region (no effect on results); this is not repeated. Directory is `n-body-methods/gem-dwarf` because `n-body-methods/gem` is legacy code. |
| tdm | `d_episodeCandidates/Intervals/Support` buffers (never read by the kernels) are not created. A zero-record launch is skipped. |
| kmeans | `rmse.c` keeps its OpenMP pragma but is built without `-mp`: serial, as in OpenCL. |

## Toolchain notes (nvc 25.7)
* `acc serial` miscompile (swat): with the trace-back written in an `acc serial` block, the GPU build walked a wrong path
  (output length about 1e6 instead of about 1200) while the DP matrices hashed equal to multicore. A separate loop
  variable alone did not help; the `routine seq` form is correct.
* nvc `-O3` vectorises floating-point reductions by reassociation. All builds use `-Mvect=noassoc` (`dwarf.mk`).
  Before the flag, hmm's log-likelihood differed from OpenCL in the last digits and far more srad pixels differed.
* CPU instruction set: no `-tp`, so nvc targets the build host (`-tp znver4` here, AMD EPYC 9454; AVX-512, FMA
  contraction on by default). The OpenCL host code is built by gcc/g++ `-O3` without `-march` (x86-64 baseline);
  PoCL compiles kernels for the host CPU at run time. GPU: `-acc=gpu -gpu=cc90`.

## Timed region audit (OpenCL file:line -> OpenACC file:line)
B/E = `dwarf_region_begin/end`. "Outside" = before B. Both trees upload and download inside the region unless noted.

| App | OpenCL: create (outside) ; B-E ; transfers inside | OpenACC: create (outside) ; B-E | Inside/outside differences |
|---|---|---|---|
| lud | lud.c 132 ; 144-222 ; W146 R216 | lud.c 159 ; 164-166 ; update 125/132 | none |
| spmv | csr.c 313-328 ; 385-431 ; W392-404 R425 | csr.c 142 ; 151-158 ; update 58/69 | none |
| bfs | bfs.cpp 166-212 ; 302-350 ; W303-313, stop R, cost R346 | bfs.cpp 100 ; 110-121 ; 111,114,117,120 | none |
| crc | crc_algo.c 343,345 ; 367-393 ; W116 R143 | crc_algo.c 151 ; 155-161 ; 81/84 | none |
| nw | needle.c 187,189 ; 199-270 ; W200,206 R265 | needle.c 197 ; 202-213 ; 203/212 | none |
| srad | srad.c 249-259 ; 282-419 ; W364 R405 per iteration | srad.c 326 ; 331-338 ; 334/336 | none |
| cfd | cfd.cpp 57 (helper), 257-265 (far field, with data), 433 ; 449-553 ; W450-452,465 R552 | cfd.cpp 502 (far field copyin) ; 511-524 ; 517,521,523 | none |
| gem | calculate_potential.cpp 339-374 ; 531-675 ; W532-596 R642 | calculate_potential.cpp 321 ; 328-339 ; 329/338 | none (pre-region first upload 377-441 not repeated) |
| tdm | tdm_ocl.cpp 192-207, events upload 209,215 ; 314-432 ; W319,382 R371,421 | tdm_ocl.cpp 277 (create + events copyin) ; 282-323 ; 288,299,304,317 | none |

Other apps: kmeans, nqueens, hmm allocate and release inside the region in both trees (kmeans inside `cluster()`,
nqueens inside `Compute()`, hmm inside `run_hmm_bwa` and per dot product). fft allocates inside and frees outside in
both. swat creates all buffers and uploads BLOSUM outside in both. Eager (`enter data`) versus the OpenCL runtime's
possibly lazy allocation was not examined.

## Validation (final binaries, NVIDIA H100; multicore with 32 cores)
Own validator: `./run_all.sh --target both --size small` -> 28/28 PASS. Against OpenCL (NVIDIA GPU runtime), same
results for GPU and multicore builds:
* identical: crc (all page CRCs), spmv (printed vector, 2 decimals), bfs (cost file), nw (traceback line; all cells also
  match the serial fill), kmeans (centroid table and RMSE), nqueens (N=8,12,14), hmm (log-likelihood, 6 decimals),
  swat (alignment text, 1K1/3K1/5K1), tdm (result file).
* numeric, none above relative 1e-3 (floating-point contraction / math library): srad max 9.3e-6, cfd 7.2e-6,
  gem 7.3e-5, fft sizes 1-7 4.6e-5.
* lud: default-build LU factors do **not** match OpenCL within 1e-3. The 256 input is ill-conditioned (double-precision
  reference: smallest pivot 0.0078, growth 434; measured before `-Mvect=noassoc`, median relative error against it was
  0.057 OpenCL, 0.040 OpenACC GPU, 0.035 multicore). With strict arithmetic on both sides (OpenCL: `FP_CONTRACT OFF`,
  correctly rounded divide; OpenACC GPU: `-Kieee -Mnofma -gpu=nofma`; multicore: `-O2 -Kieee -Mnofma -Mnovect`) the LU
  matrices are bit-identical at 256 and 2048. The application's L*U check passes for the default builds.

Outputs not printed by the stock OpenCL binaries (srad, fft, strict lud) came from a modified copy of `../opencl`
built in a scratch directory (`-DOUTPUT`, `PRINT_RESULT`, the strict flags). nqueens, hmm, nw etc. used copies of the
prebuilt binaries. No performance measurements were made.

## Not done
* Legacy benchmarks that fail to build were not ported.
