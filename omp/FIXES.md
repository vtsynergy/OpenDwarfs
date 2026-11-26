# OpenMP tree: defects found and fixed in the translations

The 14 applications that exist in both trees are direct translations of the OpenCL applications (same algorithm, inputs, command line,
validation); the five OpenMP-only applications (atax, cg, gemm, pagerank, fdtd2d) are unchanged. Fixes are marked `[FIX-OMP-nn]` in
the source. Defects that the translations inherited from the OpenCL code are fixed in both trees; the matching OpenCL entries are in
`../opencl/FIXES.md`.

| ID | App | Change | Mirrors |
|---|---|---|---|
| 01 | lud | verification uses the size-scaled tolerance from `common.c` (the local fixed 1e-4 copy failed spuriously for large matrices) | OCL-31 |
| 02 | cfd | padded elements duplicate the last element's normals (were zero) | OCL-10 |
| 03 | tdm | `MaxRecords` bound check | OCL-21 |
| 04 | nw | traceback loop `i>=0 && j>=0` | OCL-08 |
| 05 | fft | 128/256-point kernels: only row `jj == 0` stores | OCL-25 |
| 06 | fft | `pow1 = log2(fftn)` | OCL-24 |
| 07 | fft | explicit `srand(1)` | OCL-26 |
| 08 | kmeans | `reduction(+:sum_euclid)` in `rms_err` (this tree builds the file with `-fopenmp`, where the race is real) | OCL-13 |
| 09 | hmm | one extra row for `alpha`/`beta` (translation did this from the start) | OCL-16 |
| 10 | lud | host (`DEVICE=cpu`) path: a work-group was emulated by a team of real threads with barriers (66% of the time in barrier spin-waits; n=512: 1056 ms on 48 threads vs 36 ms on 1) | on a host one thread runs a work-group in barrier order and the threads are spread over work-groups (`#if DWARF_USE_TARGET` keeps the GPU path unchanged) | identical arithmetic and operation order per work-group; CPU only |
| 11 | nw | same barrier/oversubscription pathology on the host (n=1024: 320 ms on 48 threads vs 22 ms on 1) | same host path: one thread per tile, tiles of an anti-diagonal in parallel | CPU only |
| 12 | fft | clang's host runtime grants a team min(requested, OMP_NUM_THREADS) threads, so the 128-thread barrier kernels were BLOCKED_RUNTIME at the dev3 setting of 48 threads | none in the source: with the pinned 128 CPU threads (all configurations) the kernels run and validate on clang and nvhpc; pts 5-7 (256 threads per team) stay blocked (libomp caps a team at 192) | run-time configuration only |
| 13 | all GPU kernels | flat `target teams distribute parallel for` kernels had no `num_teams`/`thread_limit`, so each compiler chose its own launch shape (confound for the clang vs nvhpc and OpenMP vs OpenCL GPU comparison) | `DWARF_LAUNCH(n)` (common/dwarf_launch.h): 256-thread teams, ceil(n/256) teams, GPU builds only; algorithmic shapes kept (lud, nw, fft, swat k_match, nqueens, srad, hmm 16x16 kernels); table in omp/results/LAUNCH-CONFIGS.md | CPU code unchanged |
| 14 | lud, spmv, bfs, crc, nw, srad, cfd, gem, tdm | device memory was mapped/allocated INSIDE the timed region on every repetition, whereas OpenCL creates its buffers (and in tdm uploads events/times, in cfd the far-field values) outside it | `enter data` once outside the repeat loop, `target update` inside at the points where OpenCL writes/reads, `exit data` after the loop | timed work only; outputs checked against OpenCL |
| 15 | all (clang, DEVICE=cpu) | no `-march`: baseline x86-64 (no AVX, AVX-512 or FMA) while nvc++ and the OpenCL CPU runtimes target the host CPU | `-march=native` for clang CPU builds only (GPU-build host code is still baseline) | CPU speed and last-digit rounding; all CPU speedups measured earlier against omp-clang are invalid |
| 16 | kmeans | `rms_err` ran in parallel (OpenMP) inside the timed region, serially in OpenCL | pragma removed; the call is excluded from the timed region in every tree (`dwarf_exclude_begin/end`) | RMSE now equals OpenCL's serial value |

## Translation notes (deliberate differences, not defects)
* OpenCL work-groups/work-items become OpenMP teams/threads. OpenMP only guarantees upper bounds on both, so loops are strided over the
  teams/threads actually granted; kernels that keep private state across barriers (fft) need exactly one thread per work-item and the
  run is reported as failed on toolchains that cannot provide it (observed: `nvc++` on the V100 at 2048 points; Clang on x86 from 2048 points).
* swat: the OpenCL lock-based global barrier is replaced by one kernel launch per anti-diagonal.
* Constructs that miscompile on one of the two toolchains were rewritten in a form valid for both (no scalar `firstprivate` combined with
  array sections; plain functions instead of lambdas around OpenMP array sections; flat loops for kernels without barriers).
