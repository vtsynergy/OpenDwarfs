# OpenCL tree: defects found and fixed

Every change is marked `[FIX-OCL-nn]` in the source next to the code it touches. "Effect" says whether results can change.

| ID | App | Defect | Fix | Effect |
|---|---|---|---|---|
| 01 | crc | `-a` mismatches were only printed to stderr; exit status always 0 | counted, `SUCCESS/FAILURE` printed, reflected in `valid=` | none on correct runs |
| 02 | bfs, lud, kmeans, nw, srad, cfd | `kernel_files = malloc(sizeof(char*)*num_kernels)` allocated 8 bytes for 10-20 character names (heap overflow; glibc aborted with *buffer overflow detected*) | `malloc(256)`; the build no longer needs `-D_FORTIFY_SOURCE=0` | none |
| 03 | crc | kernel include path only valid inside the old autotools tree; `(unsigned char*)` cast to a generic-address-space pointer rejected by NVIDIA's compiler | local header staged next to the kernel; `(__private uchar*)` cast | kernel now compiles on NVIDIA |
| 04 | spmv | `float_array_comp` returned `void`; a failing `-a` check still exited 0 | returns the mismatch count | none |
| 05 | bfs | no correctness check | serial BFS compared with the device costs | validity reporting |
| 06 | nw | read-back sized with `sizeof(float)` for an int buffer | `sizeof(int)` | none |
| 07 | nw | no correctness check | serial fill of the same recurrence compared cell by cell | validity reporting |
| 08 | nw | traceback loop condition `i>=0, j>=0` (comma operator: `i` untested) | `&&` | traceback output could previously run past the matrix edge |
| 09 | srad | no correctness check; disabled CPU path stored the diffusion coefficient in `int` (truncated to 0/1) | independent float host reference added and compared; `float` array | validity reporting |
| 10 | cfd | padded elements' normals were self-assigned (uninitialised) although the comment says "duplicate the last element" | normals of the last element are copied | only the padding elements (not in the output) |
| 11 | cfd | no validity check | solution must be finite with positive density | validity reporting |
| 12 | kmeans | no validity check | centres (and RMSE) must be finite | validity reporting |
| 13 | kmeans | `rms_err` accumulated into a shared variable inside `omp parallel for` without a reduction (latent race; the OpenCL build ignores the pragma) | `reduction(+:sum_euclid)` | none in this build; correct if built with `-fopenmp` |
| 14 | nqueens | no correctness check; error-reporting `if,if,if,else` chain printed the generic message after errors 1 and 2 | compared with the original CPU solver; `else if` | validity reporting |
| 15 | hmm | a local `cl_program _cl_program;` shadowed the global, so `clReleaseProgram` released an uninitialised handle | shadow removed | none |
| 16 | hmm | `est_a_dev` reads `alpha_d`/`beta_d` one row past the allocation | one extra row allocated | removes undefined behaviour in the (unprinted) re-estimated A matrix |
| 17 | hmm | no correctness check | closed-form log-likelihood of the uniform model compared | validity reporting |
| 18 | hmm | `v_model` uninitialised without `-v`; `return;` in `int main` | initialised/diagnosed; proper status | none |
| 19 | swat | appended to `../kernelTime.txt` through an unchecked `fopen` | removed | no stray file in the parent directory |
| 20 | swat | no validity indication | best score must be positive | validity reporting |
| 21 | tdm | record buffers hold `MaxRecords` entries; larger frontiers overflowed them silently | bound check, fails loudly | none for in-range inputs |
| 22 | tdm | no validity indication | episodes mined > 0 | validity reporting |
| 23 | gem | no validity check | potentials must be finite | validity reporting |
| 24 | fft | `-D pow1` passed `fftn/16` instead of log2 (the kernels use it as a shift count) | `log2(fftn)` | none for the sizes in use (the affected group index is 0 either way) |
| 25 | fft | the 128- and 256-point kernels let every row of work-items store, though only row 0 loads data: the others overwrote the result with uninitialised registers | store guarded by `jj == 0` | **128/256-point results were garbage before** |
| 26 | fft | the input came from unseeded `rand()`; OpenCL runtimes may consume `rand()` while loading, so the input differed between implementations/runs | `srand(1)` before generating the data | results now reproducible and comparable |
| 27 | fft | no validation | sorted magnitudes compared with a double-precision DFT of the same input | validity reporting |
| 28 | lud | debug `printf` inside the block loop | removed | less output |
| 29 | lud | `create_matrix_from_random` wrote through an unallocated pointer and assigned instead of summing; `matrix_multiply` assigned instead of accumulating | allocated/accumulated | unused code paths (`-s` is unsupported) |
| 30 | lud | `lud_verify` declared `func_ret_t` but had no `return` | returns status | none |
| 31 | lud | fixed absolute verification tolerance 1e-4 on float data whose rounding error grows with n (n = 2048 failed spuriously) | tolerance scaled with matrix size and magnitude (unchanged at n <= 256) | removes false failures |
| 32 | lud | verification only with `-v`, which also dumps the matrices | quiet check on every run | validity reporting |
| 33 | kmeans | kernels, program, queue and context were released in `deallocateMemory()`, i.e. at the end of every `cluster()` call; the timer's warm-up/repeat loop calls `cluster()` again and `clSetKernelArg` used a released kernel (segfault on PoCL and NVIDIA, found in the dev3 run) | released once after the repeat loop (`releaseCL`) | none in the output (RMSE identical); makes repeated runs possible |
| 34 | spmv, bfs, crc, cfd, kmeans, nqueens, hmm, swat, tdm, gem | work-group sizes came from device queries, the CLI, or NULL (runtime choice) and differed from the OpenMP tree | pinned launch configuration of common/dwarf_launch.h (local 256, global = roundup(N,256), `gid >= n` guards added where missing; algorithmic shapes kept); details in omp/results/LAUNCH-CONFIGS.md. Also: crc kernel gets a `num_items` argument; hmm s_dot/mvm naive kernels grid-stride; hmm `init_ones` buffer sized by nsymbols; kmeansPoint launched 1-D; nqueens thread count fixed to 61440 | work-group sizes only; no algorithm change |
| 35 | nqueens | CPU devices built a different algorithm (`-D FORCE_CPU`, no atomics) than GPUs and than the OpenMP tree | CPU devices build the same persistent-thread atomics kernel (validated on PoCL and Intel, n = 12, 14) | changes the OpenCL CPU kernel; earlier CPU nqueens timings are of the other algorithm (note: the SYCL tree still selects FORCE_CPU on CPU) |
| 36 | kmeans | the RMSE quality metric was inside the timed region (`cluster()`) | excluded with `dwarf_exclude_begin/end` (same in all five trees) | timed work only |

## Common timer
All 14 applications report through `common/dwarf_timer.h` (one `DWARF_RESULT` line: median/mean/min/max of the timed region, `valid=`).
The timed region is the device work: buffer creation (where it can be repeated), host-to-device copies, kernels, device-to-host copies,
including small host steps interleaved with them. Program build/JIT, device discovery and file I/O are outside it. Where the one-time
upload of read-only data cannot be repeated cheaply (tdm), it is outside the region; see the comment in each source file.

## Not fixed (flagged)
* hmm: `mat_vec_mul('N', ...)` in `estimate_b` is a no-op (only `'n'`/`'t'` are handled) and the B matrix is indexed as
  [symbol][state] by kernels that assume the other layout; it affects only the re-estimated model, not the printed log-likelihood.
  Changing it needs the author's intent, so it is left as is and the OpenMP translation reproduces it.
* fft: the 2-D path (`--2D`) was not checked or timed.
