# SYCL tree: deliberate differences from the OpenCL tree

Algorithms, inputs, work-group shapes (`common/dwarf_launch.h` and the algorithmic shapes), validators and the timed
region follow `opencl/`, including the corrections in `opencl/FIXES.md`.

## All applications
* **Memory:** USM device allocations and explicit `memcpy` on one in-order queue stand in for `cl_mem` and the command
  queue(s). `spmv` and `crc` use three OpenCL queues; here they use the one in-order queue, waited on at the end.
  `sycl::malloc_device` may commit memory at allocation, while `clCreateBuffer` may defer it to first use. Where the
  OpenCL code allocates outside the timed region and does not upload there (lud, spmv, nw, srad, crc), the first
  repetition can therefore include a one-off cost in OpenCL that SYCL does not have (or has elsewhere). Use warm-ups.
* **Kernels** are compiled ahead of time into the binary (no `.cl` files at run time). `dsycl::prebuild` builds the
  whole kernel bundle before the timed region, in place of `clBuildProgram`. NVIDIA code is compiled for `sm_90`; CPU
  code is SPIR-V, compiled at run time by the Intel CPU OpenCL runtime.
* **Build-time `-D` options** are fixed in the source: lud and srad `BLOCK_SIZE=16` (the OpenCL halving for devices
  with a limit below 256 is dropped, as is the corresponding kmeans loop), fft sizes are template parameters, nqueens
  preprocessor branches are resolved.
* **Floating point:** `-fp-model=precise -fno-offload-fp32-prec-div -fno-offload-fp32-prec-sqrt` (icpx defaults to
  fast-math). This corresponds to the OpenCL default: contraction on, no value-unsafe optimisation.
  fft adds nothing for the OpenCL `-cl-mad-enable`.
* **Host compilers and instruction sets:** no `-march`, `-x`, `-mavx*` or `-mtune` flag anywhere. Host code targets
  baseline x86-64 (no AVX2, no FMA) for both gcc/g++ (`-O3`) and icpx (`-O3 -fp-model=precise`). No OpenMP flag is used.
  Unmodified host files are mostly built with gcc/g++ as in `opencl/` (kmeans, lud `common.c`, spmv, crc and gem host
  files). **Translated host files, and the tdm, nqueens, swat and fft-option-parser helpers, differ as follows:** all
  files that contain SYCL are built by icpx (OpenCL: gcc `-std=gnu89` or g++ `-std=gnu++14`); so are
  `nqueen_cpu.cpp`, tdm `dataio.cpp` and the swat helper `.cpp` files (the fft option parser uses g++).
  Some host work **inside the timed region** is therefore built by a different compiler: the srad ROI statistics, the hmm
  reductions, the kmeans centre accumulation in `kmeans_sycl.cpp`, the nqueens host enumeration, and the tdm prefix sums
  and conflict elimination. Account for this in toolchain comparisons.
* **C to C++** (lud, nw, srad, spmv, crc, hmm drivers): C's double-precision math is kept explicitly where C++ would pick
  a float overload (`exp` in srad input generation, `log10` in the hmm log-likelihood).
* **Device options:** `-p` is ignored. `-o` (alternative `.cl` kernels), spmv/crc `-k`, nqueens `-local` and fft `--2D` /
  `--infoDevices` stop with a message. `DWARF_DEVICE` is set by the program to the device actually used.
* **Validation dumps:** with `DWARF_DUMP_DIR` set, read-backs are written to a file, in most apps from inside the timed
  region. Never set it for timing runs.

## Timed region (allocation and transfers)
Same placement as OpenCL. Line numbers: OpenCL source | SYCL source.

| App | Outside the region | Inside the region |
|---|---|---|
| lud | buffer alloc (132 \| 159) | write, kernels, read (146-216 \| 173-186) |
| spmv | buffer alloc (313-328 \| 265-269) | 5 writes, kernel, read (392-425 \| 300-311) |
| bfs | alloc (166-212 \| 167-180); first 6 uploads (180-213 \| 183-188) | 6 uploads, level loop, read (303-346 \| 208-230) |
| crc | per-block alloc (343,345 \| 335,336) | write, kernel, read per block (367-393 \| 359-367) |
| nw | alloc (187,189 \| 194,195) | 2 writes, kernels, read (200-265 \| 206-218) |
| srad | alloc (249-259 \| 364-369) | ROI sums, write, 2 kernels, read per iteration (364-405 \| 466-471) |
| cfd | alloc, far-field/mesh uploads, 3 initialisations, step-factor zeroing (257-433 \| 455-556) | mesh uploads, initialisation, solver, read (450-552 \| 569-595) |
| gem | alloc (358-374 \| 434-450); first 17 uploads (377-441 \| 453-469) | 17 uploads, kernel loop, read (532-656 \| 485-520) |
| tdm | alloc (192-207 \| 274-289); events, times, candidates, intervals upload (209-227 \| 291-297) | per episode: start records, kernels, count/offset and record transfers (319-421 \| 373-436) |

kmeans (per `cluster()` call), hmm, nqueens and fft allocate **inside** the region in both trees. swat allocates and
uploads the scoring matrix outside it in both. hmm: the per-call `dot_production` scratch buffer is never released by
OpenCL; here it is released after the region ends.

## Per application
* **spmv:** `srand(1)` before generating x and y. This is the unseeded sequence, pinned so a runtime library that calls
  `rand()` cannot shift it.
* **nqueens:** the kernel variant follows the device exactly as the OpenCL host does. **CPU:** the `FORCE_CPU` kernels
  (`cpu_nqueen`, `cpu_nqueen1`), no atomics. **GPU:** the scalar persistent-thread kernels with `USE_ATOMICS`
  (`gpu_nqueen<true>`, `gpu_nqueen1<true>`; `-noatomics` selects `<false>`, not tested). 256 per group, 61440 threads.
  The vectorised kernels are not translated (not taken on NVIDIA or CPU); `-novec`/`-vec4` are ignored.
* **swat:** the cross-work-group spin barrier uses `atomic_ref` (acq_rel, device scope) and `atomic_fence` instead of
  OpenCL 1.2 `atom_add` + `mem_fence`. As in OpenCL, it assumes all work-groups run concurrently, which neither model
  guarantees; it completed in every run here. The block count follows `max_compute_units` as in OpenCL, so a CPU
  restricted to 14 or fewer cores gets fewer blocks (not tested).
  swat also appends its own timing summary to `../runtime.txt` (here `sycl/dynamic-programming/swat/runtime.txt`), as the
  OpenCL program does. It is run output; do not commit it.
* **cfd:** the far-field `float3` values use the 12-byte host struct on both sides (OpenCL reads them as `float3`).
* **hmm:** the flagged no-op `mat_vec_mul('N', ...)` is reproduced. The `c++11-narrowing` diagnostic is disabled for this file.
* **gem:** `native_exp` -> `sycl::native::exp`.
* **fft:** 1-D only. The 2-D kernels are not translated (that path is unverified in OpenCL). Kernel text is generated
  from `fft.cl`; `M_PI` stays a double constant, so the kernels need fp64. Host buffers use `malloc_host`.
* **kmeans:** `rms_err` (`rmse.c`) is serial, as in OpenCL: the OpenMP pragma is not enabled (gcc without `-fopenmp`;
  neither binary references OpenMP symbols).

## Maintenance
Several sources (nqueens and fft kernel headers, most host files) were derived mechanically from the OpenCL sources by
one-off scripts that are not part of the tree. The checked-in files are now the source; a change in `opencl/` must be
ported by hand. `opts/` is an unused copy and can be deleted.
