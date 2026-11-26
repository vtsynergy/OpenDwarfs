# CUDA translation notes

## Mapping
Block shapes equal the OpenCL work-groups: free-size kernels use `DWARF_LS`=128 threads and `ceil(N/128)` blocks with the
same index guard; algorithmic shapes are unchanged (lud 16, 32, 16x16; nw 16; srad 16x16; hmm 16x16; nqueens 256; swat
14 blocks of 64; gem 120 blocks of 256; fft as computed by the OpenCL host). `clFinish` -> `cudaDeviceSynchronize`,
blocking read/write -> `cudaMemcpy` (+ synchronise where OpenCL calls `clFinish`), `__local` -> `__shared__`,
`barrier` -> `__syncthreads`, global fences -> `__threadfence`, `atom_inc/atom_add` -> `atomicAdd`.
The behaviour fixes of `../opencl/FIXES.md` that concern results or validators are carried over. The HMM `estimate_b`
no-op listed there as "not fixed" is reproduced as is.

## CUDA-WA-01: nvcc miscompiles `lud_perimeter`
nvcc 12.9 (V12.9.86), `-O3`, `sm_90`: in the unrolled copy of the peri-row inner loop the stores apply the +16 column
offset twice. In the last block they wrap into the diagonal block; the LUD validator failed (about 6-7 thousand wrong
elements at n=256, over half a million at n=2048, varying between runs). The PTX shows it. `-Xptxas -O0` still fails, `-Xcicc -O0` passes, so the fault is in the
code produced before ptxas. Workaround: `#pragma unroll 1` on that loop only (`lud_kernel.cuh`).
Reproducer: `repro/cuda_wa01.cu`; command: `GPU_INDEX=1 repro/run.sh`; recorded output and PTX reading:
`repro/EVIDENCE.txt` (re-run 2026-11-25: still reproduces). Not reported to NVIDIA; other nvcc versions not tried.

## Floating point
NVIDIA OpenCL's default single-precision divide and sqrt are **not** correctly rounded (OpenCL allows 2.5 and 3 ulp unless
the program is built with `-cl-fp32-correctly-rounded-divide-sqrt`; no application here passes it; `fft` passes only
`-cl-mad-enable`). nvcc defaults to correctly rounded. The tree therefore builds with `-prec-div=false -prec-sqrt=false`.
Evidence: with these flags the LUD factors are bitwise identical to OpenCL (256, 2048); with nvcc defaults both pass the
`L*U` check but the factors differ widely (relative L2 0.94 and 1.35; elimination without pivoting amplifies rounding).
Anyone comparing CUDA with OpenCL or with another model must treat this flag as part of the configuration. `fmad` is left at
the nvcc default (on); on `cfd` this was the closest to OpenCL of the four settings tried.

## Host compiler flags
gcc/g++ 13.3.0 (Ubuntu), `-O3`, no `-march`/`-mtune`: baseline `x86-64`, generic tuning, AVX2 and FMA disabled. No
`CFLAGS`, `CXXFLAGS` or `NVCC_*_FLAGS` were set in the environment. `.c`: `-std=gnu89 -fcommon`; `.cpp`:
`-std=gnu++14 -fpermissive`; host part of `.cu` (from `nvcc --dryrun`): `-std=c++14 -fpermissive -O3 -m64`. The OpenCL
tree (`../opencl/common.mk`) also uses `-O3` without `-march`.
`kmeans` `rms_err` is serial, as in OpenCL: `rmse.c` is built without `-fopenmp` and the binary imports no OpenMP symbol.

## Timed-region audit (OpenCL vs CUDA)
Line numbers: OpenCL `../opencl/...` -> CUDA. "out" = before/after the region, "in" = between region (or accumulate-span)
begin and end, including helpers called from there. **No mismatch was found; no source change was needed.**

| App (file) | Region lines | Device allocations | Host-to-device copies | Same? |
|---|---|---|---|---|
| lud (`lud.c` -> `lud.cu`) | 144-222 -> 104-129 | matrix out (132 -> 93) | matrix in each rep (146 -> 106) | yes |
| spmv (`csr.c` -> `csr.cu`) | 385-431 -> 312-334 | 5 arrays out (324-328 -> 272-276) | 5 in (392-404 -> 319-323) | yes (1) |
| bfs (`bfs.cpp` -> `bfs.cu`) | 302-350 -> 182-213 | 7 out (166-212 -> 136-159) | initial untimed 6 out (180-213 -> 142-160); 6 in (303-308 -> 183-188) + stop flag per level (313 -> 193) | yes |
| crc (`crc_algo.c` -> `.cu`) | 367-393 -> 317-341 | per-block in/out buffers out (343,345 -> 295,296); table is static | pages in, via `enqueueCRCDevice` (116 -> 106) | yes (1) |
| nw (`needle.c` -> `.cu`) | 199-270 -> 164-195 | 2 out (187,189 -> 154,155) | 2 in (200,206 -> 165,167) | yes |
| srad (`srad.c` -> `.cu`) | 282-419 -> 194-296 | 6 out (249-259 -> 172) | J per iteration in (364 -> 275) | yes |
| cfd (`cfd.cpp` -> `.cu`) | 449-553 -> 322-359 | all out (257-397 -> 195-306) | far field and first mesh upload out; areas, neighbours, normals, step_factors in (450-452,465 -> 323-325,333) | yes |
| gem (`calculate_potential.cpp` -> `.cu`) | 531-675 -> 302-327 | 17 out (358-374 -> 264-280) | first upload of 17 out (377-441 -> 282); all 17 again in (532-596 -> 303) | yes |
| tdm (`tdm_ocl.cpp` -> `tdm_cuda.cu`) | 314-432 -> 274-344 | 9 out (192-207 -> 188-203) | events, times, candidates, intervals once, out (209-227 -> 206-209); start records and offsets in (319,382 -> 279,314) | yes |
| swat (`swat.cpp` -> `.cu`) | 401-579 -> 299-367 (accumulated spans) | all out (193-311 -> 110-210) | score matrix out (302 -> 207); sequences, diffPos, threadNum in (470-484 -> 328-331) | yes |
| hmm (`main_bwa_hmm.c` -> `.cu`) | 987-1127 -> 583-685 | 10 model buffers in (1007-1017 -> 604-613); `dot_production` allocates per call, in (313 -> 177) | a, b, pi in (1031-1043 -> 616-622); per-call zero fills in (317,655,697,831 -> 180,360,388,478) | yes (2) |
| nqueens (`nqueen_cl.cpp` -> `nqueen_cuda.cu`; region in `main`) | around `Compute()` | 4 in (354-367 -> 98-107) | all in, inside `Compute()` | yes (2) |
| kmeans (`kmeans.c`; work in `kmeans_opencl.cpp` -> `kmeans_cuda.cu`) | 260-272 -> same | all in (`allocateMemory` via `cluster()`) | all in | yes (2) |
| fft (`fft.cpp` -> `fft_cuda.cu`) | 243-253 -> 234-242 | 2 pinned host buffers out; work and temp in (245,246 -> 235,236) | in (247 -> 237) | yes (2) |

(1) OpenCL also creates command queues (and, in spmv, a kernel object) each repetition before the region; CUDA has no
counterpart. (2) Same source position, but not necessarily the same runtime work: `cudaMalloc` allocates at once, whereas
OpenCL implementations may delay the real allocation of a `clCreateBuffer` to first use. This was not measured. It matters
for apps that allocate inside the region (hmm, nqueens, kmeans, fft), and for those that allocate outside it the first
OpenCL repetition may include work that CUDA did earlier.

## Deliberate differences
* `FPFLAGS` and CUDA-WA-01 (above).
* lud, spmv, crc, nw, srad, hmm host files were C and are now compiled as C++ (`.cu`). Where C evaluates in double and C++
  would pick a float overload, the cast is explicit (`srad` `exp`; `hmm` `log10`, `pow`).
* spmv, crc: the three OpenCL queues become the default stream with `cudaMemcpyAsync` (pageable memory).
* OpenCL `__constant` buffers (kmeans centres, nqueens forbidden masks) are ordinary device memory.
* cfd uses CUDA `float3` (12 bytes, the OpenCL host's layout); only x, y, z are used.
* `native_cos/sin` -> `__cosf/__sinf` (fft); `native_exp` -> `__expf` (gem; not reached with the benchmark arguments).
* Only code the OpenCL build compiles for this GPU is translated: no FPGA or CPU variants; nqueens has the atomics,
  non-vectorised, 256-thread variant only (other `-blocksize` values are rejected; `-local`, `-novec`, `-noatomics`,
  `-vec4` are ignored); kmeans GPU-side reductions are off, as in OpenCL.
* fft: the SHOC host framework is replaced by a small `main` accepting `--pts` and `--passes`; the 2-D path is not
  translated.
* swat: `CL_DEVICE_MAX_COMPUTE_UNITS` -> `multiProcessorCount` (14 blocks are used in both).
* Kernels are compiled ahead of time; kernel-file options such as `-k` have no effect.
* `make clean` deletes nothing.

## SCALE (Spectral Compute): researched, NOT downloaded, NOT installed
* Distribution: apt/dnf repositories (need sudo), Docker image, and for users without root a tarball, so an install into
  `~/.local/apps/scale-hudson` is supported by unpacking it. Latest: 1.7.5 (2026-11-03),
  `https://pkgs.scale-lang.com/tar/scale-1.7.5-amd64.tar.xz`, about 3.74 GB; a SHA-512 file is published next to it
  (`.tar.xz.sha512`). No signature file was found (`.asc` of the "latest" tarball gives 404).
* License: "SCALE Free License". Clause 1.3: by downloading, installing and using the software you agree to all terms, so
  even the download accepts it; this was not done on the user's behalf. Clause 4.1: non-commercial purposes only (no payment
  or benefit from a third party; no use in the course of a business or for-profit undertaking); 4.1/6.2: software and
  compiled outputs may be passed on only free of charge for non-commercial use; 4.2: no reverse engineering or modification;
  3.1: revocable; 3.2: feedback is licensed to Spectral; 11.1: terms may change; 12: law of England and Wales.
  The FAQ says free for non-commercial use including research and academia. The pricing page says "free for evaluation and
  non-commercial use"; the word "evaluation" is not in the free license text. Commercial use is paid (listed from
  USD 150 per GPU per year for consumer/workstation, 800 and 1500 for datacenter classes).
* H100: documented. SCALE has an NVIDIA mode in which it replaces the nvcc compiler only and uses NVIDIA's libraries; the
  target table lists compute capability 9.0 (`sm_90`, `sm_90a`; the example GPU shown is the H200). Not verified on this
  host. Their compiler is Clang/LLVM based and does not use NVIDIA's nvcc or nvvm.
* OS: Ubuntu 22.04/24.04, Rocky/RHEL, other distributions via tarball. This host is Ubuntu 24.04.5, glibc 2.39.
* To decide: whether the intended use is non-commercial in the sense of clause 4.1 and whether to accept the license.
  If yes (the user runs this, in a new empty directory):

      curl -LO https://pkgs.scale-lang.com/tar/scale-1.7.5-amd64.tar.xz
      curl -LO https://pkgs.scale-lang.com/tar/scale-1.7.5-amd64.tar.xz.sha512
      sha512sum -c scale-1.7.5-amd64.tar.xz.sha512
      mkdir -p ~/.local/apps/scale-hudson && tar xf scale-1.7.5-amd64.tar.xz -C ~/.local/apps/scale-hudson

  The archive layout and the checksum file format were not inspected. Then `source <root>/bin/scaleenv sm_90`, check
  `nvcc --version`, and run `CUDACC=nvcc BUILD_DIR=build-scale ./run_all.sh`.
  `toolchains/install-cuda-scale.sh` was not written because nothing was installed.
* URLs: https://docs.scale-lang.com/stable/licensing/ ,
  https://docs.scale-lang.com/stable/manual/tutorials/how-to-install/ ,
  https://docs.scale-lang.com/stable/manual/tutorials/how-to-use/ , https://docs.scale-lang.com/stable/ ,
  https://docs.scale-lang.com/stable/manual/faq/ , https://scale-lang.com/pricing , https://pkgs.scale-lang.com/tar/
