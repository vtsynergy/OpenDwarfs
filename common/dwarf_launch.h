/* =============================================================================
 * dwarf_launch.h -- ONE pinned launch configuration shared by the OpenMP and OpenCL trees.
 *
 *   GPU, kernels whose size is NOT fixed by the algorithm ("free-size", barrier-free, one iteration/work-item):
 *       work-group (OpenCL local size) = team size (OpenMP thread_limit/num_threads) = DWARF_LS = 128
 *       work-groups (OpenCL global/DWARF_LS) = teams (OpenMP num_teams)             = ceil(N / 128)
 *       i.e. OpenCL NDRange {roundup(N,256)} local {256}  <=>  OpenMP num_teams(ceil(N/256)) thread_limit(256) num_threads(256)
 *   GPU, kernels whose work-group shape is part of the algorithm (lud 16 / 32 / 16x16, nw 16, fft L, kmeans 16x16, ...): the algorithmic
 *       size is used by BOTH trees, identically (documented in LAUNCH-CONFIGS.md).
 *   CPU (DWARF_USE_TARGET == 0): no team/thread clauses (they would create host "teams" of nested parallel regions); every CPU configuration
 *       uses the same 128 threads (OMP_NUM_THREADS=128, POCL_MAX_PTHREAD_COUNT=128, taskset -c 0-127).
 * Usage in OpenMP sources:   #pragma omp target teams distribute parallel for DWARF_LAUNCH(n) map(...)
 * Usage in OpenCL host code: local = DWARF_LS; global = DWARF_ROUNDUP(n);   and the kernel must guard `if (gid >= n) return;`
 * ============================================================================= */
#ifndef DWARF_LAUNCH_H
#define DWARF_LAUNCH_H
#include <stddef.h>
#ifndef DWARF_LS
#define DWARF_LS 128      /* free-size pin: 128 because NVHPC caps flat kernels at 128 threads/block (see analysis/ launch audit); Clang and OpenCL honour any value */
#endif
#ifndef DWARF_WG_ALG
#define DWARF_WG_ALG 256  /* algorithm-defined work-group size (nqueens WORK_ITEMS, 16x16 tiles of srad/hmm): identical in both trees, not subject to the free-size pin */
#endif
#define DWARF_NTEAMS(n) ((int)(((size_t)(n) + (size_t)DWARF_LS - 1) / (size_t)DWARF_LS))
#define DWARF_ROUNDUP(n) ((((size_t)(n) + (size_t)DWARF_LS - 1) / (size_t)DWARF_LS) * (size_t)DWARF_LS)
#if defined(DWARF_USE_TARGET) && DWARF_USE_TARGET
#define DWARF_LAUNCH(n) num_teams(DWARF_NTEAMS(n)) thread_limit(DWARF_LS) num_threads(DWARF_LS)
/* algorithmic work-group of t threads over n work-items */
#define DWARF_LAUNCH_T(n, t) num_teams((int)(((size_t)(n) + (size_t)(t) - 1) / (size_t)(t))) thread_limit(t) num_threads(t)
#else
#define DWARF_LAUNCH(n)
#define DWARF_LAUNCH_T(n, t)
#endif
#endif
