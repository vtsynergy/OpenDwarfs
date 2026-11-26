/* dwarf_acc.h -- launch shapes of common/dwarf_launch.h expressed as OpenACC clauses.
 * GPU builds (-DDWARF_USE_TARGET=1, set by ../dwarf.mk for ACC_TARGET=gpu):
 *   free-size kernel over n work-items:  gang vector, vector_length(DWARF_LS=128), num_gangs(ceil(n/128))
 *     <=> OpenCL global = roundup(n,128), local = 128
 *   algorithm-defined work-group of t over n work-items: vector_length(t), num_gangs(ceil(n/t))
 * Multicore builds: the clauses expand to nothing (num_gangs would set the number of host threads);
 *   the thread count is ACC_NUM_CORES / the CPU affinity mask.
 * Usage:  #pragma acc parallel loop gang vector ACC_LAUNCH(n) ...
 * Whether the compiler/runtime honoured the shape is checked with -Minfo=accel and NV_ACC_NOTIFY=2. */
#ifndef DWARF_ACC_H
#define DWARF_ACC_H
#include "dwarf_launch.h"
#if defined(DWARF_USE_TARGET) && DWARF_USE_TARGET
#define ACC_LAUNCH(n) num_gangs(DWARF_NTEAMS(n)) vector_length(DWARF_LS)
#define ACC_LAUNCH_T(n, t) num_gangs((int)(((size_t)(n) + (size_t)(t) - 1) / (size_t)(t))) vector_length(t)
#define ACC_SHAPE(g, v) num_gangs(g) vector_length(v)
#else
#define ACC_LAUNCH(n)
#define ACC_LAUNCH_T(n, t)
#define ACC_SHAPE(g, v)
#endif
#endif
