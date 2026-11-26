/* dwarf_cuda.h -- small host helpers shared by the CUDA translations (counterpart of opencl/include/common_args.*).
 *   CUDA_CHECK(call)          abort with file:line on any CUDA runtime error
 *   dwarf_cuda_init(&argc,&argv)  same CLI prefix as ocd_init(): "[-p N] [-d N] [-t N] [-o] --" (parsed only when "--" is
 *                             present; argc/argv are advanced exactly like ocd_parse() does). -d selects the CUDA device.
 *   DWARF_NBLOCKS(n)          blocks of DWARF_LS threads covering n work-items (OpenCL: DWARF_ROUNDUP(n)/DWARF_LS)
 *   CUDA_SYNC_LAUNCH()        after a launch that OpenCL follows with clFinish: reports launch errors and waits
 */
#ifndef DWARF_CUDA_H
#define DWARF_CUDA_H
#include <cuda_runtime.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dwarf_launch.h"

#define CUDA_CHECK(call) do { cudaError_t dc_e_ = (call); if (dc_e_ != cudaSuccess) { \
    fprintf(stdout, "CUDA Error %d (%s) at %s:%d: %s\n", (int)dc_e_, cudaGetErrorString(dc_e_), __FILE__, __LINE__, #call); \
    exit(1); } } while (0)
/* launch errors (bad block shape, too much shared memory, ...) are only visible through cudaGetLastError() */
#define CUDA_CHECK_LAUNCH() CUDA_CHECK(cudaGetLastError())
#define CUDA_SYNC_LAUNCH()  do { CUDA_CHECK_LAUNCH(); CUDA_CHECK(cudaDeviceSynchronize()); } while (0)
#define DWARF_NBLOCKS(n) ((unsigned int)(DWARF_ROUNDUP(n) / DWARF_LS))

static inline void dwarf_cuda_init(int *argc, char ***argv) {
  int n = *argc, i, dd = -1, dev = 0;
  char **v = *argv;
  for (i = 0; i < n; i++)
    if (strcmp(v[i], "--") == 0) { dd = i; break; }
  if (dd >= 0) {
    for (i = 1; i < dd; i++) {
      if (strcmp(v[i], "-d") == 0 && i + 1 < dd) dev = atoi(v[++i]);
      else if ((strcmp(v[i], "-p") == 0 || strcmp(v[i], "-t") == 0) && i + 1 < dd) ++i; /* OpenCL platform / device type */
    }
    *argv += dd;
    *argc -= dd;
  }
  if (dev < 0) dev = 0;
  CUDA_CHECK(cudaSetDevice(dev));
  cudaDeviceProp p;
  CUDA_CHECK(cudaGetDeviceProperties(&p, dev));
  printf("Device Chosen : %s (CUDA device %d, compute capability %d.%d)\n", p.name, dev, p.major, p.minor);
}
#endif
