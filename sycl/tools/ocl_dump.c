/* LD_PRELOAD interposer for the OpenCL reference binaries: appends every device->host transfer to
 * $DWARF_DUMP_DIR/readback.bin as {uint64 bytes, data}, the stream sycl/include/dwarf_sycl.hpp writes.
 * Non-blocking reads are recorded when their queue is clFinish()ed. No OpenCL source is modified. */
#define _GNU_SOURCE
#define CL_TARGET_OPENCL_VERSION 120
#include <CL/cl.h>
#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
static struct { cl_command_queue q; const void *p; size_t n; } pend[1 << 16];
static int npend;
static void emit(const void *p, size_t bytes) {
  static FILE *f;
  const char *dir = getenv("DWARF_DUMP_DIR");
  if (!dir) return;
  if (!f) { char path[4096]; snprintf(path, sizeof path, "%s/readback.bin", dir); f = fopen(path, "wb"); if (!f) { perror(path); exit(1); } }
  unsigned long long n = bytes;
  fwrite(&n, sizeof n, 1, f); fwrite(p, 1, bytes, f); fflush(f);
}
cl_int clEnqueueReadBuffer(cl_command_queue q, cl_mem m, cl_bool blocking, size_t off, size_t n, void *p, cl_uint ne, const cl_event *w, cl_event *ev) {
  static __typeof__(&clEnqueueReadBuffer) real;
  if (!real) real = (__typeof__(real))dlsym(RTLD_NEXT, "clEnqueueReadBuffer");
  cl_int r = real(q, m, blocking, off, n, p, ne, w, ev);
  if (r == CL_SUCCESS) {
    if (blocking) emit(p, n);
    else { if (npend >= (1 << 16)) { fprintf(stderr, "ocl_dump: too many pending reads\n"); exit(1); } pend[npend].q = q; pend[npend].p = p; pend[npend].n = n; npend++; }
  }
  return r;
}
cl_int clFinish(cl_command_queue q) {
  static __typeof__(&clFinish) real;
  if (!real) real = (__typeof__(real))dlsym(RTLD_NEXT, "clFinish");
  cl_int r = real(q);
  int i, k = 0;
  for (i = 0; i < npend; i++) { if (pend[i].q == q) emit(pend[i].p, pend[i].n); else pend[k++] = pend[i]; }
  npend = k;
  return r;
}
