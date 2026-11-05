/* Host-side helpers that the OpenCL tree gets from include/common_args.h, without any
 * OpenCL dependency. Lets host code (input readers, CPU reference code) be copied verbatim
 * into the OpenMP translations; only the #include line is rewritten to this header. */
#ifndef OCD_COMPAT_H
#define OCD_COMPAT_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#ifdef __cplusplus
extern "C" {
#endif
#define MINIMUM(i, j) ((i) < (j) ? (i) : (j))
#define ACL_ALIGNMENT 64
static int _deviceType = -1; /* referenced by some host code; the OpenMP build has no FPGA path */
static inline void check(int b, const char *msg) {
  if (!b) { fprintf(stderr, "error: %s\n\n", msg); exit(-1); }
}
static inline void *char_new_array(const size_t N, const char *m) { void *p = malloc(N * sizeof(char)); check(p != NULL, m); return p; }
static inline void *int_new_array(const size_t N, const char *m) { void *p = malloc(N * sizeof(int)); check(p != NULL, m); return p; }
static inline void *long_new_array(const size_t N, const char *m) { void *p = malloc(N * sizeof(long)); check(p != NULL, m); return p; }
static inline void *float_new_array(const size_t N, const char *m) { void *p = malloc(N * sizeof(float)); check(p != NULL, m); return p; }
static inline void *float_array_realloc(void *ptr, const size_t N, const char *m) {
  void *p = realloc(ptr, N * sizeof(float)); check(p != NULL, m); return p;
}
#ifdef __cplusplus
}
#endif
#endif
