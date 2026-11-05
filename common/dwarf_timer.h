/* =============================================================================
 * dwarf_timer.h -- common timing/reporting API for every backend tree
 * (opencl/, omp/, openacc/). Header-only, C99 and C++ compatible.
 *
 * TIMED REGION (identical definition in every tree):
 *   everything the device-offload implementation needs after the input is in
 *   host memory and before the result is in host memory: device-side allocation,
 *   host->device copies, all kernels, device->host copies, and the device
 *   synchronisation. EXCLUDED: file I/O, host-side initialisation, program/kernel
 *   JIT/compile, device/platform discovery, host verification.
 *   (One-time runtime start-up cost is therefore absorbed by the warm-up reps.)
 *
 * USAGE
 *     dwarf_init("lud", "dense-linear-algebra", "omp");   // once, after arg parsing
 *     dwarf_set_size(n, 0);                               // problem size (any units)
 *     for (int rep = dwarf_first_rep(); rep < dwarf_reps(); ++rep) {
 *         setup();                    // untimed (re-create inputs for destructive kernels)
 *         dwarf_region_begin();
 *         run();                      // the offload region; may call dwarf_phase_add(DWARF_PHASE_KERNEL, ms)
 *         dwarf_region_end(rep);      // rep < 0 are warm-ups: timed but not reported
 *     }
 * *     return dwarf_report(valid);     // prints DWARF_RESULT line, returns 0 / 2
 *
 * ENVIRONMENT (set by the run scripts)
 *     DWARF_WARMUP   warm-up reps inside one process (default 0)
 *     DWARF_REPEATS  measured reps inside one process (default 1)
 *     DWARF_TOOLCHAIN, DWARF_DEVICE   free-form labels copied into the result line
 *
 * OUTPUT (one line, stdout, parsed by scripts/):
 *   DWARF_RESULT,benchmark=..,motif=..,model=..,toolchain=..,device=..,size=..,
 *     secondary_size=..,repeats=..,mean_ms=..,median_ms=..,min_ms=..,max_ms=..,
 *     h2d_ms=..,kernel_ms=..,d2h_ms=..,valid=0|1
 * Phase fields are 0 when the backend cannot separate them (e.g. OpenMP).
 * ========================================================================== */
#ifndef DWARF_TIMER_H
#define DWARF_TIMER_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#ifdef __cplusplus
extern "C" {
#endif

#define DWARF_MAX_REPS 4096
enum { DWARF_PHASE_H2D = 0, DWARF_PHASE_KERNEL = 1, DWARF_PHASE_D2H = 2, DWARF_NPHASES = 3 };

typedef struct {
  const char *bench, *motif, *model;
  long long size, size2;
  int warmup, repeats;
  double t0, t1, acc, times[DWARF_MAX_REPS];
  int ntimes;
  double phase[DWARF_NPHASES]; /* summed over measured reps, divided at report */
  double pend[DWARF_NPHASES];  /* current rep, committed at region_end if measured */
} dwarf_state_t;

/* One instance shared by every translation unit of the program (weak definition, merged by the linker). */
#if defined(__GNUC__) || defined(__clang__) || defined(__NVCOMPILER)
__attribute__((weak)) dwarf_state_t dwarf__s;
#else
static dwarf_state_t dwarf__s;
#endif

static inline double dwarf_now_ms(void) {
  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  return ts.tv_sec * 1000.0 + ts.tv_nsec * 1e-6;
}
static inline int dwarf__env_int(const char *k, int d) {
  const char *v = getenv(k);
  return (v && *v) ? atoi(v) : d;
}
static inline void dwarf_init(const char *bench, const char *motif, const char *model) {
  memset(&dwarf__s, 0, sizeof dwarf__s);
  dwarf__s.bench = bench; dwarf__s.motif = motif; dwarf__s.model = model;
  dwarf__s.warmup = dwarf__env_int("DWARF_WARMUP", 0);
  dwarf__s.repeats = dwarf__env_int("DWARF_REPEATS", 1);
  if (dwarf__s.repeats < 1) dwarf__s.repeats = 1;
  if (dwarf__s.repeats > DWARF_MAX_REPS) dwarf__s.repeats = DWARF_MAX_REPS;
}
static inline void dwarf_set_size(long long a, long long b) { dwarf__s.size = a; dwarf__s.size2 = b; }
static inline void dwarf_set_repeats(int n) { if (n >= 1 && n <= DWARF_MAX_REPS) dwarf__s.repeats = n; }
static inline int dwarf_first_rep(void) { return -dwarf__s.warmup; }
static inline int dwarf_reps(void) { return dwarf__s.repeats; }
static inline void dwarf_region_begin(void) {
  memset(dwarf__s.pend, 0, sizeof dwarf__s.pend);
  dwarf__s.t0 = dwarf_now_ms();
}
static inline void dwarf_region_end(int rep) {
  double dt = dwarf_now_ms() - dwarf__s.t0;
  if (rep >= 0 && dwarf__s.ntimes < DWARF_MAX_REPS) {
    dwarf__s.times[dwarf__s.ntimes++] = dt;
    int p;
    for (p = 0; p < DWARF_NPHASES; ++p) dwarf__s.phase[p] += dwarf__s.pend[p];
  }
}
/* Add backend-measured time (ms) to a phase, between region_begin and region_end;
 * committed only for measured reps. */
static inline void dwarf_phase_add(int phase, double ms) {
  if (phase >= 0 && phase < DWARF_NPHASES) dwarf__s.pend[phase] += ms;
}
/* Accumulating variant for apps whose offload work is interleaved with host I/O (e.g. a database
 * search): bracket every offload span with dwarf_acc_begin()/dwarf_acc_end(); dwarf_acc_commit(rep)
 * closes the repetition with the accumulated time (same effect as region_begin/end around it). */
static inline void dwarf_acc_reset(void) { memset(dwarf__s.pend, 0, sizeof dwarf__s.pend); dwarf__s.t0 = 0; dwarf__s.acc = 0; }
static inline void dwarf_acc_begin(void) { dwarf__s.t1 = dwarf_now_ms(); }
static inline void dwarf_acc_end(void) { dwarf__s.acc += dwarf_now_ms() - dwarf__s.t1; }
static inline void dwarf_acc_commit(int rep) {
  if (rep >= 0 && dwarf__s.ntimes < DWARF_MAX_REPS) {
    int p;
    dwarf__s.times[dwarf__s.ntimes++] = dwarf__s.acc;
    for (p = 0; p < DWARF_NPHASES; ++p) dwarf__s.phase[p] += dwarf__s.pend[p];
  }
}
static inline int dwarf__cmp(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x > y) - (x < y);
}
static inline int dwarf_report(int valid) {
  int n = dwarf__s.ntimes;
  double sum = 0, mn = 0, mx = 0, med = 0, s[DWARF_MAX_REPS];
  if (n > 0) {
    memcpy(s, dwarf__s.times, n * sizeof(double));
    qsort(s, n, sizeof(double), dwarf__cmp);
    int i;
    for (i = 0; i < n; ++i) sum += s[i];
    mn = s[0]; mx = s[n - 1];
    med = (n % 2) ? s[n / 2] : 0.5 * (s[n / 2 - 1] + s[n / 2]);
  }
  const char *tc = getenv("DWARF_TOOLCHAIN"), *dev = getenv("DWARF_DEVICE");
  double div = n > 0 ? n : 1;
  printf("DWARF_RESULT,benchmark=%s,motif=%s,model=%s,toolchain=%s,device=%s,size=%lld,"
         "secondary_size=%lld,repeats=%d,mean_ms=%.6f,median_ms=%.6f,min_ms=%.6f,max_ms=%.6f,"
         "h2d_ms=%.6f,kernel_ms=%.6f,d2h_ms=%.6f,valid=%d\n",
         dwarf__s.bench, dwarf__s.motif, dwarf__s.model, tc ? tc : "unknown", dev ? dev : "unknown",
         dwarf__s.size, dwarf__s.size2, n, n ? sum / n : 0.0, med, mn, mx,
         dwarf__s.phase[0] / div, dwarf__s.phase[1] / div, dwarf__s.phase[2] / div, valid ? 1 : 0);
  fflush(stdout);
  return valid ? 0 : 2;
}

#ifdef __cplusplus
}
#endif
#endif /* DWARF_TIMER_H */
