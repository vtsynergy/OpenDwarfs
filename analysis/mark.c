/* libmark.so -- defines the timer's optional dwarf_probe_mark() hook (common/dwarf_timer.h).
 *   - logs every region/span event with the host monotonic clock to $MARK_OUT.<pid> ("B ns", "E ns rep", "S ns", "K ns rep");
 *   - MARK_PROFILER=1: calls cuProfilerStart at the first timed repetition and cuProfilerStop after the last one (for `nsys --capture-range=cudaProfilerApi`).
 * Build: gcc -shared -fPIC -O2 -o libmark.so mark.c -ldl.  Needs DWARF_WARMUP / DWARF_REPEATS in the environment (the same values given to the app). */
#define _GNU_SOURCE
#include <dlfcn.h>
#include <fcntl.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
static FILE *F; static int nbegin = 0, started = 0, stopped = 0;
static int geti(const char *k, int d) { const char *v = getenv(k); return v ? atoi(v) : d; }
static unsigned long now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return (unsigned long)t.tv_sec * 1000000000ul + (unsigned long)t.tv_nsec; }
static int ctlfd = -1, ackfd = -1;
static void perf_cmd(const char *c) {   /* perf stat --control fifo:CTL,ACK : enable/disable the counters around the timed repetitions */
    if (ctlfd == -1) { const char *p = getenv("MARK_PERF_CTL"), *a = getenv("MARK_PERF_ACK"); if (p && a) { ctlfd = open(p, O_WRONLY); ackfd = open(a, O_RDONLY); } else ctlfd = -2; }
    if (ctlfd >= 0 && ackfd >= 0) { char b[16]; if (write(ctlfd, c, strlen(c)) > 0) { ssize_t r = read(ackfd, b, sizeof b); (void)r; } }
}
static void prof(int start) {
    if (getenv("MARK_PERF_CTL")) perf_cmd(start ? "enable\n" : "disable\n");
    if (!getenv("MARK_PROFILER")) return;
    const char *dbg = getenv("MARK_DEBUG");
    int (*fn)(void) = dlsym(RTLD_DEFAULT, start ? "cuProfilerStart" : "cuProfilerStop");
    if (!fn) { void *h = dlopen("libcuda.so.1", RTLD_LAZY | RTLD_NOLOAD); if (h) fn = dlsym(h, start ? "cuProfilerStart" : "cuProfilerStop"); }
    if (!fn) { void *h = dlopen("libcudart.so", RTLD_LAZY | RTLD_NOLOAD); if (!h) h = dlopen(NULL, RTLD_LAZY); int (*f2)(void) = h ? dlsym(h, start ? "cudaProfilerStart" : "cudaProfilerStop") : NULL; if (f2) fn = f2; }
    if (dbg) fprintf(stderr, "[mark] %s -> %s\n", start ? "start" : "stop", fn ? "found" : "NOT FOUND");
    if (fn) fn();
}
void dwarf_probe_mark(int kind, int rep) {
    if (!F) { const char *b = getenv("MARK_OUT"); char p[1024]; snprintf(p, sizeof p, "%s.%d", b ? b : "mark", (int)getpid()); F = fopen(p, "w"); }
    int warm = geti("DWARF_WARMUP", 0), reps = geti("DWARF_REPEATS", 1), pr = geti("MARK_PROFILER", 0) || (getenv("MARK_PERF_CTL") != NULL);
    unsigned long t = now();
    if (kind == 1) { if (F) fprintf(F, "B %lu\n", t); if (pr && !started && nbegin == warm && getenv("MARK_ACC") == NULL) { prof(1); started = 1; } nbegin++; }
    else if (kind == 0) { if (F) fprintf(F, "E %lu %d\n", t, rep); if (pr && !stopped && rep == reps - 1) { prof(0); stopped = 1; } }
    else if (kind == 2) { if (F) fprintf(F, "S %lu\n", t); }
    else if (kind == 3) { if (F) fprintf(F, "K %lu %d\n", t, rep); if (pr && !started && rep == -1) { prof(1); started = 1; } if (pr && !stopped && rep == reps - 1) { prof(0); stopped = 1; } }
    if (F) fflush(F);
}
