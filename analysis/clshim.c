/* =============================================================================
 * clshim.c -- LD_PRELOAD interposer for the OpenCL API (any ICD: PoCL, Intel CPU runtime, NVIDIA).
 *
 * Records, without changing what the application computes:
 *   - every clEnqueueNDRangeKernel: kernel name, work_dim, global/local sizes (local = 0 when NULL), launch count;
 *   - device-side timing of every launch from OpenCL event profiling (the shim switches CL_QUEUE_PROFILING_ENABLE on and attaches
 *     its own event when the application passes none): kernel execution time and queued->start latency (dispatch latency);
 *   - host time spent inside each API class (enqueue, buffer write/read/copy/map/fill, finish/wait/flush, buffer create/release,
 *     clSetKernelArg, program build / kernel creation = "init"), with call counts and bytes moved.
 * Output: one JSON file  $CLSHIM_OUT.<pid>  written at process exit, only by processes that called the OpenCL API.
 *
 * Build:  gcc -shared -fPIC -O2 -o libclshim.so clshim.c -ldl -lpthread            Use:  LD_PRELOAD=./libclshim.so CLSHIM_OUT=/path/run ./app
 * Limits: clEnqueueReadBufferRect/WriteBufferRect/CopyBufferRect, images and SVM are not interposed (their time shows up as unaccounted host time).
 * ============================================================================= */
#define _GNU_SOURCE
#define CL_TARGET_OPENCL_VERSION 200
#define CL_USE_DEPRECATED_OPENCL_1_2_APIS
#define CL_USE_DEPRECATED_OPENCL_2_0_APIS
#include <CL/cl.h>
#include <dlfcn.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

enum { K_ENQ, K_WRITE, K_READ, K_COPY, K_MAP, K_UNMAP, K_FILL, K_FINISH, K_WAIT, K_FLUSH, K_CREATE_BUF, K_RELEASE_MEM, K_SETARG, K_INIT, K_NCLASS };
static const char *kname[K_NCLASS] = {"enqueue_ndrange", "write_buffer", "read_buffer", "copy_buffer", "map_buffer", "unmap", "fill_buffer", "finish", "wait_events", "flush",
                                      "create_buffer", "release_mem", "set_kernel_arg", "init"};
static struct { unsigned long count, ns, bytes; } cls[K_NCLASS];
static pthread_mutex_t mu = PTHREAD_MUTEX_INITIALIZER;
static int used = 0;

typedef struct { char name[112]; cl_uint dim; size_t g[3], l[3]; unsigned long count; double dev_ns, q2s_ns, submit2s_ns, idle_ns; } krec_t;
#define MAXK 512
static krec_t K[MAXK]; static int nK = 0;
typedef struct { cl_event ev; int rec; cl_command_queue q; uint64_t t0, t1; } ev_t;
typedef struct { unsigned char c; uint64_t s, e; } call_t;
static call_t *CALL = NULL; static size_t nCALL = 0, capCALL = 0;
static void push_call(int c, uint64_t s, uint64_t e) { if (nCALL == capCALL) { capCALL = capCALL ? capCALL * 2 : 65536; CALL = realloc(CALL, capCALL * sizeof *CALL); } CALL[nCALL].c = (unsigned char)c; CALL[nCALL].s = s; CALL[nCALL].e = e; nCALL++; }
static ev_t *EV = NULL; static size_t nEV = 0, capEV = 0; static unsigned long ev_fail = 0;
#define MAXKH 256
static struct { cl_kernel k; char name[112]; } KH[MAXKH]; static int nKH = 0;

static inline uint64_t now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return (uint64_t)t.tv_sec * 1000000000ull + (uint64_t)t.tv_nsec; }
#define REAL(ret, name, args) static ret (*real_##name) args; if (!real_##name) real_##name = dlsym(RTLD_NEXT, #name)
static void add(int c, uint64_t ns, unsigned long bytes) { uint64_t e = now(); pthread_mutex_lock(&mu); used = 1; cls[c].count++; cls[c].ns += ns; cls[c].bytes += bytes; push_call(c, e - ns, e); pthread_mutex_unlock(&mu); }

static const char *kernel_name(cl_kernel k) {
    REAL(cl_int, clGetKernelInfo, (cl_kernel, cl_kernel_info, size_t, void *, size_t *));
    for (int i = 0; i < nKH; i++) if (KH[i].k == k) return KH[i].name;
    char n[112] = "?"; if (real_clGetKernelInfo) real_clGetKernelInfo(k, CL_KERNEL_FUNCTION_NAME, sizeof n, n, NULL);
    pthread_mutex_lock(&mu); if (nKH < MAXKH) { KH[nKH].k = k; strncpy(KH[nKH].name, n, 111); KH[nKH].name[111] = 0; nKH++; } pthread_mutex_unlock(&mu);
    static char tmp[112]; strncpy(tmp, n, 111); tmp[111] = 0; return tmp;
}
static int krec(const char *name, cl_uint dim, const size_t *g, const size_t *l) {
    size_t gg[3] = {1, 1, 1}, ll[3] = {0, 0, 0};
    for (cl_uint i = 0; i < dim && i < 3; i++) { gg[i] = g ? g[i] : 0; ll[i] = l ? l[i] : 0; }
    for (int i = 0; i < nK; i++) if (K[i].dim == dim && !strcmp(K[i].name, name) && !memcmp(K[i].g, gg, sizeof gg) && !memcmp(K[i].l, ll, sizeof ll)) return i;
    if (nK >= MAXK) return MAXK - 1;
    krec_t *r = &K[nK]; memset(r, 0, sizeof *r); memcpy(r->name, name, strnlen(name, 110)); r->dim = dim; memcpy(r->g, gg, sizeof gg); memcpy(r->l, ll, sizeof ll); return nK++;
}

/* ---- queue creation: force profiling on ---- */
cl_command_queue clCreateCommandQueue(cl_context c, cl_device_id d, cl_command_queue_properties p, cl_int *e) {
    REAL(cl_command_queue, clCreateCommandQueue, (cl_context, cl_device_id, cl_command_queue_properties, cl_int *));
    uint64_t t0 = now(); cl_command_queue q = real_clCreateCommandQueue(c, d, p | CL_QUEUE_PROFILING_ENABLE, e); add(K_INIT, now() - t0, 0); return q;
}
cl_command_queue clCreateCommandQueueWithProperties(cl_context c, cl_device_id d, const cl_queue_properties *props, cl_int *e) {
    REAL(cl_command_queue, clCreateCommandQueueWithProperties, (cl_context, cl_device_id, const cl_queue_properties *, cl_int *));
    cl_queue_properties np[16]; int n = 0, have = 0;
    for (const cl_queue_properties *p = props; p && *p && n < 12; p += 2) {
        np[n++] = p[0]; np[n++] = p[1];
        if (p[0] == CL_QUEUE_PROPERTIES) { np[n - 1] |= CL_QUEUE_PROFILING_ENABLE; have = 1; }
    }
    if (!have) { np[n++] = CL_QUEUE_PROPERTIES; np[n++] = CL_QUEUE_PROFILING_ENABLE; }
    np[n++] = 0; uint64_t t0 = now(); cl_command_queue q = real_clCreateCommandQueueWithProperties(c, d, np, e); add(K_INIT, now() - t0, 0); return q;
}

/* ---- launches ---- */
cl_int clEnqueueNDRangeKernel(cl_command_queue q, cl_kernel k, cl_uint dim, const size_t *off, const size_t *g, const size_t *l, cl_uint nw, const cl_event *wl, cl_event *ev) {
    REAL(cl_int, clEnqueueNDRangeKernel, (cl_command_queue, cl_kernel, cl_uint, const size_t *, const size_t *, const size_t *, cl_uint, const cl_event *, cl_event *));
    REAL(cl_int, clRetainEvent, (cl_event));
    cl_event mine = NULL; cl_event *evp = ev ? ev : &mine;
    uint64_t t0 = now(); cl_int r = real_clEnqueueNDRangeKernel(q, k, dim, off, g, l, nw, wl, evp); uint64_t t1 = now();
    const char *nm = kernel_name(k);
    pthread_mutex_lock(&mu); used = 1; cls[K_ENQ].count++; cls[K_ENQ].ns += t1 - t0; push_call(K_ENQ, t0, t1);
    int ri = krec(nm, dim, g, l); K[ri].count++;
    if (r == CL_SUCCESS && *evp) {
        if (ev && real_clRetainEvent) real_clRetainEvent(*ev);
        if (nEV == capEV) { capEV = capEV ? capEV * 2 : 4096; EV = realloc(EV, capEV * sizeof *EV); }
        EV[nEV].ev = *evp; EV[nEV].rec = ri; EV[nEV].q = q; EV[nEV].t0 = t0; EV[nEV].t1 = t1; nEV++;
    }
    pthread_mutex_unlock(&mu); return r;
}

/* ---- transfers / sync / memory ---- */
cl_int clEnqueueWriteBuffer(cl_command_queue q, cl_mem m, cl_bool b, size_t o, size_t n, const void *p, cl_uint nw, const cl_event *wl, cl_event *ev) {
    REAL(cl_int, clEnqueueWriteBuffer, (cl_command_queue, cl_mem, cl_bool, size_t, size_t, const void *, cl_uint, const cl_event *, cl_event *));
    uint64_t t0 = now(); cl_int r = real_clEnqueueWriteBuffer(q, m, b, o, n, p, nw, wl, ev); add(K_WRITE, now() - t0, n); return r; }
cl_int clEnqueueReadBuffer(cl_command_queue q, cl_mem m, cl_bool b, size_t o, size_t n, void *p, cl_uint nw, const cl_event *wl, cl_event *ev) {
    REAL(cl_int, clEnqueueReadBuffer, (cl_command_queue, cl_mem, cl_bool, size_t, size_t, void *, cl_uint, const cl_event *, cl_event *));
    uint64_t t0 = now(); cl_int r = real_clEnqueueReadBuffer(q, m, b, o, n, p, nw, wl, ev); add(K_READ, now() - t0, n); return r; }
cl_int clEnqueueCopyBuffer(cl_command_queue q, cl_mem s, cl_mem d, size_t so, size_t dof, size_t n, cl_uint nw, const cl_event *wl, cl_event *ev) {
    REAL(cl_int, clEnqueueCopyBuffer, (cl_command_queue, cl_mem, cl_mem, size_t, size_t, size_t, cl_uint, const cl_event *, cl_event *));
    uint64_t t0 = now(); cl_int r = real_clEnqueueCopyBuffer(q, s, d, so, dof, n, nw, wl, ev); add(K_COPY, now() - t0, n); return r; }
cl_int clEnqueueFillBuffer(cl_command_queue q, cl_mem m, const void *pat, size_t ps, size_t o, size_t n, cl_uint nw, const cl_event *wl, cl_event *ev) {
    REAL(cl_int, clEnqueueFillBuffer, (cl_command_queue, cl_mem, const void *, size_t, size_t, size_t, cl_uint, const cl_event *, cl_event *));
    uint64_t t0 = now(); cl_int r = real_clEnqueueFillBuffer(q, m, pat, ps, o, n, nw, wl, ev); add(K_FILL, now() - t0, n); return r; }
void *clEnqueueMapBuffer(cl_command_queue q, cl_mem m, cl_bool b, cl_map_flags f, size_t o, size_t n, cl_uint nw, const cl_event *wl, cl_event *ev, cl_int *e) {
    REAL(void *, clEnqueueMapBuffer, (cl_command_queue, cl_mem, cl_bool, cl_map_flags, size_t, size_t, cl_uint, const cl_event *, cl_event *, cl_int *));
    uint64_t t0 = now(); void *r = real_clEnqueueMapBuffer(q, m, b, f, o, n, nw, wl, ev, e); add(K_MAP, now() - t0, n); return r; }
cl_int clEnqueueUnmapMemObject(cl_command_queue q, cl_mem m, void *p, cl_uint nw, const cl_event *wl, cl_event *ev) {
    REAL(cl_int, clEnqueueUnmapMemObject, (cl_command_queue, cl_mem, void *, cl_uint, const cl_event *, cl_event *));
    uint64_t t0 = now(); cl_int r = real_clEnqueueUnmapMemObject(q, m, p, nw, wl, ev); add(K_UNMAP, now() - t0, 0); return r; }
cl_int clFinish(cl_command_queue q) { REAL(cl_int, clFinish, (cl_command_queue)); uint64_t t0 = now(); cl_int r = real_clFinish(q); add(K_FINISH, now() - t0, 0); return r; }
cl_int clFlush(cl_command_queue q) { REAL(cl_int, clFlush, (cl_command_queue)); uint64_t t0 = now(); cl_int r = real_clFlush(q); add(K_FLUSH, now() - t0, 0); return r; }
cl_int clWaitForEvents(cl_uint n, const cl_event *l) { REAL(cl_int, clWaitForEvents, (cl_uint, const cl_event *)); uint64_t t0 = now(); cl_int r = real_clWaitForEvents(n, l); add(K_WAIT, now() - t0, 0); return r; }
cl_mem clCreateBuffer(cl_context c, cl_mem_flags f, size_t n, void *h, cl_int *e) {
    REAL(cl_mem, clCreateBuffer, (cl_context, cl_mem_flags, size_t, void *, cl_int *));
    uint64_t t0 = now(); cl_mem r = real_clCreateBuffer(c, f, n, h, e); add(K_CREATE_BUF, now() - t0, n); return r; }
cl_int clReleaseMemObject(cl_mem m) { REAL(cl_int, clReleaseMemObject, (cl_mem)); uint64_t t0 = now(); cl_int r = real_clReleaseMemObject(m); add(K_RELEASE_MEM, now() - t0, 0); return r; }
cl_int clSetKernelArg(cl_kernel k, cl_uint i, size_t s, const void *v) { REAL(cl_int, clSetKernelArg, (cl_kernel, cl_uint, size_t, const void *)); uint64_t t0 = now(); cl_int r = real_clSetKernelArg(k, i, s, v); add(K_SETARG, now() - t0, 0); return r; }
/* ---- init class (excluded from per-repetition budgets by the collector) ---- */
cl_int clBuildProgram(cl_program p, cl_uint n, const cl_device_id *d, const char *o, void (CL_CALLBACK *cb)(cl_program, void *), void *u) {
    REAL(cl_int, clBuildProgram, (cl_program, cl_uint, const cl_device_id *, const char *, void (CL_CALLBACK *)(cl_program, void *), void *));
    uint64_t t0 = now(); cl_int r = real_clBuildProgram(p, n, d, o, cb, u); add(K_INIT, now() - t0, 0); return r; }
cl_kernel clCreateKernel(cl_program p, const char *nm, cl_int *e) { REAL(cl_kernel, clCreateKernel, (cl_program, const char *, cl_int *)); uint64_t t0 = now(); cl_kernel r = real_clCreateKernel(p, nm, e); add(K_INIT, now() - t0, 0); return r; }

static int cmp_i64_impl(const void *a, const void *b) { int64_t x = *(const int64_t *)a, y = *(const int64_t *)b; return x < y ? -1 : x > y; }
int cmp_i64(const void *a, const void *b) { return cmp_i64_impl(a, b); }
/* ---- report ---- */
__attribute__((destructor)) static void report(void) {
    if (!used) return;
    REAL(cl_int, clGetEventProfilingInfo, (cl_event, cl_profiling_info, size_t, void *, size_t *));
    REAL(cl_int, clReleaseEvent, (cl_event));
    cl_command_queue lq[16]; cl_ulong lend[16]; int nlq = 0;
    int64_t *offs = malloc((nEV + 1) * sizeof *offs); cl_ulong *ES = malloc((nEV + 1) * sizeof *ES), *EE = malloc((nEV + 1) * sizeof *EE); size_t noff = 0; char *okv = calloc(nEV + 1, 1);
    for (size_t i = 0; i < nEV; i++) {
        cl_ulong q = 0, sub = 0, s = 0, e = 0;
        if (real_clGetEventProfilingInfo && real_clGetEventProfilingInfo(EV[i].ev, CL_PROFILING_COMMAND_QUEUED, 8, &q, NULL) == CL_SUCCESS &&
            real_clGetEventProfilingInfo(EV[i].ev, CL_PROFILING_COMMAND_SUBMIT, 8, &sub, NULL) == CL_SUCCESS &&
            real_clGetEventProfilingInfo(EV[i].ev, CL_PROFILING_COMMAND_START, 8, &s, NULL) == CL_SUCCESS &&
            real_clGetEventProfilingInfo(EV[i].ev, CL_PROFILING_COMMAND_END, 8, &e, NULL) == CL_SUCCESS) {
            int sl = -1; for (int j = 0; j < nlq; j++) if (lq[j] == EV[i].q) sl = j;
            if (sl < 0 && nlq < 16) { sl = nlq++; lq[sl] = EV[i].q; lend[sl] = 0; }
            if (sl >= 0) { if (lend[sl] && s > lend[sl]) K[EV[i].rec].idle_ns += (double)(s - lend[sl]); if (e > lend[sl]) lend[sl] = e; }
            ES[i] = s; EE[i] = e; okv[i] = 1; offs[noff++] = (int64_t)((EV[i].t0 + EV[i].t1) / 2) - (int64_t)q;
            K[EV[i].rec].dev_ns += (double)(e - s); K[EV[i].rec].q2s_ns += (double)(s - q); K[EV[i].rec].submit2s_ns += (double)(s - sub);
        } else ev_fail++;
        if (real_clReleaseEvent) real_clReleaseEvent(EV[i].ev);
    }
    const char *base = getenv("CLSHIM_OUT"); char path[1024];
    snprintf(path, sizeof path, "%s.%d", base ? base : "clshim", (int)getpid());
    FILE *f = fopen(path, "w"); if (!f) return;
    fprintf(f, "{\"api\":{");
    for (int c = 0; c < K_NCLASS; c++) fprintf(f, "%s\"%s\":{\"count\":%lu,\"total_ns\":%lu,\"bytes\":%lu}", c ? "," : "", kname[c], cls[c].count, cls[c].ns, cls[c].bytes);
    fprintf(f, "},\"events_without_profile\":%lu,\"kernels\":[", ev_fail);
    for (int i = 0; i < nK; i++)
        fprintf(f, "%s{\"name\":\"%s\",\"dim\":%u,\"global\":[%zu,%zu,%zu],\"local\":[%zu,%zu,%zu],\"count\":%lu,\"dev_ns\":%.0f,\"queued_to_start_ns\":%.0f,\"submit_to_start_ns\":%.0f,\"idle_before_ns\":%.0f}",
                i ? "," : "", K[i].name, K[i].dim, K[i].g[0], K[i].g[1], K[i].g[2], K[i].l[0], K[i].l[1], K[i].l[2], K[i].count, K[i].dev_ns, K[i].q2s_ns, K[i].submit2s_ns, K[i].idle_ns);
    fprintf(f, "]}\n"); fclose(f);
    /* timeline file: host-call intervals and event intervals on the host monotonic clock (device clock aligned via median(enqueue time - QUEUED stamp)) */
    int cmp_i64(const void *a, const void *b);
    qsort(offs, noff, sizeof *offs, cmp_i64); int64_t off = noff ? offs[noff / 2] : 0;
    char p2[1100]; snprintf(p2, sizeof p2, "%s.tl", path); FILE *g = fopen(p2, "w");
    if (g) {
        fprintf(g, "# classes:"); for (int c = 0; c < K_NCLASS; c++) fprintf(g, " %d=%s", c, kname[c]); fprintf(g, "\n");
        for (size_t i = 0; i < nCALL; i++) fprintf(g, "C %d %lu %lu\n", CALL[i].c, (unsigned long)CALL[i].s, (unsigned long)CALL[i].e);
        for (size_t i = 0; i < nEV; i++) if (okv[i]) fprintf(g, "E %ld %ld\n", (long)((int64_t)ES[i] + off), (long)((int64_t)EE[i] + off));
        fclose(g);
    }
}
