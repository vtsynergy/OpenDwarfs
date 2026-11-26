/* Baum-Welch algorithm for HMMs (graphical-models dwarf) -- OpenACC translation of
 * opencl/graphical-models/hmm/main_bwa_hmm.c + bwa_hmm_opencl.cl.
 *
 * Same structure: the host drives the algorithm step by step (calc_alpha, calc_beta, calc_gamma_sum,
 * calc_xi_sum, estimate_a/b/pi) and launches the 17 small kernels of bwa_hmm_opencl.cl, one per step;
 * dot products are formed on the device as element-wise products and summed serially on the host
 * (dot_production), exactly as in the OpenCL code. Quirks of the original are preserved on purpose:
 *   - est_a_dev reads alpha_d/beta_d at row `length` (one row past the data): the buffers get one
 *     extra row here so the read is in bounds (the OpenCL code reads past the allocation);
 *   - mat_vec_mul(trans='N') in estimate_b is a no-op (only 'n' and 't' are handled);
 *   - the observation sequence is all zeros; the run prints the log-likelihood.
 * CLI as the OpenCL app:  bwa_hmm [-n states | -s symbols | -t observations] -v n|s|t
 * Kernels are independent per element (no barriers/local memory) so each is a flat parallel loop.
 * Timed region: run_hmm_bwa() = device allocation + uploads + all iterations + downloads
 * (OpenCL's tic/toc also contains program build, which is excluded here).
 * Validation (added): the log-likelihood has a closed form for the uniform HMM used here
 * (T * log10(1/S) for an all-zero observation sequence).
 */
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include <string.h>
#include <sys/time.h>
#include <getopt.h>
#include <malloc.h>
#include <openacc.h>
#include "ocd_compat.h"
#include "dwarf_timer.h"
#include "dwarf_acc.h"   /* [FIX-OCL-34] */

#define T           1000        /* Number of static observations */
#define S           2           /* Number of static symbols */
#define N           60          /* Number of static states */
#define ITERATIONS  1           /* Number of iterations */
#define MAX_THREADS_PER_BLOCK 256
#define BLOCK_DIM 16
#define EXIT_ERROR 1
#define SDOT_BLOCK_SIZE (128)
#define SDOT_BLOCK_NUM (80)
#define MVMUL_BLOCK_SIZE (128)
#define MVMUL_BLOCK_NUM  (64)

typedef struct { int nstates; int nsymbols; float *a; float *b; float *pi; } Hmm;
typedef struct { int length; int *data; } Obs;

/* Global variables for device (host arrays that name the OpenACC present-table entries) */
static int nstates, nsymbols, length;
static int *obs;
static float *scale;
static float *a_d, *b_d, *pi_d, *alpha_d, *beta_d, *gamma_sum_d, *xi_sum_d, *c_d, *ones_n_d, *ones_s_d;
static size_t LA, LB;  /* element counts of a_d / b_d (square) and alpha_d / beta_d (with the extra row) */

#define ALPHA_LEN ((size_t)nstates * (length + 1))
/* [FIX-OCL-34] pinned launch. 1-D kernels: vector length DWARF_LS=128, ceil(n/128) gangs (OpenCL local 128, global roundup(n,128)).
 * 2-D kernels (16x16 work-groups in OpenCL, no barriers): one gang per tile, vector length 256 over the tile's work-items. */
#define TILES(n) (((n) + BLOCK_DIM - 1) / BLOCK_DIM)

/* ======================= kernels of bwa_hmm_opencl.cl (one parallel iteration = one work-item) ===== */
static void k_init_ones(float *ones_s, int nsym) {
#pragma acc parallel loop gang vector ACC_LAUNCH(nsym) present(ones_s[0 : nsym])
  for (int idx = 0; idx < nsym; idx++) ones_s[idx] = 1.0f;
}
static void k_init_alpha(float *b, float *pi, int ns, float *alpha, float *ones_n, int obs_t, int nsym) {
  size_t al = ALPHA_LEN;
#pragma acc parallel loop gang vector ACC_LAUNCH(ns) present(b[0 : (size_t)ns * nsym], pi[0 : ns], alpha[0 : al], ones_n[0 : ns])
  for (int idx = 0; idx < ns; idx++) {
    alpha[idx] = pi[idx] * b[(obs_t * ns) + idx];
    ones_n[idx] = 1.0f;
  }
}
static void k_calc_alpha(int ns, float *alpha, int offset, float *b, int obs_t, int nsym) {
  size_t al = ALPHA_LEN;
#pragma acc parallel loop gang vector ACC_LAUNCH(ns) present(alpha[0 : al], b[0 : (size_t)ns * nsym])
  for (int idx = 0; idx < ns; idx++) alpha[offset + idx] = alpha[offset + idx] * b[(obs_t * ns) + idx];
}
static void k_scale_alpha(int ns, float *alpha, int offset, float sc) {
  size_t al = ALPHA_LEN;
#pragma acc parallel loop gang vector ACC_LAUNCH(ns) present(alpha[0 : al])
  for (int idx = 0; idx < ns; idx++) alpha[offset + idx] = alpha[offset + idx] / sc;
}
static void k_init_beta(int ns, float *beta, int offset, float sc) {
  size_t al = ALPHA_LEN;
#pragma acc parallel loop gang vector ACC_LAUNCH(ns) present(beta[0 : al])
  for (int idx = 0; idx < ns; idx++) beta[offset + idx] = 1.0f / sc;
}
static void k_calc_beta(float *beta, float *b, float scale_t, int ns, int obs_t, int t, int nsym) {
  size_t al = ALPHA_LEN;
#pragma acc parallel loop gang vector ACC_LAUNCH(ns) present(beta[0 : al], b[0 : (size_t)ns * nsym])
  for (int idx = 0; idx < ns; idx++) beta[(t * ns) + idx] = beta[((t + 1) * ns) + idx] * b[(obs_t * ns) + idx] / scale_t;
}
static void k_calc_gamma(float *gamma_sum, float *alpha, float *beta, int ns, int t) {
  size_t al = ALPHA_LEN;
#pragma acc parallel loop gang vector ACC_LAUNCH(ns) present(gamma_sum[0 : ns], alpha[0 : al], beta[0 : al])
  for (int idx = 0; idx < ns; idx++) gamma_sum[idx] += alpha[(t * ns) + idx] * beta[(t * ns) + idx];
}
static void k_calc_xi(float *xi_sum, float *a, float *b, float *alpha, float *beta, float sum_ab, int ns, int obs_t, int t, int nsym) {
  size_t al = ALPHA_LEN;
  int gx = TILES(ns), gy = TILES(ns);
#pragma acc parallel loop gang collapse(2) ACC_SHAPE(gx * gy, 256) present(xi_sum[0 : (size_t)ns * ns], a[0 : (size_t)ns * ns], b[0 : (size_t)ns * nsym], alpha[0 : al], beta[0 : al])
  for (int by = 0; by < gy; by++)
    for (int bx = 0; bx < gx; bx++) {
#pragma acc loop vector collapse(2) independent
      for (int ly = 0; ly < BLOCK_DIM; ly++)
        for (int lx = 0; lx < BLOCK_DIM; lx++) {
          int idy = by * BLOCK_DIM + ly, idx = bx * BLOCK_DIM + lx;
          if (idx < ns && idy < ns)
            xi_sum[(idy * ns) + idx] += alpha[(t * ns) + idy] * a[(idy * ns) + idx] * b[(obs_t * ns) + idx] * beta[((t + 1) * ns) + idx] / sum_ab;
        }
    }
}
static void k_est_a(float *a, float *alpha, float *beta, float *xi_sum, float *gamma_sum, float sum_ab, int ns, int len) {
  size_t al = ALPHA_LEN;
  int gx = TILES(ns), gy = TILES(ns);
#pragma acc parallel loop gang collapse(2) ACC_SHAPE(gx * gy, 256) present(a[0 : (size_t)ns * ns], alpha[0 : al], beta[0 : al], xi_sum[0 : (size_t)ns * ns], gamma_sum[0 : ns])
  for (int by = 0; by < gy; by++)
    for (int bx = 0; bx < gx; bx++) {
#pragma acc loop vector collapse(2) independent
      for (int ly = 0; ly < BLOCK_DIM; ly++)
        for (int lx = 0; lx < BLOCK_DIM; lx++) {
          int idy = by * BLOCK_DIM + ly, idx = bx * BLOCK_DIM + lx;
          if (idx < ns && idy < ns)
            a[(idy * ns) + idx] = xi_sum[(idy * ns) + idx] / (gamma_sum[idy] - alpha[(len * ns) + idy] * beta[(len * ns) + idy] / sum_ab);
        }
    }
}
static void k_scale_a(float *a, float *c, int ns) {
  int gx = TILES(ns), gy = TILES(ns);
#pragma acc parallel loop gang collapse(2) ACC_SHAPE(gx * gy, 256) present(a[0 : (size_t)ns * ns], c[0 : ns])
  for (int by = 0; by < gy; by++)
    for (int bx = 0; bx < gx; bx++) {
#pragma acc loop vector collapse(2) independent
      for (int ly = 0; ly < BLOCK_DIM; ly++)
        for (int lx = 0; lx < BLOCK_DIM; lx++) {
          int idy = by * BLOCK_DIM + ly, idx = bx * BLOCK_DIM + lx;
          if (idx < ns && idy < ns) a[(idy * ns) + idx] = a[(idy * ns) + idx] / c[idy];
        }
    }
}
static void k_acc_b(float *b, float *alpha, float *beta, float sum_ab, int ns, int nsym, int obs_t, int t) {
  size_t al = ALPHA_LEN;
  int gx = TILES(ns), gy = TILES(nsym);
#pragma acc parallel loop gang collapse(2) ACC_SHAPE(gx * gy, 256) present(b[0 : (size_t)ns * nsym], alpha[0 : al], beta[0 : al])
  for (int by = 0; by < gy; by++)
    for (int bx = 0; bx < gx; bx++) {
#pragma acc loop vector collapse(2) independent
      for (int ly = 0; ly < BLOCK_DIM; ly++)
        for (int lx = 0; lx < BLOCK_DIM; lx++) {
          int idy = by * BLOCK_DIM + ly, idx = bx * BLOCK_DIM + lx;
          if (idy < nsym && idx < ns && obs_t == idy) b[(idy * ns) + idx] += alpha[(t * ns) + idx] * beta[(t * ns) + idx] / sum_ab;
        }
    }
}
static void k_est_b(float *b, float *gamma_sum, int ns, int nsym) {
  int gx = TILES(ns), gy = TILES(nsym);
#pragma acc parallel loop gang collapse(2) ACC_SHAPE(gx * gy, 256) present(b[0 : (size_t)ns * nsym], gamma_sum[0 : ns])
  for (int by = 0; by < gy; by++)
    for (int bx = 0; bx < gx; bx++) {
#pragma acc loop vector collapse(2) independent
      for (int ly = 0; ly < BLOCK_DIM; ly++)
        for (int lx = 0; lx < BLOCK_DIM; lx++) {
          int idy = by * BLOCK_DIM + ly, idx = bx * BLOCK_DIM + lx;
          if (idy < nsym && idx < ns) b[(idy * ns) + idx] = b[(idy * ns) + idx] / gamma_sum[idx];
        }
    }
}
static void k_scale_b(float *b, float *c, int ns, int nsym) {
  int gx = TILES(ns), gy = TILES(nsym);
#pragma acc parallel loop gang collapse(2) ACC_SHAPE(gx * gy, 256) present(b[0 : (size_t)ns * nsym], c[0 : ns])
  for (int by = 0; by < gy; by++)
    for (int bx = 0; bx < gx; bx++) {
#pragma acc loop vector collapse(2) independent
      for (int ly = 0; ly < BLOCK_DIM; ly++)
        for (int lx = 0; lx < BLOCK_DIM; lx++) {
          int idy = by * BLOCK_DIM + ly, idx = bx * BLOCK_DIM + lx;
          if (idx < ns && idy < nsym) {
            if (b[(idy * ns) + idx] == 0) b[(idy * ns) + idx] = 1e-10;
            else b[(idy * ns) + idx] = b[(idy * ns) + idx] / c[idx];
          }
        }
    }
}
static void k_est_pi(float *pi, float *alpha, float *beta, float sum_ab, int ns) {
  size_t al = ALPHA_LEN;
#pragma acc parallel loop gang vector ACC_LAUNCH(ns) present(pi[0 : ns], alpha[0 : al], beta[0 : al])
  for (int idx = 0; idx < ns; idx++) pi[idx] = alpha[idx] * beta[idx] / sum_ab;
}
/* s_dot_kernel_naive: element-wise products into partialSum (grid-stride loop == one pass over i) */
static void k_sdot(int n, float *A, int offA, float *B, int offB, float *partial) {
  size_t al = ALPHA_LEN;
#pragma acc parallel loop gang vector ACC_LAUNCH(n) present(A[0 : (size_t)(offA + n > 0 ? offA + n : 1)], B[0 : (size_t)(offB + n > 0 ? offB + n : 1)], partial[0 : n])
  for (int i = 0; i < n; i++) partial[i] = A[offA + i] * B[offB + i];
  (void)al;
}
/* mvm_non_kernel_naive / mvm_trans_kernel_naive */
static void k_mvm(int trans, int m, int n, float *A, int lda, float *x, int offX, float *y, int offY, size_t xylen, size_t alen) {
  int n_size, m_size;
  if (lda == m) { n_size = n; m_size = m; } else { n_size = m; m_size = n; }
  if (!trans) {
#pragma acc parallel loop gang vector ACC_LAUNCH(m_size) present(A[0 : alen], x[0 : xylen], y[0 : xylen])
    for (int i = 0; i < m_size; i++) {
      float sum = 0.f;
#pragma acc loop seq
      for (int j = 0; j < n_size; j++) sum += A[i * n_size + j] * x[j + offX];
      y[i + offY] = sum;
    }
  } else {
#pragma acc parallel loop gang vector ACC_LAUNCH(m_size) present(A[0 : alen], x[0 : xylen], y[0 : xylen])
    for (int i = 0; i < m_size; i++) {
      float sum = 0.f;
#pragma acc loop seq
      for (int j = 0; j < n_size; j++) sum += A[j * n_size + i] * x[j + offX];
      y[i + offY] = sum;
    }
  }
}

/* ======================= host algorithm (as in main_bwa_hmm.c) ======================= */
/* Dot product: products on the device, serial sum on the host. Buffers are identified by pointer so the
 * array lengths for the mapping are looked up from the global buffers. */
static size_t buflen(const float *p) {
  if (p == alpha_d || p == beta_d) return ALPHA_LEN;
  if (p == ones_n_d || p == pi_d || p == c_d || p == gamma_sum_d) return (size_t)nstates;
  if (p == ones_s_d) return (size_t)nsymbols;
  if (p == a_d || p == xi_sum_d) return LA;
  return LB;
}
static float dot_production(int n, float *paramA, int offsetA, float *paramB, int offsetB) {
  float result = 0.0f;
  int i;
  if (n <= 0) return result;
  float *partialSum = (float *)calloc(n, sizeof(float));  /* set partialSum_d to all zeros */
#pragma acc enter data create(partialSum[0 : n])
#pragma acc update device(partialSum[0 : n])
  (void)buflen;
  k_sdot(n, paramA, offsetA, paramB, offsetB, partialSum);
#pragma acc update self(partialSum[0 : n])
#pragma acc exit data delete(partialSum[0 : n])
  for (i = 0; i < n; i++) result += partialSum[i];
  free(partialSum);
  return result;
}

/* matrix x vector multiplication; trans 'n' -> non-transposed kernel, 't' -> transposed, anything else: no-op */
static void mat_vec_mul(char trans, int m, int n, float *A, int lda, float *x, int offsetX, float *y, int offsetY) {
  if ((trans != 'n') && (trans != 't')) return;
  size_t xylen = buflen(x), alen = buflen(A);
  k_mvm(trans == 't', m, n, A, lda, x, offsetX, y, offsetY, xylen, alen);
}

static float calc_alpha(void) {
  int offset_cur, offset_prev, t;
  float log_lik;
  /* Initialize alpha variables */
  k_init_alpha(b_d, pi_d, nstates, alpha_d, ones_n_d, *obs, nsymbols);
  /* Sum alpha values to get scaling factor */
  scale[0] = dot_production(nstates, alpha_d, 0, ones_n_d, 0);
  /* Scale alpha values */
  k_scale_alpha(nstates, alpha_d, 0, scale[0]);
  /* Initialize log likelihood */
  log_lik = log10(scale[0]);
  /* Calculate the rest of the alpha variables */
  for (t = 1; t < length; t++) {
    offset_prev = (t - 1) * nstates;
    offset_cur = t * nstates;
    /* Multiply A matrix by alpha(t-1) */
    mat_vec_mul('n', nstates, nstates, a_d, nstates, alpha_d, offset_prev, alpha_d, offset_cur);
    /* Calculate alpha(t) */
    k_calc_alpha(nstates, alpha_d, offset_cur, b_d, *(obs + t), nsymbols);
    /* Sum alpha values to get scaling factor */
    scale[t] = dot_production(nstates, alpha_d, offset_cur, ones_n_d, 0);
    /* Scale alpha values */
    k_scale_alpha(nstates, alpha_d, offset_cur, scale[t]);
    /* Update log likelihood */
    log_lik += log10(scale[t]);
  }
  return log_lik;
}

static int calc_beta(void) {
  int t, offset;
  /* Initialize beta variables */
  offset = ((length - 1) * nstates);
  k_init_beta(nstates, beta_d, offset, *(scale + length - 1));
  /* Calculate the rest of the beta variables */
  for (t = length - 2; t >= 0; t--) {
    /* Calculate first step of beta: B.*beta/scale */
    k_calc_beta(beta_d, b_d, *(scale + t), nstates, *(obs + t + 1), t, nsymbols);
    /* Multiply A matrix by beta(t) */
    mat_vec_mul('n', nstates, nstates, a_d, nstates, beta_d, t * nstates, beta_d, t * nstates);
  }
  return 0;
}

static void calc_gamma_sum(void) {
  int t;
  /* clEnqueueWriteBuffer(gamma_sum_d, zeros) */
  float *g = gamma_sum_d; int ns = nstates;
  memset(g, 0, sizeof(float) * ns);
#pragma acc update device(g[0 : ns])
  for (t = 0; t < length; t++) k_calc_gamma(gamma_sum_d, alpha_d, beta_d, nstates, t);
}

static int calc_xi_sum(void) {
  float sum_ab;
  int t;
  /* clEnqueueWriteBuffer(xi_sum_d, zeros) */
  float *x = xi_sum_d; size_t xl = LA;
  memset(x, 0, sizeof(float) * xl);
#pragma acc update device(x[0 : xl])
  /* Find the sum of xi variables */
  for (t = 0; t < length - 1; t++) {
    sum_ab = dot_production(nstates, alpha_d, t * nstates, beta_d, t * nstates);
    k_calc_xi(xi_sum_d, a_d, b_d, alpha_d, beta_d, sum_ab, nstates, *(obs + t + 1), t, nsymbols);
  }
  return 0;
}

static int estimate_a(void) {
  float sum_ab;
  /* Calculate denominator */
  sum_ab = dot_production(nstates, alpha_d, (length - 1) * nstates, beta_d, (length - 1) * nstates);
  /* Calculate new value of A */
  k_est_a(a_d, alpha_d, beta_d, xi_sum_d, gamma_sum_d, sum_ab, nstates, length);
  /* Sum rows of A to get scaling values */
  mat_vec_mul('t', nstates, nstates, a_d, nstates, ones_n_d, 0, c_d, 0);
  /* Normalize A matrix */
  k_scale_a(a_d, c_d, nstates);
  return 0;
}

static int estimate_b(void) {
  float sum_ab;
  int t;
  /* clEnqueueWriteBuffer(b_d, zeros) */
  float *bz = b_d; size_t bl = LB;
  memset(bz, 0, sizeof(float) * bl);
#pragma acc update device(bz[0 : bl])
  for (t = 0; t < length; t++) {
    sum_ab = dot_production(nstates, alpha_d, t * nstates, beta_d, t * nstates);
    k_acc_b(b_d, alpha_d, beta_d, sum_ab, nstates, nsymbols, *(obs + t), t);
  }
  /* Re-estimate B values */
  k_est_b(b_d, gamma_sum_d, nstates, nsymbols);
  /* Sum rows of B to get scaling values (the original passes trans='N': a no-op, kept) */
  mat_vec_mul('N', nstates, nsymbols, b_d, nstates, ones_s_d, 0, c_d, 0);
  /* Normalize B matrix */
  k_scale_b(b_d, c_d, nstates, nsymbols);
  return 0;
}

static int estimate_pi(void) {
  float sum_ab = dot_production(nstates, alpha_d, 0, beta_d, 0);
  k_est_pi(pi_d, alpha_d, beta_d, sum_ab, nstates);
  return 0;
}

/* Runs the Baum-Welch Algorithm on the supplied HMM and observation sequence */
static float run_hmm_bwa(Hmm *hmm, Obs *in_obs, int iterations, float threshold) {
  float *a, *b, *pi;
  float new_log_lik = 0, old_log_lik = 0;
  int iter;

  a = hmm->a; b = hmm->b; pi = hmm->pi;
  nsymbols = hmm->nsymbols;
  nstates = hmm->nstates;
  obs = in_obs->data;
  length = in_obs->length;
  LA = (size_t)nstates * nstates;
  LB = (size_t)nstates * nsymbols;

  scale = (float *)malloc(sizeof(float) * length);
  if (scale == 0) { fprintf(stderr, "ERROR: Host memory allocation error (scale)\n"); return EXIT_ERROR; }

  /* Allocate device memory (alpha/beta get one extra row, see the header comment) */
  a_d = (float *)calloc(LA, sizeof(float));
  b_d = (float *)calloc(LB, sizeof(float));
  pi_d = (float *)calloc(nstates, sizeof(float));
  alpha_d = (float *)calloc(ALPHA_LEN, sizeof(float));
  beta_d = (float *)calloc(ALPHA_LEN, sizeof(float));
  gamma_sum_d = (float *)calloc(nstates, sizeof(float));
  xi_sum_d = (float *)calloc(LA, sizeof(float));
  c_d = (float *)calloc(nstates, sizeof(float));
  ones_n_d = (float *)calloc(nstates, sizeof(float));
  ones_s_d = (float *)calloc(nsymbols, sizeof(float));
  size_t AL = ALPHA_LEN;
  int ns = nstates, nsym = nsymbols;
#pragma acc enter data create(a_d[0 : LA], b_d[0 : LB], pi_d[0 : ns], alpha_d[0 : AL], beta_d[0 : AL], gamma_sum_d[0 : ns], xi_sum_d[0 : LA], c_d[0 : ns], ones_n_d[0 : ns], ones_s_d[0 : nsym])

  /* Transfer device data */
  memcpy(a_d, a, sizeof(float) * LA); memcpy(b_d, b, sizeof(float) * LB); memcpy(pi_d, pi, sizeof(float) * nstates);
  {
    float *ad = a_d, *bd = b_d, *pd = pi_d; size_t la = LA, lb = LB;
#pragma acc update device(ad[0 : la], bd[0 : lb], pd[0 : ns])
  }
  /* Initialize ones array */
  k_init_ones(ones_s_d, nsymbols);

  /* Run BWA for either max iterations or until threshold is reached */
  for (iter = 0; iter < iterations; iter++) {
    new_log_lik = calc_alpha();
    if (new_log_lik == EXIT_ERROR) return EXIT_ERROR;
    if (calc_beta() == EXIT_ERROR) return EXIT_ERROR;
    calc_gamma_sum();
    if (calc_xi_sum() == EXIT_ERROR) return EXIT_ERROR;
    if (estimate_a() == EXIT_ERROR) return EXIT_ERROR;
    if (estimate_b() == EXIT_ERROR) return EXIT_ERROR;
    if (estimate_pi() == EXIT_ERROR) return EXIT_ERROR;
    /* check log_lik vs. threshold */
    if (threshold > 0 && iter > 0) {
      if (fabs(pow(10, new_log_lik) - pow(10, old_log_lik)) < threshold) break;
    }
    old_log_lik = new_log_lik;
  }

  /* Copy device variables back to host */
  {
    float *ad = a_d, *bd = b_d, *pd = pi_d; size_t la = LA, lb = LB;
#pragma acc update self(ad[0 : la], bd[0 : lb], pd[0 : ns])
  }
  memcpy(a, a_d, sizeof(float) * LA); memcpy(b, b_d, sizeof(float) * LB); memcpy(pi, pi_d, sizeof(float) * nstates);

  /* Free memory */
#pragma acc exit data delete(a_d[0 : LA], b_d[0 : LB], pi_d[0 : ns], alpha_d[0 : AL], beta_d[0 : AL], gamma_sum_d[0 : ns], xi_sum_d[0 : LA], c_d[0 : ns], ones_n_d[0 : ns], ones_s_d[0 : nsym])
  free(scale);
  free(a_d); free(b_d); free(pi_d); free(alpha_d); free(beta_d); free(gamma_sum_d); free(xi_sum_d); free(c_d); free(ones_n_d); free(ones_s_d);
  return new_log_lik;
}

static void free_vars(Hmm *hmm, Obs *o) {
  if (hmm != NULL) { if (hmm->a) free(hmm->a); if (hmm->b) free(hmm->b); if (hmm->pi) free(hmm->pi); free(hmm); }
  if (o != NULL) { if (o->data) free(o->data); free(o); }
}

static struct option size_opts[] = {{"state number", 1, NULL, 'n'}, {"symbol number", 1, NULL, 's'},
                                    {"observation number", 1, NULL, 't'}, {"varying mode", 1, NULL, 'v'}, {0, 0, 0, 0}};

int main(int argc, char *argv[]) {
  int s = S, t = T, n = N, i, opt, opt_index = 0;
  char v_model = 'n';
  printf("Starting bwa_hmm\n");
  while ((opt = getopt_long(argc, argv, "::n:s:t:v:", size_opts, &opt_index)) != -1) {
    switch (opt) {
      case 'v': v_model = optarg[0]; break;
      case 'n': n = atoi(optarg); break;
      case 's': s = atoi(optarg); break;
      case 't': t = atoi(optarg); break;
      default: fprintf(stderr, "Usage %s [-n state number | -s symbol number | -t observation number] [-v varying model]\n", argv[0]); exit(EXIT_FAILURE);
    }
  }
  /* parameters of the run according to the varying mode (as the three branches of the OpenCL main) */
  int nst = N, nsy = S, olen = T, varied = 0;
  if (v_model == 'n') { if (n >= 8000) return 0; nst = n; varied = n; }
  else if (v_model == 's') { if (s >= 8000) return 0; nsy = s; varied = s; }
  else if (v_model == 't') { if (t >= 10000) return 0; olen = t; varied = t; }
  else { fprintf(stderr, "-v n|s|t required\n"); return 1; }

  float log_lik = 0;
  acc_init(acc_get_device_type()); /* device initialisation, outside the timed region */
  dwarf_init("hmm", "graphical-models", "openacc");
  dwarf_set_size(nst, olen);
  for (int rep = dwarf_first_rep(); rep < dwarf_reps(); ++rep) {
    /* untimed: create the HMM and the all-zero observation sequence */
    Hmm *hmm = (Hmm *)malloc(sizeof(Hmm));
    Obs *o = (Obs *)malloc(sizeof(Obs));
    o->length = olen;
    o->data = (int *)malloc(sizeof(int) * olen);
    for (i = 0; i < olen; i++) o->data[i] = 0;
    hmm->nstates = nst; hmm->nsymbols = nsy;
    hmm->a = (float *)malloc(sizeof(float) * nst * nst);
    for (i = 0; i < nst * nst; i++) hmm->a[i] = 1.0f / (float)nst;
    hmm->b = (float *)malloc(sizeof(float) * nst * nsy);
    for (i = 0; i < nst * nsy; i++) hmm->b[i] = 1.0f / (float)nsy;
    hmm->pi = (float *)malloc(sizeof(float) * nst);
    for (i = 0; i < nst; i++) hmm->pi[i] = 1.0f / (float)nst;

    dwarf_region_begin();
    log_lik = run_hmm_bwa(hmm, o, ITERATIONS, 0);
    dwarf_region_end(rep);
    free_vars(hmm, o);
  }
  printf("Observations\tTime (s)\tLog_likelihood\n");
  printf("%i\t", varied);
  printf("%f\n", log_lik);

  /* validation: closed form for the uniform HMM with all-zero observations */
  float expect = (float)olen * log10f(1.0f / (float)nsy);
  int ok = fabsf(log_lik - expect) <= 1e-3f * fabsf(expect) + 1e-3f;
  printf("%s: log-likelihood %f vs analytic %f\n", ok ? "SUCCESS" : "FAILURE", log_lik, expect);
  return dwarf_report(ok);
}
