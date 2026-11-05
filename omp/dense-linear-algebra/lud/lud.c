/*
 * LU Decomposition -- OpenMP target-offload translation of opencl/dense-linear-algebra/lud
 * (Rodinia-derived blocked LUD).
 *
 * Fidelity: same host flow, same three kernels (lud_diagonal / lud_perimeter / lud_internal)
 * with the same BLOCK_SIZE=16 work-group structure and indexing. Each OpenCL work-group maps to
 * an OpenMP team and each work-item to a thread of that team; barrier(...) becomes
 * `#pragma omp barrier`. Same input format (-i file), same -v verification output.
 * The OpenCL device/platform plumbing (ocd_*) is the only part that is not carried over.
 *
 * Timed region (common/dwarf_timer.h): device allocation + host->device copy + all kernels +
 * device->host copy -- the same span as the OpenCL stopwatch (write buffer .. read buffer).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <omp.h>
#include "common.h"
#include "dwarf_timer.h"

int BLOCK_SIZE = 16;

static struct option long_options[] = {
    {"input", 1, NULL, 'i'}, {"verify", 0, NULL, 'v'}, {"no-verify", 0, NULL, 'n'}, {0, 0, 0, 0}};

/* ---- kernel 1: lud_diagonal  (OpenCL: global=local=BLOCK_SIZE, 1 work-group) -------------- */
static void lud_diagonal(float *m, int matrix_dim, int offset, int B) {
  size_t total = (size_t)matrix_dim * matrix_dim;
#pragma omp target teams num_teams(1) thread_limit(B) map(tofrom : m[0 : total]) firstprivate(matrix_dim, offset, B)
  {
#pragma omp parallel num_threads(B)
    {
      /* A work-group of B work-items is emulated by however many threads the team really gets
       * (OpenMP may grant fewer than requested): each thread covers lid = tid, tid+nt, ... .
       * Both phases of an iteration are data-parallel over lid and are separated by the same
       * barriers as in the OpenCL kernel, so the result does not depend on the thread count. */
      int tid = omp_get_thread_num(), nt = omp_get_num_threads();
      int lid, i, j;
      int array_offset = offset * matrix_dim + offset;
      for (i = 0; i < B; i++) {
        for (lid = tid; lid < B; lid += nt) {
          if (lid > i) {
            for (j = 0; j < i; j++)
              m[array_offset + lid * matrix_dim + i] -=
                  m[array_offset + lid * matrix_dim + j] * m[array_offset + j * matrix_dim + i];
            m[array_offset + lid * matrix_dim + i] /= m[array_offset + i * matrix_dim + i];
          }
        }
#pragma omp barrier
        for (lid = tid; lid < B; lid += nt) {
          if (lid > i) {
            for (j = 0; j < i + 1; j++)
              m[array_offset + (i + 1) * matrix_dim + lid] -=
                  m[array_offset + (i + 1) * matrix_dim + j] * m[array_offset + j * matrix_dim + lid];
          }
        }
#pragma omp barrier
      }
    }
  }
}

/* ---- kernel 2: lud_perimeter  (OpenCL: (dim-offset)/B-1 work-groups of 2*B work-items) ---- */
static void lud_perimeter(float *m, int matrix_dim, int offset, int B) {
  size_t total = (size_t)matrix_dim * matrix_dim;
  int groups = (matrix_dim - offset) / B - 1;
  if (groups <= 0) return;
#pragma omp target teams num_teams(groups) thread_limit(2 * B) map(tofrom : m[0 : total]) \
    firstprivate(matrix_dim, offset, B)
  {
    /* num_teams() is only an upper bound in OpenMP: cover work-group ids by striding over the
     * teams actually created, exactly like the thread striding below. */
    for (int gid = omp_get_team_num(); gid < groups; gid += omp_get_num_teams())
#pragma omp parallel num_threads(2 * B)
    {
      int tid = omp_get_thread_num(), nt = omp_get_num_threads();
      int lid, i, j, idx, array_offset, array_offset_peri;
      for (lid = tid; lid < 2 * B; lid += nt) { /* work-item lid of the 2*B-item work-group */
        if (lid < B) { /* peri-row */
          idx = lid;
          array_offset = offset * matrix_dim + offset;
          for (i = 1; i < B; i++) {
            for (j = 0; j < i; j++)
              m[array_offset + i * matrix_dim + (gid + 1) * B + idx] -=
                  m[array_offset + i * matrix_dim + j] * m[array_offset + j * matrix_dim + (gid + 1) * B + idx];
          }
        } else { /* peri-col */
          idx = lid - B;
          array_offset_peri = (offset + (gid + 1) * B) * matrix_dim + offset;
          array_offset = offset * matrix_dim + offset;
          for (i = 0; i < B; i++) {
            for (j = 0; j < i; j++)
              m[array_offset_peri + idx * matrix_dim + i] -=
                  m[array_offset_peri + idx * matrix_dim + j] * m[array_offset + j * matrix_dim + i];
            m[array_offset_peri + idx * matrix_dim + i] /= m[array_offset + i * matrix_dim + i];
          }
        }
      }
    }
  }
}

/* ---- kernel 3: lud_internal  (OpenCL: G x G work-groups of B x B work-items) -------------- */
static void lud_internal(float *m, int matrix_dim, int offset, int B) {
  size_t total = (size_t)matrix_dim * matrix_dim;
  int G = (matrix_dim - offset) / B - 1;
  if (G <= 0) return;
#pragma omp target teams distribute collapse(2) map(tofrom : m[0 : total]) firstprivate(matrix_dim, offset, B, G)
  for (int gy = 0; gy < G; gy++) {
    for (int gx = 0; gx < G; gx++) {
      int global_row_id = offset + (gy + 1) * B;
      int global_col_id = offset + (gx + 1) * B;
#pragma omp parallel for collapse(2) num_threads(B * B)
      for (int ly = 0; ly < B; ly++) {
        for (int lx = 0; lx < B; lx++) {
          float sum = 0;
          for (int i = 0; i < B; i++)
            sum += m[(global_row_id + ly) * matrix_dim + offset + i] * m[(offset + i) * matrix_dim + global_col_id + lx];
          m[(global_row_id + ly) * matrix_dim + global_col_id + lx] -= sum;
        }
      }
    }
  }
}

/* The offloaded region: H2D, kernel loop of the OpenCL host code, D2H. */
static void lud_offload(float *m, int matrix_dim, int B) {
  size_t total = (size_t)matrix_dim * matrix_dim;
  int i;
#pragma omp target enter data map(to : m[0 : total])
  for (i = 0; i < matrix_dim - B; i += B) {
    lud_diagonal(m, matrix_dim, i, B);
    lud_perimeter(m, matrix_dim, i, B);
    lud_internal(m, matrix_dim, i, B);
  }
  lud_diagonal(m, matrix_dim, i, B);
#pragma omp target exit data map(from : m[0 : total])
}

/* Quiet equivalent of lud_verify() in common.c (same L*U reconstruction and 1e-4 tolerance). */
static int lud_check(const float *m, const float *lu, int dim) {
  int mismatches = 0, i, j, k;
  for (i = 0; i < dim; i++)
    for (j = 0; j < dim; j++) {
      float sum = 0, l, u;
      for (k = 0; k <= MIN(i, j); k++) {
        l = (i == k) ? 1 : lu[i * dim + k];
        u = lu[k * dim + j];
        sum += l * u;
      }
      if (fabs(m[i * dim + j] - sum) > 0.0001) mismatches++;
    }
  return mismatches;
}

int main(int argc, char *argv[]) {
  int matrix_dim = 32, opt, option_index = 0, do_verify = 0, quiet_check = 1;
  const char *input_file = NULL;
  float *m, *mm, *work;

  while ((opt = getopt_long(argc, argv, "::vni:", long_options, &option_index)) != -1) {
    switch (opt) {
      case 'i': input_file = optarg; break;
      case 'v': do_verify = 1; break;
      case 'n': quiet_check = 0; break;
      default: fprintf(stderr, "Usage: %s [-v] [--no-verify] -i input_file\n", argv[0]); exit(EXIT_FAILURE);
    }
  }
  if (!input_file) { fprintf(stderr, "Usage: %s [-v] [--no-verify] -i input_file\n", argv[0]); exit(EXIT_FAILURE); }
  printf("Reading matrix from file %s\n", input_file);
  if (create_matrix_from_file(&m, input_file, &matrix_dim) != RET_SUCCESS) {
    fprintf(stderr, "error create matrix from file %s\n", input_file);
    exit(EXIT_FAILURE);
  }
  matrix_duplicate(m, &mm, matrix_dim); /* pristine copy: reference for verification & per-rep restore */
  work = (float *)malloc(sizeof(float) * (size_t)matrix_dim * matrix_dim);

  /* Same as the OpenCL host code: shrink BLOCK_SIZE until BLOCK_SIZE^2 fits the device work-group limit. */
  {
    int max_wg = 1024; /* NVIDIA / typical limit; CPU target has no hard limit */
    while (BLOCK_SIZE * BLOCK_SIZE > max_wg) BLOCK_SIZE /= 2;
  }

  omp_set_max_active_levels(2); /* host fallback: teams -> parallel needs two active levels (work-group -> work-item) */
  dwarf_init("lud", "dense-linear-algebra", "omp");
  dwarf_set_size(matrix_dim, 0);
  for (int rep = dwarf_first_rep(); rep < dwarf_reps(); ++rep) {
    memcpy(work, mm, sizeof(float) * (size_t)matrix_dim * matrix_dim); /* untimed setup */
    dwarf_region_begin();
    lud_offload(work, matrix_dim, BLOCK_SIZE);
    dwarf_region_end(rep);
  }
  memcpy(m, work, sizeof(float) * (size_t)matrix_dim * matrix_dim);

  int valid = 1;
  if (do_verify) { printf(">>>Verify<<<<\nmatrix_dim: %d\n", matrix_dim); lud_verify(mm, m, matrix_dim); }
  if (quiet_check) {
    int mism = lud_check(mm, m, matrix_dim);
    printf("%s: %d differences of more than 0.0001f (absolute)\n", mism ? "FAILURE" : "SUCCESS", mism);
    valid = (mism == 0);
  }
  free(mm); free(m); free(work);
  return dwarf_report(valid);
}
