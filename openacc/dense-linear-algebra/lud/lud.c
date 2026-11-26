/*
 * LU Decomposition -- OpenACC translation of opencl/dense-linear-algebra/lud
 * (Rodinia-derived blocked LUD; lud.c + lud_kernel.cl).
 *
 * Same host flow, same three kernels with BLOCK_SIZE=16 and the same indexing and work-group shapes:
 *   lud_diagonal  : 1 work-group of 16          -> 1 gang, vector length 16; barrier = end of a vector loop
 *   lud_perimeter : groups of 2*16 (no barrier) -> gang per group, vector length 32
 *   lud_internal  : G x G groups of 16 x 16     -> gang per tile, vector length 256
 * On multicore the gangs are host threads and a whole work-group runs in one thread, wave by wave.
 * Same input format (-i file), same -v verification output (including the LU matrix).
 * Timed region: host->device copy + all kernels + device->host copy (the device buffer is created
 * outside it, as in the OpenCL host).
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <getopt.h>
#include <openacc.h>
#include "common.h"
#include "dwarf_timer.h"
#include "dwarf_acc.h"

#define BS 16
int BLOCK_SIZE = BS;

static struct option long_options[] = {
    {"input", 1, NULL, 'i'}, {"verify", 0, NULL, 'v'}, {"no-verify", 0, NULL, 'n'}, {0, 0, 0, 0}};

/* ---- kernel 1: lud_diagonal (OpenCL global = local = BLOCK_SIZE) ---------------------------- */
static void lud_diagonal(float *m, int matrix_dim, int offset) {
  size_t total = (size_t)matrix_dim * matrix_dim;
#pragma acc parallel num_gangs(1) vector_length(BS) present(m[0 : total])
  {
    int array_offset = offset * matrix_dim + offset;
#pragma acc loop seq
    for (int i = 0; i < BS; i++) {
#pragma acc loop vector independent
      for (int lid = 0; lid < BS; lid++) {
        if (lid > i) {
#pragma acc loop seq
          for (int j = 0; j < i; j++)
            m[array_offset + lid * matrix_dim + i] -=
                m[array_offset + lid * matrix_dim + j] * m[array_offset + j * matrix_dim + i];
          m[array_offset + lid * matrix_dim + i] /= m[array_offset + i * matrix_dim + i];
        }
      }
#pragma acc loop vector independent
      for (int lid = 0; lid < BS; lid++) {
        if (lid > i) {
#pragma acc loop seq
          for (int j = 0; j < i + 1; j++)
            m[array_offset + (i + 1) * matrix_dim + lid] -=
                m[array_offset + (i + 1) * matrix_dim + j] * m[array_offset + j * matrix_dim + lid];
        }
      }
    }
  }
}

/* ---- kernel 2: lud_perimeter ((dim-offset)/B-1 work-groups of 2*B work-items) --------------- */
static void lud_perimeter(float *m, int matrix_dim, int offset) {
  size_t total = (size_t)matrix_dim * matrix_dim;
  int groups = (matrix_dim - offset) / BS - 1;
  if (groups <= 0) return;
#pragma acc parallel loop gang ACC_SHAPE(groups, 2 * BS) present(m[0 : total])
  for (int gid = 0; gid < groups; gid++) {
#pragma acc loop vector independent
    for (int lid = 0; lid < 2 * BS; lid++) {
      if (lid < BS) { /* peri-row */
        int idx = lid;
        int array_offset = offset * matrix_dim + offset;
#pragma acc loop seq
        for (int i = 1; i < BS; i++) {
#pragma acc loop seq
          for (int j = 0; j < i; j++)
            m[array_offset + i * matrix_dim + (gid + 1) * BS + idx] -=
                m[array_offset + i * matrix_dim + j] * m[array_offset + j * matrix_dim + (gid + 1) * BS + idx];
        }
      } else { /* peri-col */
        int idx = lid - BS;
        int array_offset_peri = (offset + (gid + 1) * BS) * matrix_dim + offset;
        int array_offset = offset * matrix_dim + offset;
#pragma acc loop seq
        for (int i = 0; i < BS; i++) {
#pragma acc loop seq
          for (int j = 0; j < i; j++)
            m[array_offset_peri + idx * matrix_dim + i] -=
                m[array_offset_peri + idx * matrix_dim + j] * m[array_offset + j * matrix_dim + i];
          m[array_offset_peri + idx * matrix_dim + i] /= m[array_offset + i * matrix_dim + i];
        }
      }
    }
  }
}

/* ---- kernel 3: lud_internal (G x G work-groups of B x B work-items) ------------------------- */
static void lud_internal(float *m, int matrix_dim, int offset) {
  size_t total = (size_t)matrix_dim * matrix_dim;
  int G = (matrix_dim - offset) / BS - 1;
  if (G <= 0) return;
#pragma acc parallel loop gang collapse(2) ACC_SHAPE(G * G, BS * BS) present(m[0 : total])
  for (int gy = 0; gy < G; gy++) {
    for (int gx = 0; gx < G; gx++) {
#pragma acc loop vector collapse(2) independent
      for (int ly = 0; ly < BS; ly++) {
        for (int lx = 0; lx < BS; lx++) {
          int global_row_id = offset + (gy + 1) * BS;
          int global_col_id = offset + (gx + 1) * BS;
          float sum = 0;
#pragma acc loop seq
          for (int i = 0; i < BS; i++)
            sum += m[(global_row_id + ly) * matrix_dim + offset + i] * m[(offset + i) * matrix_dim + global_col_id + lx];
          m[(global_row_id + ly) * matrix_dim + global_col_id + lx] -= sum;
        }
      }
    }
  }
}

/* The offloaded region: H2D, kernel loop of the OpenCL host code, D2H. */
static void lud_offload(float *m, int matrix_dim) {
  size_t total = (size_t)matrix_dim * matrix_dim;
  int i;
#pragma acc update device(m[0 : total])
  for (i = 0; i < matrix_dim - BS; i += BS) {
    lud_diagonal(m, matrix_dim, i);
    lud_perimeter(m, matrix_dim, i);
    lud_internal(m, matrix_dim, i);
  }
  lud_diagonal(m, matrix_dim, i);
#pragma acc update self(m[0 : total])
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

  acc_init(acc_get_device_type()); /* device initialisation, outside the timed region */
  size_t total = (size_t)matrix_dim * matrix_dim;
#pragma acc enter data create(work[0 : total])
  dwarf_init("lud", "dense-linear-algebra", "openacc");
  dwarf_set_size(matrix_dim, 0);
  for (int rep = dwarf_first_rep(); rep < dwarf_reps(); ++rep) {
    memcpy(work, mm, sizeof(float) * total); /* untimed setup */
    dwarf_region_begin();
    lud_offload(work, matrix_dim);
    dwarf_region_end(rep);
  }
#pragma acc exit data delete(work[0 : total])
  memcpy(m, work, sizeof(float) * total);

  int valid = 1;
  if (do_verify) {
    printf("After LUD\n");
    print_matrix(m, matrix_dim);
    printf(">>>Verify<<<<\nmatrix_dim: %d\n", matrix_dim);
    lud_verify(mm, m, matrix_dim);
  }
  if (quiet_check) {
    int mism = lud_check(mm, m, matrix_dim);
    printf("%s: %d differences beyond the size-scaled tolerance\n", mism ? "FAILURE" : "SUCCESS", mism);
    valid = (mism == 0);
  }
  free(mm); free(m); free(work);
  return dwarf_report(valid);
}
