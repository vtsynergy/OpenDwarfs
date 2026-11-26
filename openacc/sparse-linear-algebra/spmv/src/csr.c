/* Main File for SPMV (aka CSR), Sparse-Linear-Algebra dwarf -- OpenACC translation of
 * opencl/sparse-linear-algebra/SPMV/src/csr.c + spmv_kernel.cl.
 *
 * out = A*x + y with A in CSR. Same input reader (read_csr), same random x/y generation, same CLI
 * (-i -v -p -a -r -w; -k is accepted and ignored, the launch shape is pinned), same kernel body
 * (one work-item per row), same serial CPU check (-a).
 * Device buffers are created outside the timed region (clCreateBuffer in the OpenCL host); the timed
 * region is, per matrix, H2D of Ap/Aj/Ax/x/y, the kernel and D2H of y on one async queue, then wait.
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <assert.h>
#include <errno.h>
#include <string.h>
#include <openacc.h>

#include "ocd_compat.h"
#include "../inc/common.h"
#include "../inc/sparse_formats.h"
#include "dwarf_timer.h"
#include "dwarf_acc.h" /* [FIX-OCL-34] pinned launch configuration */

static struct option long_options[] = {{"verbose", 0, NULL, 'v'}, {"input_file", 1, NULL, 'i'}, {"print", 0, NULL, 'p'},
                                       {"affirm", 0, NULL, 'a'}, {"repeat", 1, NULL, 'r'}, {"kernel_file", 1, NULL, 'k'},
                                       {"wg_size", 1, NULL, 'w'}, {0, 0, 0, 0}};

/* Compares N float values; counts entries differing by more than .001 (as the OpenCL app) */
static unsigned int float_array_comp(const float *control, const float *experimental, const unsigned int N,
                                     const unsigned int exec_num) {
  unsigned int j, bad = 0;
  for (j = 0; j < N; j++) {
    float diff = experimental[j] - control[j];
    if (fabsf(diff) > .001) {
      float perc = fabsf(diff / control[j]) * 100;
      fprintf(stderr, "Possible error on exec #%u, difference of %.3f (%.1f%% error) [control=%.3f, experimental=%.3f] at row %u \n",
              exec_num, diff, perc, control[j], experimental[j], j);
      bad++;
    }
  }
  return bad;
}

/* Sparse Matrix-Vector Multiply: multiplies csr matrix by x, adds y, stores in out */
static void spmv_csr_cpu(const csr_matrix *csr, const float *x, const float *y, float *out) {
  unsigned int row, jj;
  for (row = 0; row < csr->num_rows; row++) {
    float sum = y[row];
    for (jj = csr->Ap[row]; jj < csr->Ap[row + 1]; jj++) sum += csr->Ax[jj] * x[csr->Aj[jj]];
    out[row] = sum;
  }
}

/* The OpenCL kernel `csr`: one work-item per row, with the writes and the read around it. */
static void csr_enqueue(const unsigned int num_rows, unsigned int *Ap, unsigned int *Aj, float *Ax,
                        float *x, float *y, unsigned int nnz, unsigned int ncols) {
#pragma acc update device(Ap[0 : num_rows + 1], Aj[0 : nnz], Ax[0 : nnz], x[0 : ncols], y[0 : num_rows]) async(1)
#pragma acc parallel loop gang vector ACC_LAUNCH(num_rows) \
    present(Ap[0 : num_rows + 1], Aj[0 : nnz], Ax[0 : nnz], x[0 : ncols], y[0 : num_rows]) async(1)
  for (unsigned int row = 0; row < num_rows; row++) {
    float sum = y[row];
    const unsigned int row_start = Ap[row];
    const unsigned int row_end = Ap[row + 1];
#pragma acc loop seq
    for (unsigned int jj = row_start; jj < row_end; jj++) sum += Ax[jj] * x[Aj[jj]];
    y[row] = sum;
  }
#pragma acc update self(y[0 : num_rows]) async(1)
}

int main(int argc, char **argv) {
  int verbosity = 0, do_print = 0, do_affirm = 0, opt, option_index = 0;
  unsigned int num_execs = 1, num_matrices, i, ii, k, j, num_wg_sizes = 0;
  char *file_path = NULL;
  size_t *wg_sizes = NULL;
  void *tmp;
  const char *usage = "Usage: %s -i <file_path> [-v] [-p] [-a] [-r <num_execs>] [-w <wg_size>]\n";

  while ((opt = getopt_long(argc, argv, "::vmw:k:i:par:::", long_options, &option_index)) != -1) {
    switch (opt) {
      case 'v': verbosity++; break;
      case 'i': file_path = optarg ? optarg : argv[optind]; printf("Reading Input from '%s'\n", file_path); break;
      case 'p': do_print = 1; break;
      case 'a': do_affirm = 1; break;
      case 'r': num_execs = atoi(optarg ? optarg : argv[optind]); printf("Executing %d times\n", num_execs); break;
      case 'k': printf("Note: -k ignored (single OpenACC kernel)\n"); break;
      case 'w':
        num_wg_sizes++;
        tmp = realloc(wg_sizes, sizeof(size_t) * num_wg_sizes);
        check(tmp != NULL, "csr.main() - Heap Overflow! Cannot allocate space for wg_sizes");
        wg_sizes = tmp;
        wg_sizes[num_wg_sizes - 1] = atoi(optarg ? optarg : argv[optind]);
        break;
      default: fprintf(stderr, usage, argv[0]); exit(EXIT_FAILURE);
    }
  }
  if (!file_path) { fprintf(stderr, "-i Option must be supplied\n\n"); fprintf(stderr, usage, argv[0]); exit(EXIT_FAILURE); }

  csr_matrix *csr = read_csr(&num_matrices, file_path);
  if (do_print) print_csr_arr_std(csr, num_matrices, stdout);
  else if (verbosity) { printf("Number of input matrices: %d\nMatrix 0 Metadata:\n", num_matrices); print_csr_metadata(&csr[0], stdout); }

  float *x_host = NULL, *y_host = NULL, **device_out = malloc(sizeof(float *) * num_matrices), *host_out = NULL;
  unsigned int max_row_len = 0, max_col_len = 0;
  for (ii = 0; ii < num_matrices; ii++) {
    device_out[ii] = float_new_array(csr[ii].num_rows, "csr.main() - Heap Overflow! Cannot Allocate Space for device_out");
    if (max_row_len < csr[ii].num_rows) {
      max_row_len = csr[ii].num_rows;
      y_host = float_array_realloc(y_host, csr[ii].num_rows, "csr.main() - Heap Overflow! Cannot Allocate Space for y_host");
      if (do_affirm) { host_out = realloc(host_out, sizeof(float) * max_row_len); check(host_out != NULL, "Heap Overflow! host_out"); }
    }
    if (max_col_len < csr[ii].num_cols) {
      max_col_len = csr[ii].num_cols;
      x_host = float_array_realloc(x_host, csr[ii].num_cols, "csr.main() - Heap Overflow! Cannot Allocate Space for x_host");
    }
  }
  for (ii = 0; ii < max_col_len; ii++) { x_host[ii] = rand() / (RAND_MAX + 1.0); if (do_print) printf("x[%d] = %6.2f\n", ii, x_host[ii]); }
  for (ii = 0; ii < max_row_len; ii++) { y_host[ii] = rand() / (RAND_MAX + 2.0); if (do_print) printf("y[%d] = %6.2f\n", ii, y_host[ii]); }
  if (verbosity) printf("Input Generated.\n");

  /* Default work-group size: the largest that divides the global size into the fewest work-groups
   * (same rule as default_wg_sizes() in the OpenCL host code, with the device limit taken as 1024). */
  if (!wg_sizes) {
    size_t max_wg = 1024, global_size = csr[0].num_rows;
    num_wg_sizes = 1;
    wg_sizes = malloc(sizeof(size_t));
    wg_sizes[0] = max_wg;
    size_t num_wg = global_size / wg_sizes[0];
    if (global_size % wg_sizes[0] != 0) { num_wg++; wg_sizes[0] = global_size / num_wg; }
    if (wg_sizes[0] < 1) wg_sizes[0] = 1;
  }

  unsigned int bad_total = 0;
  acc_init(acc_get_device_type()); /* device initialisation, outside the timed region */
  dwarf_init("spmv", "sparse-linear-algebra", "openacc");
  dwarf_set_size(csr[0].num_rows, csr[0].num_nonzeros);
  /* device buffers, outside the timed region */
  for (k = 0; k < num_matrices; k++) {
    unsigned int *Ap = csr[k].Ap, *Aj = csr[k].Aj, nr = csr[k].num_rows, nnz = csr[k].num_nonzeros, nc = csr[k].num_cols;
    float *Ax = csr[k].Ax, *y = device_out[k];
#pragma acc enter data create(Ap[0 : nr + 1], Aj[0 : nnz], Ax[0 : nnz], x_host[0 : nc], y[0 : nr])
  }
  printf("Kernel #1: 'csr'\n\n");
  for (ii = 0; ii < num_wg_sizes; ii++) {
    printf("Executing with WG Size #%d of %d: %d...\n", ii + 1, num_wg_sizes, (int)wg_sizes[ii]);
    /* -r N : same data N times (the OpenCL host code's repeat loop); DWARF_REPEATS overrides if larger */
    if (num_execs > 1) dwarf_set_repeats((int)num_execs);
    for (int rep = dwarf_first_rep(); rep < dwarf_reps(); ++rep) {
      for (k = 0; k < num_matrices; k++) memcpy(device_out[k], y_host, sizeof(float) * csr[k].num_rows); /* untimed */
      dwarf_region_begin();
      for (k = 0; k < num_matrices; k++) {
        if (verbosity >= 2) printf("Enqueuing Matrix #%d of %d into pipeline...\n", k + 1, num_matrices);
        csr_enqueue(csr[k].num_rows, csr[k].Ap, csr[k].Aj, csr[k].Ax, x_host, device_out[k], csr[k].num_nonzeros,
                    csr[k].num_cols);
      }
#pragma acc wait
      dwarf_region_end(rep);
      if (rep >= 0 || 1) {
        for (k = 0; k < num_matrices && do_print; k++) {
          printf("\nMatrix #%d of %d:\n", k + 1, num_matrices);
          for (j = 0; j < csr[k].num_rows; j++) printf("\trow: %d	output: %6.2f \n", j, device_out[k][j]);
        }
        if (do_affirm) {
          if (verbosity) printf("Validating results with serial C code on CPU...\n");
          for (k = 0; k < num_matrices; k++) {
            spmv_csr_cpu(&csr[k], x_host, y_host, host_out);
            bad_total += float_array_comp(host_out, device_out[k], csr[k].num_rows, (unsigned)(rep + 1));
          }
        }
      }
    }
  }
  if (do_affirm) printf("%s: %u rows differ by more than 0.001\n", bad_total ? "FAILURE" : "SUCCESS", bad_total);
  (void)i;
  return dwarf_report(bad_total == 0);
}
