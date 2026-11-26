#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#include "dwarf_sycl.hpp"
#include "common.h"
#include "dwarf_timer.h"
static int lud_valid = 1;

#define BLOCK_SIZE 16   /* the OpenCL host passes -D BLOCK_SIZE=16 (halved only if the device limit is < 256) */
static int do_verify = 0;

/* kernels (translated from lud_kernel.cl) */
static void lud_diagonal(sycl::queue &q, float *m, int matrix_dim, int offset)
{
	q.parallel_for<class lud_diag>(sycl::nd_range<1>(BLOCK_SIZE, BLOCK_SIZE), [=](sycl::nd_item<1> it) {
		int i,j;
		int array_offset = offset*matrix_dim+offset;
		const int lid = it.get_local_id(0);
		for(i=0; i < BLOCK_SIZE; i++) {
			if (lid>i){
				for(j=0; j < i; j++)
					m[array_offset+lid*matrix_dim+i] -= m[array_offset+lid*matrix_dim+j]*m[array_offset+j*matrix_dim+i];
				m[array_offset+lid*matrix_dim+i] /= m[array_offset+i*matrix_dim+i];
			}
			sycl::group_barrier(it.get_group());
			if (lid>i){
				for(j=0; j < i+1; j++)
					m[array_offset+(i+1)*matrix_dim+lid] -= m[array_offset+(i+1)*matrix_dim+j]*m[array_offset+j*matrix_dim+lid];
			}
			sycl::group_barrier(it.get_group());
		}
	}).wait_and_throw();
}

static void lud_perimeter(sycl::queue &q, size_t groups, float *m, int matrix_dim, int offset)
{
	q.parallel_for<class lud_peri>(sycl::nd_range<1>(groups*BLOCK_SIZE*2, BLOCK_SIZE*2), [=](sycl::nd_item<1> it) {
		int i,j, array_offset;
		int array_offset_peri;
		int idx;
		const int lid = it.get_local_id(0), gid = it.get_group(0);
		if (lid < BLOCK_SIZE) { //peri-row
			idx=lid;
			array_offset = offset*matrix_dim+offset;
			for(i=1; i < BLOCK_SIZE; i++){
				for (j=0; j < i; j++)
					m[array_offset+i*matrix_dim+(gid+1)*BLOCK_SIZE+idx] -= m[array_offset+i*matrix_dim+j]*m[array_offset+j*matrix_dim+(gid+1)*BLOCK_SIZE+idx];
			}
		} else { //peri-col
			idx=lid - BLOCK_SIZE;
			array_offset_peri = (offset+(gid+1)*BLOCK_SIZE)*matrix_dim+offset;
			array_offset = offset*matrix_dim+offset;
			for(i=0; i < BLOCK_SIZE; i++){
				for(j=0; j < i; j++)
					m[array_offset_peri+idx*matrix_dim+i] -= m[array_offset_peri+idx*matrix_dim+j] * m[array_offset+j*matrix_dim+i];
				m[array_offset_peri+idx*matrix_dim+i] /= m[array_offset+i*matrix_dim+i];
			}
		}
	}).wait_and_throw();
}

/* OpenCL NDRange {n, n} local {16, 16}: OpenCL dim 0 = SYCL dim 1 */
static void lud_internal(sycl::queue &q, size_t groups, float *m, int matrix_dim, int offset)
{
	q.parallel_for<class lud_int>(sycl::nd_range<2>({groups*BLOCK_SIZE, groups*BLOCK_SIZE}, {BLOCK_SIZE, BLOCK_SIZE}), [=](sycl::nd_item<2> it) {
		int i;
		float sum;
		const int lx = it.get_local_id(1), ly = it.get_local_id(0);
		int global_row_id = offset + ((int)it.get_group(0)+1)*BLOCK_SIZE;
		int global_col_id = offset + ((int)it.get_group(1)+1)*BLOCK_SIZE;
		sum = 0;
		for (i=0; i < BLOCK_SIZE; i++)
			sum += m[(global_row_id+ly)*matrix_dim+offset+i] * m[(offset+i)*matrix_dim+global_col_id+lx];
		m[(global_row_id+ly)*matrix_dim+global_col_id+lx] -= sum;
	}).wait_and_throw();
}


static struct option long_options[] = {
	/* name, has_arg, flag, val */
	{"input", 1, NULL, 'i'},
	{"platform", 1, NULL, 'p'},
	{"device", 1, NULL, 'd'},
	{"size", 1, NULL, 's'},
	{"verify", 0, NULL, 'v'},
	{0,0,0,0}
};

	int
main ( int argc, char *argv[] )
{
	int matrix_dim = 32; /* default matrix_dim */
	int opt, option_index=0;
	func_ret_t ret;
	const char *input_file = NULL;
	float *m, *mm;
	stopwatch sw;

	float *d_m;

	dsycl::parse_args(&argc, &argv);
	sycl::queue q = dsycl::make_queue();

	while ((opt = getopt_long(argc, argv, "::vs:i:", 
					long_options, &option_index)) != -1 ) {
		switch(opt){
			case 'i':
				input_file = optarg;
				break;
			case 'v':
				do_verify = 1;
				break;
			case 's':
				matrix_dim = atoi(optarg);
				fprintf(stderr, "Currently not supported, use -i instead\n");
				fprintf(stderr, "Usage: %s [-v] [-s matrix_size|-i input_file]\n", argv[0]);
				exit(EXIT_FAILURE);
			case '?':
				fprintf(stderr, "invalid option\n");
				break;
			case ':':
				fprintf(stderr, "missing argument\n");
				//break;
			default:
				fprintf(stderr, "Usage: %s [-v] [-s matrix_size|-i input_file]\n",
						argv[0]);
				exit(EXIT_FAILURE);
		}
	}

	if ( (optind < argc) || (optind == 1)) {
		fprintf(stderr, "Usage: %s [-v] [-s matrix_size|-i input_file]\n", argv[0]);
		exit(EXIT_FAILURE);
	}

	if (input_file) {
		printf("Reading matrix from file %s\n", input_file);
		ret = create_matrix_from_file(&m, input_file, &matrix_dim);
		if (ret != RET_SUCCESS) {
			m = NULL;
			fprintf(stderr, "error create matrix from file %s\n", input_file);
			exit(EXIT_FAILURE);
		}
	} else {
		printf("No input file specified!\n");
		exit(EXIT_FAILURE);
	}

	if (do_verify){
		printf("Before LUD\n");
		print_matrix(m, matrix_dim);
		matrix_duplicate(m, &mm, matrix_dim);
	}

	d_m = dsycl::alloc<float>(q, (size_t)matrix_dim*matrix_dim);
	dsycl::prebuild(q);   /* kernel JIT, as clBuildProgram: outside the timed region */

	/* Common timer API (common/dwarf_timer.h): the timed region is the same span as the original
	 * stopwatch -- H2D write, all kernels, D2H read; build/JIT happened above and is excluded. */
	float *m_orig = NULL; int rep; int i;
	matrix_duplicate(m, &m_orig, matrix_dim);
	dwarf_init("lud", "dense-linear-algebra", "sycl");
	dwarf_set_size(matrix_dim, 0);
	for (rep = dwarf_first_rep(); rep < dwarf_reps(); ++rep) {
	memcpy(m, m_orig, (size_t)matrix_dim*matrix_dim*sizeof(float));
	stopwatch_start(&sw);
	dwarf_region_begin();

	dsycl::write(q, d_m, m, (size_t)matrix_dim*matrix_dim*sizeof(float));

	i=0;
#ifdef START_POWER
	for( int iter = 0; iter < 1000; iter++)
#endif
		for (i=0; i < matrix_dim-BLOCK_SIZE; i += BLOCK_SIZE) {
			lud_diagonal(q, d_m, matrix_dim, i);
			lud_perimeter(q, (matrix_dim-i)/BLOCK_SIZE-1, d_m, matrix_dim, i);
			lud_internal(q, (matrix_dim-i)/BLOCK_SIZE-1, d_m, matrix_dim, i);
		}
	lud_diagonal(q, d_m, matrix_dim, i);

	dsycl::read(q, m, d_m, (size_t)matrix_dim*matrix_dim*sizeof(float));
		/* end of timing point */
		stopwatch_stop(&sw);
	dwarf_region_end(rep);
	printf("Time consumed(ms): %lf\n", 1000*get_interval_by_sec(&sw));
	}

	dsycl::release(q, d_m);

	if (do_verify){
		printf("After LUD\n");
		print_matrix(m, matrix_dim);
		printf(">>>Verify<<<<\n");
		printf("matrix_dim: %d\n",matrix_dim);
		lud_verify(mm, m, matrix_dim); 
		free(mm);
	}

	/* [FIX-OCL-32] validity check on every run (previously only with -v, which also prints the whole matrices): quiet L*U check */
	{
		int mism = lud_check(m_orig, m, matrix_dim);
		printf("%s: %d differences beyond the size-scaled tolerance\n", mism ? "FAILURE" : "SUCCESS", mism);
		lud_valid = (mism == 0);
	}
	free(m);
	free(m_orig);
	return dwarf_report(lud_valid) == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}				/* ----------  end of function main  ---------- */
