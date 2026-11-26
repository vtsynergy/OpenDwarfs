#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "dwarf_cuda.h"
#include "lud_kernel.cuh"
#include "common.h"
#include "dwarf_timer.h"
static int lud_valid = 1;

//#define USEGPU 1
int BLOCK_SIZE = BLOCK_SIZE_K; /* OpenCL halves this only if the device max work-item size is < 256 */
static int do_verify = 0;


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

	dwarf_cuda_init(&argc, &argv);

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

	CUDA_CHECK(cudaMalloc((void**)&d_m, (size_t)matrix_dim*matrix_dim*sizeof(float)));

	/* Common timer API (common/dwarf_timer.h): the timed region is the same span as the original
	 * stopwatch -- H2D write, all kernels, D2H read; build/JIT happened above and is excluded. */
	float *m_orig = NULL; int rep; int i;
	matrix_duplicate(m, &m_orig, matrix_dim);
	dwarf_init("lud", "dense-linear-algebra", "cuda");
	dwarf_set_size(matrix_dim, 0);
	for (rep = dwarf_first_rep(); rep < dwarf_reps(); ++rep) {
	memcpy(m, m_orig, (size_t)matrix_dim*matrix_dim*sizeof(float));
	stopwatch_start(&sw);
	dwarf_region_begin();

	CUDA_CHECK(cudaMemcpy(d_m, m, (size_t)matrix_dim*matrix_dim*sizeof(float), cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaDeviceSynchronize());

	i=0;
#ifdef START_POWER
	for( int iter = 0; iter < 1000; iter++)
#endif
		for (i=0; i < matrix_dim-BLOCK_SIZE; i += BLOCK_SIZE) {
			lud_diagonal<<<1, BLOCK_SIZE>>>(d_m, matrix_dim, i);
			CUDA_SYNC_LAUNCH();
			unsigned int nb = (matrix_dim-i)/BLOCK_SIZE-1;
			lud_perimeter<<<nb, BLOCK_SIZE*2>>>(d_m, matrix_dim, i);
			CUDA_SYNC_LAUNCH();
			lud_internal<<<dim3(nb, nb), dim3(BLOCK_SIZE, BLOCK_SIZE)>>>(d_m, matrix_dim, i);
			CUDA_SYNC_LAUNCH();
		}
	lud_diagonal<<<1, BLOCK_SIZE>>>(d_m, matrix_dim, i);
	CUDA_SYNC_LAUNCH();

	CUDA_CHECK(cudaMemcpy(m, d_m, (size_t)matrix_dim*matrix_dim*sizeof(float), cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaDeviceSynchronize());
		/* end of timing point */
		stopwatch_stop(&sw);
	dwarf_region_end(rep);
	printf("Time consumed(ms): %lf\n", 1000*get_interval_by_sec(&sw));
	}

	CUDA_CHECK(cudaFree(d_m));

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
