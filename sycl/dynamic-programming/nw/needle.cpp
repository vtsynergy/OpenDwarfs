#define LIMIT -999
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "needle.h"
#ifndef BLOCK_SIZE
#define BLOCK_SIZE 16   /* needle_kernel.cl */
#endif
#include <sys/time.h>
#include "dwarf_sycl.hpp"
#include "dwarf_timer.h"
static int nw_valid = 1;
#include <malloc.h>
#define AOCL_ALIGNMENT 64


//#define TRACE

void runTest( int argc, char** argv);

int blosum62[24][24] = {
	{ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4},
	{-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4},
	{-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4},
	{-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4},
	{ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4},
	{-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4},
	{-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
	{ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4},
	{-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4},
	{-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4},
	{-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4},
	{-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4},
	{-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4},
	{-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4},
	{-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4},
	{ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4},
	{ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4},
	{-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4},
	{-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4},
	{ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4},
	{-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4},
	{-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
	{ 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4},
	{-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1}
};

double gettime() {
	struct timeval t;
	gettimeofday(&t,NULL);
	return t.tv_sec+t.tv_usec*1e-6;
}

int 
maximum( int a,
		int b,
		int c){

	int k;
	if( a <= b )
		k = b;
	else 
		k = a;

	if( k <=c )
		return(c);
	else
		return(k);

}

/* kernels (translated from needle_kernel.cl). OpenCL NDRange {i*16, 1} local {16, 1}: OpenCL dim 0 = SYCL dim 1.
 * which = 1: needle_opencl_shared_1 (top-left), which = 2: needle_opencl_shared_2 (bottom-right). */
template <int which>
static void needle_kernel(sycl::queue &q, size_t groups, int* referrence, int* matrix_opencl, int cols, int penalty, int i, int block_width)
{
	q.parallel_for(sycl::nd_range<2>({1, groups*BLOCK_SIZE}, {1, BLOCK_SIZE}), [=](sycl::nd_item<2> it) {
		int bx = it.get_group(1);
		int tx = it.get_local_id(1);

		int b_index_x = which == 1 ? bx : bx + block_width - i;
		int b_index_y = which == 1 ? i - 1 - bx : block_width - bx - 1;
		int index   = cols * BLOCK_SIZE * b_index_y + BLOCK_SIZE * b_index_x + tx + ( cols + 1 );

		for( int m = 0 ; m < BLOCK_SIZE ; m++){
			if ( tx <= m ){
				int ref_x=index+(m-tx)*cols;
				matrix_opencl[ref_x] = maximum( matrix_opencl[ref_x-(cols+1)] + referrence[ref_x],
						matrix_opencl[ref_x-1]  - penalty,
						matrix_opencl[ref_x-cols]  - penalty);
			}
			sycl::group_barrier(it.get_group());
		}

		for( int m = BLOCK_SIZE - 2 ; m >=0 ; m--){
			if ( tx <= m){
				int ref_x=index+(m-tx)*cols+(cols+1)*(BLOCK_SIZE-1-m);
				matrix_opencl[ref_x] = maximum( matrix_opencl[ref_x-(cols+1)] + referrence[ref_x],
						matrix_opencl[ref_x-1]  - penalty,
						matrix_opencl[ref_x-cols]  - penalty);
			}
			sycl::group_barrier(it.get_group());
		}
	}).wait_and_throw();
}

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
	int
main( int argc, char** argv) 
{
	dsycl::parse_args(&argc, &argv);
	runTest( argc, argv);
	return dwarf_report(nw_valid) == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}

void usage(int argc, char **argv)
{
	fprintf(stderr, "Usage: %s <max_rows/max_cols> <penalty> \n", argv[0]);
	fprintf(stderr, "\t<dimension>  - x and y dimensions\n");
	fprintf(stderr, "\t<penalty> - penalty(positive integer)\n");
	exit(1);
}

void runTest( int argc, char** argv) 
{
	int max_rows, max_cols, penalty;
	int *input_itemsets, *output_itemsets, *referrence;
	int *matrix_cuda,*referrence_cuda;
	int size;
	int i, j;

	// the lengths of the two sequences should be able to divided by 16.
	// And at current stage  max_rows needs to equal max_cols
	if (argc == 3)
	{
		max_rows = atoi(argv[1]);
		max_cols = atoi(argv[1]);
		penalty = atoi(argv[2]);
	}
	else{
		usage(argc, argv);
	}

	if(atoi(argv[1])%16!=0){
		fprintf(stderr,"The dimension values must be a multiple of 16\n");
		exit(1);
	}

	max_rows = max_rows + 1;
	max_cols = max_cols + 1;
	referrence      = (int *)malloc( max_rows * max_cols * sizeof(int) );
	input_itemsets  = (int *)malloc(max_rows * max_cols * sizeof(int) );
	output_itemsets = (int *)malloc(max_rows * max_cols * sizeof(int) );

	if (!input_itemsets)
		fprintf(stderr, "error: can not allocate memory");

	srand ( 7 );
	//Initialize to zero
	for (i = 0 ; i < max_cols; i++){
		for (j = 0 ; j < max_rows; j++){
			input_itemsets[i*max_cols+j] = 0;
		}
	}

	printf("Start Needleman-Wunsch\n");
	//Used for input itemsets, then after the referrence matrix is constructed it is used/rewritten to store scores.
	for(i=1; i< max_rows ; i++){    //please define your own sequence. 
		input_itemsets[i*max_cols] = rand() % 10 + 1;
	}
	for(j=1; j< max_cols ; j++){    //please define your own sequence.
		input_itemsets[j] = rand() % 10 + 1;
	}


	for (i = 1 ; i < max_cols; i++){
		for (j = 1 ; j < max_rows; j++){
			referrence[i*max_cols+j] = blosum62[input_itemsets[i*max_cols]][input_itemsets[j]];
		}
	}

	//Fill first row and first column with initial scores -10, -20itouch, -30, ... (upper left corner is set to zero).
	for(i = 1; i< max_rows ; i++)
		input_itemsets[i*max_cols] = -i * penalty;
	for(j = 1; j< max_cols ; j++)
		input_itemsets[j] = -j * penalty;

	sycl::queue q = dsycl::make_queue();

	size = max_cols * max_rows;
	referrence_cuda = dsycl::alloc<int>(q, size);
	matrix_cuda = dsycl::alloc<int>(q, size);
	dsycl::prebuild(q);   /* kernel JIT, as clBuildProgram: outside the timed region */

	/* common timer: repetitions around the whole H2D -> wavefront kernels -> D2H pipeline (build/JIT is above).
	 * input_itemsets is never modified, so every repetition starts from the same state. */
	dwarf_init("nw", "dynamic-programming", "sycl");
	dwarf_set_size(max_cols - 1, penalty);
	int rep;
	for(rep = dwarf_first_rep(); rep < dwarf_reps(); rep++)
	{
	dwarf_region_begin();
	dsycl::write(q, referrence_cuda, referrence, sizeof(int)*size);
	dsycl::write(q, matrix_cuda, input_itemsets, sizeof(int)*size);

	int block_width = ( max_cols - 1 )/BLOCK_SIZE;

	printf("Processing top-left matrix\n");
	for(i = 1 ; i <= block_width ; i++)
		needle_kernel<1>(q, i, referrence_cuda, matrix_cuda, max_cols, penalty, i, block_width);
	printf("Processing bottom-right matrix\n");
	for(i = block_width - 1  ; i >= 1 ; i--)
		needle_kernel<2>(q, i, referrence_cuda, matrix_cuda, max_cols, penalty, i, block_width);

	dsycl::read(q, output_itemsets, matrix_cuda, sizeof(int)*size);
	dwarf_region_end(rep);
	}

	/* [FIX-OCL-07] the app had no correctness check. Added: serial fill of the same recurrence compared cell by cell. */
	{
		int *ref = (int *)malloc(sizeof(int) * size); long bad = 0;
		memcpy(ref, input_itemsets, sizeof(int) * size);
		for (i = 1; i < max_rows; i++)
			for (j = 1; j < max_cols; j++)
				ref[i*max_cols+j] = maximum(ref[(i-1)*max_cols + j-1] + referrence[i*max_cols+j], ref[i*max_cols + j-1] - penalty, ref[(i-1)*max_cols + j] - penalty);
		for (i = 0; i < size; i++) if (ref[i] != output_itemsets[i]) bad++;
		printf("%s: %ld cells differ from serial Needleman-Wunsch\n", bad ? "FAILURE" : "SUCCESS", bad);
		nw_valid = (bad == 0);
		free(ref);
	}


//#ifdef TRACE

	printf("print traceback value GPU:\n");

	for (i = max_rows - 2,  j = max_rows - 2; i>=0 && j>=0;){ /* [FIX-OCL-08] was `i>=0, j>=0` (comma operator: only j>=0 was tested) */

		int nw, n, w, traceback;

		if ( i == max_rows - 2 && j == max_rows - 2 )
			printf("%d ", output_itemsets[ i * max_cols + j]); //print the first element


		if ( i == 0 && j == 0 )
			break;


		if ( i > 0 && j > 0 ){
			nw = output_itemsets[(i - 1) * max_cols + j - 1];
			w  = output_itemsets[ i * max_cols + j - 1 ];
			n  = output_itemsets[(i - 1) * max_cols + j];
		}
		else if ( i == 0 ){
			nw = n = LIMIT;
			w  = output_itemsets[ i * max_cols + j - 1 ];
		}
		else if ( j == 0 ){
			nw = w = LIMIT;
			n  = output_itemsets[(i - 1) * max_cols + j];
		}
		else{
		}

		traceback = maximum(nw, w, n);

		printf("%d ", traceback);

		if(traceback == nw )
		{i--; j--; continue;}

		else if(traceback == w )
		{j--; continue;}

		else if(traceback == n )
		{i--; continue;}

		else
			;
	}
	printf("\n");

//#endif

	dsycl::release(q, referrence_cuda);
	dsycl::release(q, matrix_cuda);

}

