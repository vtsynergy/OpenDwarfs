#define LIMIT -999
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "needle.h"
#include <sys/time.h>
#include "dwarf_cuda.h"
#include "dwarf_timer.h"
static int nw_valid = 1;
#include <malloc.h>
#include "needle_kernel.cuh"


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

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
	int
main( int argc, char** argv) 
{
	dwarf_cuda_init(&argc, &argv);
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
	int *matrix_cuda, *referrence_cuda;
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

	size = max_cols * max_rows;
	CUDA_CHECK(cudaMalloc((void**)&referrence_cuda, sizeof(int)*size));
	CUDA_CHECK(cudaMalloc((void**)&matrix_cuda, sizeof(int)*size));

	/* common timer: repetitions around the whole H2D -> wavefront kernels -> D2H pipeline (build/JIT is above).
	 * input_itemsets is never modified, so every repetition starts from the same state. */
	dwarf_init("nw", "dynamic-programming", "cuda");
	dwarf_set_size(max_cols - 1, penalty);
	int rep;
	for(rep = dwarf_first_rep(); rep < dwarf_reps(); rep++)
	{
	dwarf_region_begin();
	CUDA_CHECK(cudaMemcpy(referrence_cuda, referrence, sizeof(int)*size, cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaDeviceSynchronize());
	CUDA_CHECK(cudaMemcpy(matrix_cuda, input_itemsets, sizeof(int)*size, cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaDeviceSynchronize());

	int block_width = ( max_cols - 1 )/BLOCK_SIZE;

	printf("Processing top-left matrix\n");
	//process top-left matrix
	//Does what the 1st kernel loop does in a higher (block) level. i.e., takes care of blocks of BLOCK_SIZExBLOCK_SIZE in a wave-front pattern upwards
	//the main anti-diagonal (on block-level).
	//Each iteration takes care of 1, 2, 3, ... blocks that can be computed in parallel w/o dependencies
	//E.g. first block [0][0], then blocks [0][1] and [1][0], then [0][2], [1][1], [2][0], etc.
	for(i = 1 ; i <= block_width ; i++){
		needle_opencl_shared_1<<<i, BLOCK_SIZE>>>(referrence_cuda, matrix_cuda, max_cols, penalty, i, block_width);
		CUDA_SYNC_LAUNCH();
	}
	printf("Processing bottom-right matrix\n");
	//process bottom-right matrix
	//Does what the 2nd kernel loop does in a higher (block) level. i.e., takes care of blocks of BLOCK_SIZExBLOCK_SIZE in a wave-front pattern downwards
	//the main anti-diagonal.
	//Each iteration takes care of ..., 3, 2, 1 blocks that can be computed in parallel w/o dependencies
	for(i = block_width - 1  ; i >= 1 ; i--){
		needle_opencl_shared_2<<<i, BLOCK_SIZE>>>(referrence_cuda, matrix_cuda, max_cols, penalty, i, block_width);
		CUDA_SYNC_LAUNCH();
	}

	/* [FIX-OCL-06] was sizeof(float)*size: same width as int on all supported targets, but the buffer holds ints */
	CUDA_CHECK(cudaMemcpy(output_itemsets, matrix_cuda, sizeof(int)*size, cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaDeviceSynchronize());
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

	CUDA_CHECK(cudaFree(referrence_cuda));
	CUDA_CHECK(cudaFree(matrix_cuda));

}

