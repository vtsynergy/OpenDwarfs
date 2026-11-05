#define LIMIT -999
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "needle.h"
#include <stdint.h>
#include <sys/time.h>
#include <omp.h>
#include "dwarf_timer.h"
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
#pragma omp declare target
static inline int maximum_d(int a, int b, int c) {
	int k;
	if (a <= b) k = b; else k = a;
	if (k <= c) return c; else return k;
}
#pragma omp end declare target

/* needle_opencl_shared_1 / needle_opencl_shared_2 (needle_kernel.cl).
 * One work-group = one BLOCK_SIZE x BLOCK_SIZE tile processed in a wavefront of BLOCK_SIZE
 * work-items separated by barriers. phase 1: tiles on/above the anti-diagonal (i grows),
 * phase 2: tiles below it. A work-group is a team, a work-item a thread (strided over however many
 * threads/teams the runtime grants; both the thread loop and the group loop are barrier-uniform). */
static void needle_kernel(const int *referrence, int *matrix, int cols, int penalty, int i, int block_width,
                          int phase, size_t size) {
#pragma omp target teams num_teams(i) thread_limit(BLOCK_SIZE) map(to : referrence[0 : size]) map(tofrom : matrix[0 : size])
	{
#pragma omp parallel num_threads(BLOCK_SIZE)
		{
			int tid = omp_get_thread_num(), nt = omp_get_num_threads();
			for (int bx = omp_get_team_num(); bx < i; bx += omp_get_num_teams()) {
				int b_index_x, b_index_y;
				if (phase == 1) { b_index_x = bx; b_index_y = i - 1 - bx; }
				else { b_index_x = bx + block_width - i; b_index_y = block_width - bx - 1; }
				for (int m = 0; m < BLOCK_SIZE; m++) {
					for (int tx = tid; tx < BLOCK_SIZE; tx += nt) {
						if (tx <= m) {
							int index = cols * BLOCK_SIZE * b_index_y + BLOCK_SIZE * b_index_x + tx + (cols + 1);
							int ref_x = index + (m - tx) * cols;
							matrix[ref_x] = maximum_d(matrix[ref_x - (cols + 1)] + referrence[ref_x],
							                          matrix[ref_x - 1] - penalty, matrix[ref_x - cols] - penalty);
						}
					}
#pragma omp barrier
				}
				for (int m = BLOCK_SIZE - 2; m >= 0; m--) {
					for (int tx = tid; tx < BLOCK_SIZE; tx += nt) {
						if (tx <= m) {
							int index = cols * BLOCK_SIZE * b_index_y + BLOCK_SIZE * b_index_x + tx + (cols + 1);
							int ref_x = index + (m - tx) * cols + (cols + 1) * (BLOCK_SIZE - 1 - m);
							matrix[ref_x] = maximum_d(matrix[ref_x - (cols + 1)] + referrence[ref_x],
							                          matrix[ref_x - 1] - penalty, matrix[ref_x - cols] - penalty);
						}
					}
#pragma omp barrier
				}
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int main( int argc, char** argv) 
{
	runTest( argc, argv);
	return EXIT_SUCCESS;
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
	int *matrix = (int *)malloc(sizeof(int) * size);   /* "matrix_cuda" of the OpenCL code, host side */
	int block_width = ( max_cols - 1 )/BLOCK_SIZE;

	dwarf_init("nw", "dynamic-programming", "omp");
	dwarf_set_size(max_cols - 1, penalty);
	for (int rep = dwarf_first_rep(); rep < dwarf_reps(); ++rep) {
		memcpy(matrix, input_itemsets, sizeof(int) * size);   /* untimed: restore initial scores */
		dwarf_region_begin();
#pragma omp target enter data map(to : matrix[0 : size], referrence[0 : size])
		printf("Processing top-left matrix\n");
		//process top-left matrix: blocks in a wave-front pattern up to the main (block-level) anti-diagonal
		for(i = 1 ; i <= block_width ; i++)
			needle_kernel(referrence, matrix, max_cols, penalty, i, block_width, 1, (size_t)size);
		printf("Processing bottom-right matrix\n");
		//process bottom-right matrix: blocks in a wave-front pattern downwards the main anti-diagonal
		for(i = block_width - 1  ; i >= 1 ; i--)
			needle_kernel(referrence, matrix, max_cols, penalty, i, block_width, 2, (size_t)size);
#pragma omp target exit data map(from : matrix[0 : size]) map(delete : referrence[0 : size])
		dwarf_region_end(rep);
	}
	memcpy(output_itemsets, matrix, sizeof(int) * size);

	/* validation (not in the OpenCL app): serial fill of the same recurrence */
	{
		int *ref = (int *)malloc(sizeof(int) * size);
		memcpy(ref, input_itemsets, sizeof(int) * size);
		for (i = 1; i < max_rows; i++)
			for (j = 1; j < max_cols; j++)
				ref[i*max_cols+j] = maximum(ref[(i-1)*max_cols + j-1] + referrence[i*max_cols+j],
				                            ref[i*max_cols + j-1] - penalty, ref[(i-1)*max_cols + j] - penalty);
		long bad = 0;
		for (i = 0; i < size; i++) if (ref[i] != output_itemsets[i]) bad++;
		printf("%s: %ld cells differ from serial Needleman-Wunsch\n", bad ? "FAILURE" : "SUCCESS", bad);
		free(ref);
		dwarf_set_size(max_cols - 1, penalty);
		printf("print traceback value GPU:\n");
		for (i = max_rows - 2,  j = max_rows - 2; i>=0, j>=0;){
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
		exit(dwarf_report(bad == 0));
	}
}
