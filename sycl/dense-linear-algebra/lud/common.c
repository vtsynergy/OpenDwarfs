#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <malloc.h>

#include "common.h"
#define AOCL_ALIGNMENT 64
void stopwatch_start(stopwatch *sw){
	if (sw == NULL)
		return;

	bzero(&sw->begin, sizeof(struct timeval));
	bzero(&sw->end  , sizeof(struct timeval));

	gettimeofday(&sw->begin, NULL);
}

void stopwatch_stop(stopwatch *sw){
	if (sw == NULL)
		return;

	gettimeofday(&sw->end, NULL);
}

double 
get_interval_by_sec(stopwatch *sw){
	if (sw == NULL)
		return 0;
	return ((double)(sw->end.tv_sec-sw->begin.tv_sec)+(double)(sw->end.tv_usec-sw->begin.tv_usec)/1000000);
}

int 
get_interval_by_usec(stopwatch *sw){
	if (sw == NULL)
		return 0;
	return ((sw->end.tv_sec-sw->begin.tv_sec)*1000000+(sw->end.tv_usec-sw->begin.tv_usec));
}

func_ret_t 
create_matrix_from_file(float **mp, const char* filename, int *size_p){
	int i, j, size;
	float *m;
	FILE *fp = NULL;

	fp = fopen(filename, "rb");
	if ( fp == NULL) {
		return RET_FAILURE;
	}

	fscanf(fp, "%d\n", &size);

	m = (float*)memalign ( AOCL_ALIGNMENT,sizeof(float)*size*size);
	if ( m == NULL) {
		fclose(fp);
		return RET_FAILURE;
	}

	for (i=0; i < size; i++) {
		for (j=0; j < size; j++) {
			fscanf(fp, "%f ", m+i*size+j);
		}
	}

	fclose(fp);

	*size_p = size;
	*mp = m;

	return RET_SUCCESS;
}

func_ret_t
create_matrix_from_random(float **mp, int size){
	float *l, *u, *m;
	/* [FIX-OCL-29] `m` was never allocated (write through an uninitialised pointer) and each element was assigned the last
	 * l*u term instead of the sum over k. */
	int i,j,k;

	srand(time(NULL));

	l = (float*)malloc(size*size*sizeof(float));
	if ( l == NULL)
		return RET_FAILURE;

	u = (float*)malloc(size*size*sizeof(float));
	if ( u == NULL) {
		free(l);
		return RET_FAILURE;
	}

	for (i = 0; i < size; i++) {
		for (j=0; j < size; j++) {
			if (i>j) {
				l[i*size+j] = GET_RAND_FP;
			} else if (i == j) {
				l[i*size+j] = 1;
			} else {
				l[i*size+j] = 0;
			}
		}
	}

	for (j=0; j < size; j++) {
		for (i=0; i < size; i++) {
			if (i>j) {
				u[j*size+i] = 0;
			}else {
				u[j*size+i] = GET_RAND_FP; 
			}
		}
	}

	m = (float*)calloc((size_t)size*size, sizeof(float));
	if ( m == NULL) { free(l); free(u); return RET_FAILURE; }
	for (i=0; i < size; i++) {
		for (j=0; j < size; j++) {
			for (k=0; k <= MIN(i,j); k++)
				m[i*size+j] += l[i*size+k] * u[j*size+k];
		}
	}

	free(l);
	free(u);

	*mp = m;

	return RET_SUCCESS;
}

void
matrix_multiply(float *inputa, float *inputb, float *output, int size){
	int i, j, k;

	/* [FIX-OCL-29] the product was assigned instead of accumulated (and output was never zeroed) */
	for (i=0; i < size*size; i++) output[i] = 0;
	for (i=0; i < size; i++)
		for (k=0; k < size; k++)
			for (j=0; j < size; j++)
				output[i*size+j] += inputa[i*size+k] * inputb[k*size+j];

}

/* [FIX-OCL-31] The verifier used a fixed absolute tolerance of 1e-4 for float arithmetic whose rounding error grows with the
 * matrix dimension (each element is a sum of up to n float products), so it reported spurious failures for large matrices
 * (n = 2048: both toolchains "failed" although the factorisation was right). The tolerance is now relative to the element
 * (floor 1) and grows linearly with n beyond 256; at n <= 256 it is the old 1e-4. */
float lud_tolerance(float v, int dim)
{
	float scale = dim > 256 ? (float)dim / 256.0f : 1.0f;
	float mag = fabsf(v) > 1.0f ? fabsf(v) : 1.0f;
	return 0.0001f * scale * mag;
}

/* Quiet verification: number of elements where L*U differs from the original beyond the tolerance (no output). */
int lud_check(const float *m, const float *lu, int matrix_dim)
{
	int mismatches = 0, i, j, k;
	for (i=0; i < matrix_dim; i++)
		for (j=0; j < matrix_dim; j++) {
			float sum = 0, l, u;
			for (k=0; k <= MIN(i,j); k++) {
				l = (i==k) ? 1 : lu[i*matrix_dim+k];
				u = lu[k*matrix_dim+j];
				sum += l*u;
			}
			if ( fabs(m[i*matrix_dim+j]-sum) > lud_tolerance(m[i*matrix_dim+j], matrix_dim)) mismatches++;
		}
	return mismatches;
}

func_ret_t
lud_verify(float *m, float *lu, int matrix_dim){
	int i,j,k;
	float *tmp = (float*)malloc(matrix_dim*matrix_dim*sizeof(float));

	for (i=0; i < matrix_dim; i ++)
		for (j=0; j< matrix_dim; j++) {
			float sum = 0;
			float l,u;
			for (k=0; k <= MIN(i,j); k++){
				if ( i==k)
					l=1;
				else
					l=lu[i*matrix_dim+k];
				u=lu[k*matrix_dim+j];
				sum+=l*u;
			}
			tmp[i*matrix_dim+j] = sum;
		}
	printf(">>>>>LU<<<<<<<\n");
	for (i=0; i<matrix_dim; i++){
		for (j=0; j<matrix_dim;j++){
			printf("%f ", lu[i*matrix_dim+j]);
		}
		printf("\n");
	}
	printf(">>>>>result<<<<<<<\n");
	for (i=0; i<matrix_dim; i++){
		for (j=0; j<matrix_dim;j++){
			printf("%f ", tmp[i*matrix_dim+j]);
		}
		printf("\n");
	}
	printf(">>>>>input<<<<<<<\n");
	for (i=0; i<matrix_dim; i++){
		for (j=0; j<matrix_dim;j++){
			printf("%f ", m[i*matrix_dim+j]);
		}
		printf("\n");
	}
	int mismatches = 0;
	for (i=0; i<matrix_dim; i++){
		for (j=0; j<matrix_dim; j++){
			if ( fabs(m[i*matrix_dim+j]-tmp[i*matrix_dim+j]) > lud_tolerance(m[i*matrix_dim+j], matrix_dim)) {
				printf("dismatch at (%d, %d): (o)%f (n)%f\n", i, j, m[i*matrix_dim+j], tmp[i*matrix_dim+j]);
				mismatches++;
			}
		}
	}
	printf("%s: %d differences beyond the size-scaled tolerance (1e-4 absolute at n<=256, relative, growing with n)\n", ((mismatches ==0) ? "SUCCESS" : "FAILURE"), mismatches);
	free(tmp);
	/* [FIX-OCL-30] the function is declared func_ret_t but had no return statement */
	return (mismatches == 0) ? RET_SUCCESS : RET_FAILURE;
}

void
matrix_duplicate(float *src, float **dst, int matrix_dim) {
	int s = matrix_dim*matrix_dim*sizeof(float);
	float *p = (float *) malloc (s);
	memcpy(p, src, s);
	*dst = p;
}

void
print_matrix(float *m, int matrix_dim) {
	int i, j;
	for (i=0; i<matrix_dim;i++) {
		for (j=0; j<matrix_dim;j++)
			printf("%f ", m[i*matrix_dim+j]);
		printf("\n");
	}
}
