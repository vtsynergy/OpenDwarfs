/**
 * twodconv.c: This file is part of the PolyBench/GPU 1.0 test suite.
 *
 *
 * Contact: Scott Grauer-Gray <sgrauerg@gmail.com>
 * Louis-Noel Pouchet <pouchet@cse.ohio-state.edu>
 * Web address: http://www.cse.ohio-state.edu/~pouchet/software/polybench/GPU
 */

#include <unistd.h>
#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <accel.h>

#include "polybenchUtilFuncts.h"
#include "ctsar.h"

//define the error threshold for the results "not matching"
#define PERCENT_DIFF_ERROR_THRESHOLD 0.05

/* Problem dimensions */
/* #define NI 4096
#define NJ 4096 */

#define NI 8192 * 3
#define NJ 8192 * 3

/* Can switch DATA_TYPE between float and double */
typedef float DATA_TYPE;



void conv2D(DATA_TYPE **A, DATA_TYPE **B)
{
	int i, j, k;
	DATA_TYPE c11, c12, c13, c21, c22, c23, c31, c32, c33;

	c11 = +2;  c21 = +5;  c31 = -8;
	c12 = -3;  c22 = +6;  c32 = -9;
	c13 = +4;  c23 = +7;  c33 = +10;

#pragma omp parallel for private(i,j) default(shared)
	for (i = 1; i < NI - 1; ++i) // 0
	{
		for (j = 1; j < NJ - 1; ++j) // 1
		{
			B[i][j] = c11 * A[i - 1][j - 1]  +  c12 * A[i + 0][j - 1]  +  c13 * A[i + 1][j - 1]
				+ c21 * A[i - 1][j + 0]  +  c22 * A[i + 0][j + 0]  +  c23 * A[i + 1][j + 0] 
				+ c31 * A[i - 1][j + 1]  +  c32 * A[i + 0][j + 1]  +  c33 * A[i + 1][j + 1];
		}
	}
}

void conv2Dg(ctsar *s, DATA_TYPE **Aa, DATA_TYPE **Ba)
{
    int i, j, k;
    DATA_TYPE c11, c12, c13, c21, c22, c23, c31, c32, c33;

    c11 = +2;  c21 = +5;  c31 = -8;
    c12 = -3;  c22 = +6;  c32 = -9;
    c13 = +4;  c23 = +7;  c33 = +10;
#pragma omp parallel default(shared) private(i,j)
    {
        int tid = omp_get_thread_num();
        do{
            ctsar_next(s,NI);
            int gts = CSTART(s,tid) ? CSTART(s,tid) : 1;
            int gte = CEND(s,tid) == NI ? NI-1 : CEND(s,tid);
            /* {
                ctsar_start(s);
                for (i = gts; i < gte; ++i) // 0
                {
                    for (j = 1; j < NJ - 1; ++j) // 1
                    {
                        B[i][j] = c11 * A[i - 1][j - 1]  +  c12 * A[i + 0][j - 1]  +  c13 * A[i + 1][j - 1]
                            + c21 * A[i - 1][j + 0]  +  c22 * A[i + 0][j + 0]  +  c23 * A[i + 1][j + 0] 
                            + c31 * A[i - 1][j + 1]  +  c32 * A[i + 0][j + 1]  +  c33 * A[i + 1][j + 1];
                    }
                }
                ctsar_end(s);
            }else */
            {
                float *A = ctsar_reg_mem_2d(s, Aa[0], sizeof(float),(NI), (NJ),
                        CTSAR_MEM_PARTIAL | CTSAR_MEM_INPUT,0,1,NULL);
                float *B = ctsar_reg_mem(s, Ba[0], sizeof(float)*(NI), (NJ),
                        CTSAR_MEM_PARTIAL | CTSAR_MEM_INPUT | CTSAR_MEM_OUTPUT);

                ctsar_start(s);
#pragma acc data region deviceptr(A,B) if(ctsar_get_type(s)==CTSAR_DEV_GPU)
                {
#pragma acc region for independent private(i,j) if(ctsar_get_type(s)==CTSAR_DEV_GPU)
                    for (i = gts; i < gte; ++i) // 0
                    {
#pragma acc for independent
                        for (j = 1; j < NJ - 1; ++j) // 1
                        {
			B[((i) * NI) + (j)] = c11 * A[((i - 1) * NI) + (j - 1)]  +  c12 * A[((i + 0) * NI) + (j - 1)]  +  c13 * A[((i + 1) * NI) + (j - 1)]
				+ c21 * A[((i - 1) * NI) + (j + 0)]  +  c22 * A[((i + 0) * NI) + (j + 0)]  +  c23 * A[((i + 1) * NI) + (j + 0)] 
				+ c31 * A[((i - 1) * NI) + (j + 1)]  +  c32 * A[((i + 0) * NI) + (j + 1)]  +  c33 * A[((i + 1) * NI) + (j + 1)];
/* #define B(X,Y) lB[(X)*NJ+(Y)]
#define A(X,Y) lA[(X)*NJ+(Y)]
                            B(i,j) = c11 * A(i - 1,j - 1)  +  c12 * A(i + 0,j - 1)  +  c13 * A(i + 1,j - 1)
                                + c21 * A(i - 1,j + 0)  +  c22 * A(i + 0,j + 0)  +  c23 * A(i + 1,j + 0) 
                                + c31 * A(i - 1,j + 1)  +  c32 * A(i + 0,j + 1)  +  c33 * A(i + 1,j + 1);
#undef B
#undef A */
                        }
                    }
                }
                ctsar_end(s);
            }
            /* printf("here: %d\n", omp_get_thread_num()); */
        }while(ctsar_loop(s));
    }
}

void init(DATA_TYPE **A)
{
	int i, j;

	for (i = 0; i < NI; ++i)
	{
		for (j = 0; j < NJ; ++j)
		{
			A[i][j] = i % 12 + 2 * (j % 7);
		}
	}
}


void compareResults(DATA_TYPE **B, DATA_TYPE **B_outputFromGpu)
{
	int i, j, fail;
	fail = 0;
	
	// Compare a and b
	for (i=1; i < (NI-1); i++) 
	{
		for (j=1; j < (NJ-1); j++) 
		{
			if (percentDiff(B[i][j], B_outputFromGpu[i][j]) > PERCENT_DIFF_ERROR_THRESHOLD) 
			{
				fail++;
			}
		}
	}
	
	// Print results
	printf("Non-Matching CPU-GPU Outputs Beyond Error Threshold of %4.2f Percent: %d\n", PERCENT_DIFF_ERROR_THRESHOLD, fail);
	
}


int main(int argc, char *argv[])
{
    double t_start, t_end;

    DATA_TYPE **A = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NI,NJ);
    DATA_TYPE **B = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NI,NJ);  // CPU target results
    DATA_TYPE **B_outputFromGpu = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NI,NJ);  // GPU exec results

    //test to confirm the memory layout, allocated as contiguous buffer
    ctsar_pre_init();
    init(A);

    int iters,i;
    char * env;
    if(env = getenv("TEST_ITERATIONS")){
        iters = atoi(env);
    }else{
        iters = 1;
    }
    t_start = rtclock();
    ctsar * s = NULL;
    ctsar_init(&s,NI,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
    for(i=0; i < iters; i++)
    {
        double lt_start, lt_end;
        lt_start = rtclock();
        ctsar_init(&s,NI,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
        conv2Dg(s,A, B_outputFromGpu);
        lt_end = rtclock();
        fprintf(stderr, "GPU Runtime:it%d: %0.6lfs\n", i, lt_end - lt_start);
    }
    t_end = rtclock();
    fprintf(stderr, "GPU Runtime: %0.6lf\n", t_end - t_start);

    t_start = rtclock();
    conv2D(A, B);
    t_end = rtclock();
    fprintf(stderr, "CPU Runtime: %0.6lf\n", t_end - t_start);

    compareResults(B, B_outputFromGpu);

    return 0;
}
