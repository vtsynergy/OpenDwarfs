/**
 * gesummv.c: This file is part of the PolyBench/GPU 1.0 test suite.
 *
 *
 * Contact: Scott Grauer-Gray <sgrauerg@gmail.com>
 * Louis-Noel Pouchet <pouchet@cse.ohio-state.edu>
 * Web address: http://www.cse.ohio-state.edu/~pouchet/software/polybench/GPU
 */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdarg.h>

#include "polybenchUtilFuncts.h"
#include "ctsar.h"
ctsar * s = NULL;

//define the error threshold for the results "not matching"
#define PERCENT_DIFF_ERROR_THRESHOLD 0.05

/* Problem size. */
#define N 4096 * 5
/* #define N 16384 */

/* Can switch DATA_TYPE between float and double */
typedef float DATA_TYPE;


void runGesummv_gpu(DATA_TYPE **ai, DATA_TYPE **bi, DATA_TYPE *x1i, DATA_TYPE *y1i, DATA_TYPE *tmp1i)
{
    int i, j;

    DATA_TYPE alpha = 43532;
    DATA_TYPE beta = 12313;

#pragma omp parallel default(shared) private(j,i)
    {
        int tid = omp_get_thread_num();

        do{
            ctsar_next(s,N);
            int gts = CSTART(s,tid);
            int gte = CEND(s,tid);
            DATA_TYPE *a, *b, *x1, *y1;
            a = ctsar_reg_mem(s, *ai, sizeof(DATA_TYPE)*N, N,
                    CTSAR_MEM_PARTIAL | CTSAR_MEM_INPUT);
            b = ctsar_reg_mem(s, *bi, sizeof(DATA_TYPE)*N, N,
                    CTSAR_MEM_PARTIAL | CTSAR_MEM_INPUT);
            x1 = ctsar_reg_mem(s, x1i, sizeof(DATA_TYPE), N,
                    CTSAR_MEM_INPUT);
            y1 = ctsar_reg_mem(s, y1i, sizeof(DATA_TYPE), N,
                    CTSAR_MEM_PARTIAL | CTSAR_MEM_OUTPUT);
            /* tmp1 = ctsar_reg_mem(s, tmp1i, sizeof(DATA_TYPE), N,
                    CTSAR_MEM_NONE); */

            ctsar_start(s);
#pragma acc data region deviceptr(a,b,x1,y1)  if(ctsar_get_type(s) == CTSAR_DEV_GPU)
            {
#pragma acc region for private(i,j) independent  if(ctsar_get_type(s) == CTSAR_DEV_GPU)
                for (i = gts; i < gte; i++)
                {
                    DATA_TYPE tmp1 = 0.0;
                    DATA_TYPE tmp2 = 0.0;

#pragma acc for seq
                    for (j = 0; j < N; j++)
                    {
                        tmp1 += a[(i*N) + j] * x1[j];
                        tmp2 += b[(i*N) + j] * x1[j];
                    }
                    y1[i] = alpha * tmp1 + beta * tmp2;
                }
            }
            ctsar_end(s);
            /* printf("here: %d\n", omp_get_thread_num()); */
        }while(ctsar_loop(s));
    }
}

void runGesummv(DATA_TYPE **a, DATA_TYPE **b, DATA_TYPE *x1, DATA_TYPE *y1, DATA_TYPE *tmp1)
{
    int i, j;

    DATA_TYPE alpha = 43532;
    DATA_TYPE beta = 12313;


#pragma omp parallel for private(i,j) schedule(runtime)
    for (i = 0; i < N; i++)
    {
        tmp1[i] = 0;
        y1[i] = 0;

        for (j = 0; j < N; j++)
        {
            tmp1[i] = a[i][j] * x1[j] + tmp1[i];
            y1[i] = b[i][j] * x1[j] + y1[i];
        }
        y1[i] = alpha * tmp1[i] + beta * y1[i];
    }
}


void init(DATA_TYPE **A, DATA_TYPE *x)
{
    int i, j;

    for (i = 0; i < N; i++)
    {
        x[i] = ((DATA_TYPE) i) / N;

        for (j = 0; j < N; j++) 
        {
            A[i][j] = ((DATA_TYPE) i*j) / N;
        }
    }
}


void compareResults(DATA_TYPE *y, DATA_TYPE *y_outputFromGpu)
{
    int i, fail;
    fail = 0;

    for (i=0; i< N; i++) 
    {
        if (percentDiff(y[i], y_outputFromGpu[i]) > PERCENT_DIFF_ERROR_THRESHOLD) 
        {
            fail++;
        }
    }

    // Print results
    printf("Non-Matching CPU-GPU Outputs Beyond Error Threshold of %4.2f Percent: %d\n", PERCENT_DIFF_ERROR_THRESHOLD, fail);

}


int main(int argc, char** argv)
{
    double t_start, t_end;

    /* Array declaration */
    DATA_TYPE **A = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE), N, N);
    DATA_TYPE **B = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE), N, N);
    DATA_TYPE *x = (DATA_TYPE *)ctsar_calloc_test(sizeof(DATA_TYPE), N);
    DATA_TYPE *y = (DATA_TYPE *)ctsar_calloc_test(sizeof(DATA_TYPE), N);
    DATA_TYPE *y_outputFromGpu = (DATA_TYPE *)ctsar_calloc_test(sizeof(DATA_TYPE), N);
    DATA_TYPE *tmp = (DATA_TYPE *)ctsar_calloc_test(sizeof(DATA_TYPE), N);

    /* Initialize array. */
    init(A, x);
    ctsar_pre_init();

    int iters;
    char * env;
    if(env = getenv("TEST_ITERATIONS")){
        iters = atoi(env);
    }else{
        iters = 1;
    }
    t_start = rtclock();

    ctsar_init(&s,N,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
    for(int i=0; i < iters; i++)
    {
        double lt_start, lt_end;
        lt_start = rtclock();
        /* grunCorr(symmat_outputFromGpu, symmat_outputFromGpu, stddev_Gpu, mean_Gpu, float_n, eps); */
        if((env = getenv("OMP_CTSAR_FALLBACK")) && atoi(env)){
            runGesummv(A, B, x, y_outputFromGpu, tmp);
        }else{
            ctsar_init(&s,N,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
            runGesummv_gpu(A, B, x, y_outputFromGpu, tmp);
        }
        lt_end = rtclock();
        fprintf(stderr, "GPU Runtime:it%d: %0.6lfs\n", i, lt_end - lt_start);
    }


    t_end = rtclock();
    fprintf(stderr, "GPU Runtime: %0.6lfs\n", t_end - t_start);

    if((env = getenv("TEST_VERIFY")) && atoi(env)){
        t_start = rtclock();

        runGesummv(A, B, x, y, tmp);

        t_end = rtclock();
        fprintf(stderr, "CPU Runtime: %0.6lfs\n", t_end - t_start);

        compareResults(y, y_outputFromGpu);
    }

    return 0;
}
