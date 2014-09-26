/**
 * syrk.c: This file is part of the PolyBench/GPU 1.0 test suite.
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

/* Problem size */
#define N 2048
#define M 2048

/* Can switch DATA_TYPE between float and double */
typedef float DATA_TYPE;



void grunSyrTwoK(DATA_TYPE **a_a, DATA_TYPE **c_a)
{
    int i, j, k, n, m;
    DATA_TYPE alpha, beta;

    alpha = 12435;
    beta = 4546;

    /*    C := alpha*A*B' + alpha*B*A' + beta*C */

#pragma omp parallel default(shared) private(j,i,k)
    {
        int tid = omp_get_thread_num();

#pragma omp for private(i,j) schedule(runtime)
        for (i = 0; i < N; ++i) // 0
        {
            for (j = 0; j < N; j++)
            {
                c_a[0][(i * N) + j] *= beta;
            }
        }
        do{
            ctsar_next(s,N);
            int gts = CSTART(s,tid);
            int gte = CEND(s,tid);
            DATA_TYPE *a,*c;
            a = ctsar_reg_mem(s, *a_a, sizeof(DATA_TYPE)*N, N,
                    CTSAR_MEM_INPUT);
            c = ctsar_reg_mem(s, *c_a, sizeof(DATA_TYPE)*N, N,
                    CTSAR_MEM_PARTIAL | CTSAR_MEM_INPUT | CTSAR_MEM_OUTPUT);

            ctsar_start(s);
#pragma acc data region deviceptr(a,c)  if(ctsar_get_type(s) == CTSAR_DEV_GPU)
            {

#pragma acc region for independent if(ctsar_get_type(s) == CTSAR_DEV_GPU) \
                private(i,j,k)
                for (i = gts; i < gte; ++i) // 0
                {
                    DATA_TYPE * ia = &a[i*N];
#pragma acc for independent
                    for (j = 0; j < N; j++)
                    {
                        DATA_TYPE l1 = c[(i*N) + j];
                        DATA_TYPE * ja = &a[j*N];
#pragma acc for reduction(+: l1)
                        for (k = 0; k < M; k++)
                        {
                            l1 += alpha * ia[k] * ja[k];
                        }
                        c[(i*N) + j] = l1;
                    }
                }
            }
            ctsar_end(s);
            /* printf("here: %d\n", omp_get_thread_num()); */
        }while(ctsar_loop(s));
    }
}

void runSyrTwoK(DATA_TYPE **a, DATA_TYPE **c)
{
    int i, j, k, n, m;
    DATA_TYPE alpha, beta;

    alpha = 12435;
    beta = 4546;

    /*    C := alpha*A*B' + alpha*B*A' + beta*C */

#pragma omp parallel 
    {
#pragma omp for private(i,j) schedule(runtime)
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                c[i][j] *= beta;
            }
        }

#pragma omp for private(i,j,k) schedule(runtime)
        for (i = 0; i < N; i++)
        {
            for (j = 0; j < N; j++)
            {
                for (k = 0; k < M; k++)
                {
                    c[i][j] += alpha * a[i][k] * a[j][k];
                }
            }
        }
    }
}


void init_arrays(DATA_TYPE **A, DATA_TYPE **C, DATA_TYPE **C_Gpu)
{
    int i, j;

    for (i = 0; i < N; i++)
    {
        for (j = 0; j < N; j++)
        {
            C[i][j] = ((DATA_TYPE) i*j + 2) / N;
            C_Gpu[i][j] = ((DATA_TYPE) i*j + 2) / N;
        }

        for (j = 0; j < M; j++)
        {
            A[i][j] = ((DATA_TYPE) i*j) / N;
        }
    }
}


void compareResults(DATA_TYPE **C, DATA_TYPE **C_outputFromGpu)
{
    int i,j,fail;
    fail = 0;

    // Compare C with D
    for (i=0; i<N; i++)
    {
        for (j=0; j<N; j++)
        {
            if (percentDiff(C[i][j], C_outputFromGpu[i][j]) > PERCENT_DIFF_ERROR_THRESHOLD)
            { 
                fail++;
            }
        }
    }

    // print results
    printf("Non-Matching CPU-GPU Outputs Beyond Error Threshold of %4.2f Percent: %d\n", PERCENT_DIFF_ERROR_THRESHOLD, fail);

}


int main(int argc, char** argv)
{
    double t_start, t_end;

    /* Array declaration */
    DATA_TYPE **A = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),N,M);
    DATA_TYPE **C = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),N,N);
    DATA_TYPE **C_outputFromGpu = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),N,N);

    /* Initialize array. */
    init_arrays(A, C, C_outputFromGpu);

    ctsar_pre_init();

    // Run GPU code

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
            runSyrTwoK(A, C_outputFromGpu);
        }else{
            ctsar_init(&s,N,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
            grunSyrTwoK(A, C_outputFromGpu);
        }
        lt_end = rtclock();
        fprintf(stderr, "GPU Runtime:it%d: %0.6lfs\n", i, lt_end - lt_start);

    }


    t_end = rtclock();
    fprintf(stderr, "GPU Runtime: %0.6lfs\n", t_end - t_start);


    if((env = getenv("TEST_VERIFY")) && atoi(env)){
        t_start = rtclock();
        runSyrTwoK(A, C);
        t_end = rtclock();
        fprintf(stderr, "CPU Runtime: %0.6lfs\n", t_end - t_start);

        compareResults(C, C_outputFromGpu);
    }

    return 0;
}
