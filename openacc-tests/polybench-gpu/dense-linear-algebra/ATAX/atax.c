/**
 * atax.c: This file is part of the PolyBench/GPU 1.0 test suite.
 *
 *
 * Contact: Scott Grauer-Gray <sgrauerg@gmail.com>
 * Louis-Noel Pouchet <pouchet@cse.ohio-state.edu>
 * Web address: http://www.cse.ohio-state.edu/~pouchet/software/polybench/GPU
 */

#include <stdio.h>
#include <unistd.h>
#include <string.h>
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
int NX=4096*4;
int NY=4096*4;

/* Constant for pi */
/* #define M_PI 3.14159 */

/* Can switch DATA_TYPE between float and double */
typedef float DATA_TYPE;




void gataxLoopb(DATA_TYPE **aa, DATA_TYPE *ya, DATA_TYPE *tmpa)
{
    int i, j;


#pragma omp parallel default(shared) private(j,i)
    {
        int tid = omp_get_thread_num();

        do{
            ctsar_next(s,NX);
            int gts = CSTART(s,tid);
            int gte = CEND(s,tid);
            int col_base = gts > -1 ? 0 : CSTART(s,tid);
            int col_width = NX;
            DATA_TYPE *a,*x,*y;

            a = ctsar_reg_mem_2d(s, aa[0], sizeof(DATA_TYPE), NX, NY,
                    CTSAR_MEM_INPUT | CTSAR_MEM_COLUMN, 0, 0, NULL);
            y = ctsar_reg_mem(s, ya, sizeof(DATA_TYPE), NX,
                    CTSAR_MEM_OUTPUT | CTSAR_MEM_PARTIAL);
            DATA_TYPE *tmp = ctsar_reg_mem(s, tmpa, sizeof(DATA_TYPE), NX,
                    CTSAR_MEM_INPUT);

            ctsar_start(s);
#pragma acc data region  deviceptr(a,x,y,tmp) if(ctsar_get_type(s) == CTSAR_DEV_GPU)
            {

#pragma acc kernels for independent private(i,j) if(ctsar_get_type(s) == CTSAR_DEV_GPU)
                for (j = gts; j < gte; j++)
                {
                    DATA_TYPE t = 0;
#pragma acc for reduction(+:t)
                    for(i = 0; i < NY; i++)
                    {
                        t = t + a[(i*col_width) + j - col_base] * tmp[i];
                    }
                    y[j] = t;
                }
            }
            ctsar_end(s);
            /* printf("here: %d\n", omp_get_thread_num()); */
        }while(ctsar_loop(s));
    }
}


void ataxLoopa(DATA_TYPE **a, DATA_TYPE *x, DATA_TYPE *y, DATA_TYPE *tmp)
{
    int i, j;


    for (i= 0; i < NX; i++)
    {
        y[i] = 0;
    }

#pragma omp parallel for private(i,j) schedule(runtime)
    for (i = 0; i < NX; i++)
    {
        DATA_TYPE t = 0;

        for (j = 0; j < NY; j++)
        {
            t = t + a[i][j] * x[j];
        }

        tmp[i] = t;
    }
}

void ataxLoopb(DATA_TYPE **a, DATA_TYPE *y, DATA_TYPE *tmp)
{
    int i, j;


#pragma omp parallel for private(i,j) schedule(runtime)
    for (j = 0; j < NX; j++)
    {
        for(i = 0; i < NY; i++)
        {
            y[j] = y[j] + a[i][j] * tmp[i];
        }
    }
}


void init_array(DATA_TYPE **A, DATA_TYPE *x)
{
    int i, j;

    for (i = 0; i < NX; i++)
    {
        x[i] = i * M_PI;
        for (j = 0; j < NY; j++)
        {
            A[i][j] = ((DATA_TYPE) i*j) / NX;
        }
    }
}


void compareResults(DATA_TYPE *y, DATA_TYPE *y_outputFromGpu)
{
    int i, fail;
    fail = 0;

    for (i=0; i<NY; i++)
    {
        if (percentDiff(y[i], y_outputFromGpu[i]) > PERCENT_DIFF_ERROR_THRESHOLD)
        {
                if(fail < 10)
                    printf("non-matching at i=%d, cpu=%f gpu=%f\n", i, y[i], y_outputFromGpu[i]);
            fail++;
        }
    }

    // print results
    printf("Non-Matching CPU-GPU Outputs Beyond Error Threshold of %4.2f Percent: %d\n", PERCENT_DIFF_ERROR_THRESHOLD, fail);

}



int main(int argc, char** argv)
{
    double t_start, t_end;
    char * env;
    if((env = getenv("TEST_SIZE")) && atoi(env) > 0){
      NX = atoi(env);
      NY = atoi(env);
    }

    /* Array declaration */
    DATA_TYPE **A = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NX,NY);
    DATA_TYPE *tmp = (DATA_TYPE *)ctsar_calloc_test(sizeof(DATA_TYPE),NX);
    DATA_TYPE *tmp_Gpu = (DATA_TYPE *)ctsar_calloc_test(sizeof(DATA_TYPE),NX);
    DATA_TYPE *x = (DATA_TYPE *)ctsar_calloc_test(sizeof(DATA_TYPE),NY);
    DATA_TYPE *y = (DATA_TYPE *)ctsar_calloc_test(sizeof(DATA_TYPE),NY);
    DATA_TYPE *y_outputFromGpu = (DATA_TYPE *)ctsar_calloc_test(sizeof(DATA_TYPE),NY);


    /* Initialize array. */
    init_array(A, x);
    ctsar_pre_init();


    t_start = rtclock();
    int iters;
    if(env = getenv("TEST_ITERATIONS")){
        iters = atoi(env);
    }else{
        iters = 1;
    }
    t_start = rtclock();

    ctsar_init(&s,NX,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
    for(int i=0; i < iters; i++)
    {
        double lt_start, lt_end;
        lt_start = rtclock();
        /* grunCorr(symmat_outputFromGpu, symmat_outputFromGpu, stddev_Gpu, mean_Gpu, float_n, eps); */
        if((env = getenv("OMP_CTSAR_FALLBACK")) && atoi(env)){
            ataxLoopa(A, x, y, tmp);
            ataxLoopb(A, y, tmp);
        }else{
            ctsar_init(&s,NX,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
            ataxLoopa(A, x, y_outputFromGpu, tmp_Gpu);
            gataxLoopb(A, y_outputFromGpu, tmp_Gpu);
        }
        lt_end = rtclock();
        fprintf(stderr, "GPU Runtime:it%d: %0.6lfs\n", i, lt_end - lt_start);

    }


    t_end = rtclock();
    fprintf(stderr, "GPU Runtime: %0.6lfs\n", t_end - t_start);

    if((env = getenv("TEST_VERIFY")) && atoi(env)){
        t_start = rtclock();

        ataxLoopa(A, x, y, tmp);
        ataxLoopb(A, y, tmp);

        t_end = rtclock();
        fprintf(stderr, "CPU Runtime: %0.6lfs\n", t_end - t_start);

        compareResults(y, y_outputFromGpu);
    }

    return 0;
}
