/**
 * twomm.c: This file is part of the PolyBench/GPU 1.0 test suite.
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
ctsar * ct = NULL;

//define the error threshold for the results "not matching"
#define PERCENT_DIFF_ERROR_THRESHOLD 0.05

/* Problem size */
/* #define NI 512
#define NJ 512
#define NK 512
#define NL 512 */

#define NI 2046
#define NJ 2046
#define NK 2046
#define NL 2046

/* Can switch DATA_TYPE between float and double */
typedef float DATA_TYPE;

void gtwoMMloopa(DATA_TYPE **aa, DATA_TYPE **ba, DATA_TYPE **ca)
{
    /* Determine mean of column vectors of input data matrix */
#pragma omp parallel default(shared)
    {
    int i, j, k;
        int tid = omp_get_thread_num();
        /* if(s == NULL)
            perror("what the hell?"); */

        do{
            ctsar_next(s,NI);
            int gts = CSTART(s,tid);
            int gte = CEND(s,tid);
            int ni = NI;
            DATA_TYPE * a, *b, *c;
            a = ctsar_reg_mem(s, aa[0], sizeof(DATA_TYPE)*(NI), (NJ),
                    CTSAR_MEM_PARTIAL | CTSAR_MEM_INPUT);
            b = ctsar_reg_mem(s, ba[0], sizeof(DATA_TYPE)*(NI), (NJ),
                    CTSAR_MEM_INPUT);
            c = ctsar_reg_mem(s, ca[0], sizeof(DATA_TYPE)*(NI), (NJ),
                    CTSAR_MEM_PARTIAL | CTSAR_MEM_INPUT | CTSAR_MEM_OUTPUT);

            ctsar_start(s);
#pragma acc data region deviceptr(a,b,c)  if(ctsar_get_type(s) == CTSAR_DEV_GPU)
            {
#pragma acc region for independent private(i,j,k) if(ctsar_get_type(s) == CTSAR_DEV_GPU)
                for(i = gts; i < gte; i++)
                {
#pragma acc for independent
                    for(j = 0; j < NJ; j++)
                    {
                        DATA_TYPE lc = c[((i) * NI) + (j)];
/* #pragma acc for seq */
                        for (k = 0; k < NK; k++)
                        {
                            lc += a[((i) * NI) + (k)] * b[((k) * NI) + (j)];
                        }
                        c[((i) * NI) + (j)] = lc;
                    }
                } 
            }
            ctsar_end(s);
            /* printf("here: %d\n", omp_get_thread_num()); */
        }while(ctsar_loop(s));
    }
}

void gtwoMMloopb(DATA_TYPE **ca, DATA_TYPE **da, DATA_TYPE **ea)
{
    /* Determine mean of column vectors of input data matrix */
#pragma omp parallel default(shared)
    {
    int i, j, k;
        int tid = omp_get_thread_num();
        /* if(s == NULL)
            perror("what the hell?"); */

        do{
            ctsar_next(ct,NI);
            int gts = CSTART(ct,tid);
            int gte = CEND(ct,tid);
            int ni = NI;
            DATA_TYPE * d, *e, *cg;
            cg = ctsar_reg_mem(ct, ca[0], sizeof(DATA_TYPE)*(NI), (NJ),
                    CTSAR_MEM_PARTIAL | CTSAR_MEM_INPUT);
            d = ctsar_reg_mem(ct, da[0], sizeof(DATA_TYPE)*(NI), (NJ),
                    CTSAR_MEM_INPUT);
            e = ctsar_reg_mem(ct, ea[0], sizeof(DATA_TYPE)*(NI), (NJ),
                    CTSAR_MEM_PARTIAL | CTSAR_MEM_INPUT | CTSAR_MEM_OUTPUT);

            ctsar_start(ct);
#pragma acc data region deviceptr(e,d,cg)  if(ctsar_get_type(ct) == CTSAR_DEV_GPU)
            {
#pragma acc region for independent private(i,j,k)  if(ctsar_get_type(ct) == CTSAR_DEV_GPU)
                for(i = gts; i < gte; i++) 
                {
#pragma acc for independent
                    for (j = 0; j < NL; j++)
                    {
                        DATA_TYPE le = e[((i) * NI) + (j)];
/* #pragma acc for seq */
                        for (k = 0; k < NJ; ++k)
                        {
                            le += cg[((i) * NI) + (k)] * d[((k) * NI) + (j)]; 
                        }
                        e[((i) * NI) + (j)] = le;
                    }
                }
            }
            ctsar_end(ct);
            /* printf("here: %d\n", omp_get_thread_num()); */
        }while(ctsar_loop(ct));
    }
}

void twoMMloopa(DATA_TYPE **a, DATA_TYPE **b, DATA_TYPE **c)
{
    int i, j, k;


#pragma omp parallel for private(i,j,k) schedule(runtime)
    for(i = 0; i < NI; i++)
    {
        for(j = 0; j < NJ; j++)
        {
            for (k = 0; k < NK; k++)
            {
                c[i][j] += a[i][k] * b[k][j];
            }
        }
    } 
}

void twoMMloopb(DATA_TYPE **c, DATA_TYPE **d, DATA_TYPE **e)
{
    int i, j, k;


#pragma omp parallel for private(i,j,k) schedule(runtime)
    for(i = 0; i < NI; i++) 
    {
        for (j = 0; j < NL; j++)
        {
            for (k = 0; k < NJ; ++k)
            {
                e[i][j] += c[i][k] * d[k][j]; 
            }
        }
    }
}

void init_array(DATA_TYPE **A, DATA_TYPE **B, DATA_TYPE **C, DATA_TYPE **C_gpu, DATA_TYPE **D,
        DATA_TYPE **E, DATA_TYPE **E_outputFromGpu)
{
    int i, j;

#pragma omp parallel for private(i,j)
    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NK; j++)
        {
            A[i][j] = ((DATA_TYPE) i*j)/NI;
        }
    }

    for (i = 0; i < NK; i++)
    {
        for (j = 0; j < NJ; j++)
        {
            B[i][j] = ((DATA_TYPE) i*j + 1)/NJ;
        }
    }

    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NJ; j++)
        {
            C[i][j] = ((DATA_TYPE) i*j + 2)/NJ;
            C_gpu[i][j] = ((DATA_TYPE) i*j + 2)/NJ;
        }
    }

    for (i = 0; i < NJ; i++)
    {
        for (j = 0; j < NL; j++)
        {
            D[i][j] = ((DATA_TYPE) i*j + 2)/NJ;
        }
    }

    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NL; j++)
        {
            E[i][j] = ((DATA_TYPE) i*j + 2)/NJ;
            E_outputFromGpu[i][j] = ((DATA_TYPE) i*j + 2)/NJ;
        }
    }
}



void compareResults(DATA_TYPE **E, DATA_TYPE **E_outputFromGpu)
{
    int i,j,fail;
    fail = 0;

    for (i=0; i < NI; i++)
    {
        for (j=0; j < NL; j++)
        {
            if (percentDiff(E[i][j], E_outputFromGpu[i][j]) > PERCENT_DIFF_ERROR_THRESHOLD)
            {
                fail++;
                if(fail < 20)
                    printf("i: %d j: %d ec: %f eg: %f\n", i, j, E[i][j], E_outputFromGpu[i][j]);
            }
        }
    }

    // print results
    printf("Non-Matching CPU-GPU Outputs Beyond Error Threshold of %4.2f Percent: %d\n", PERCENT_DIFF_ERROR_THRESHOLD, fail);

}


int main(int argc, char** argv)
{
    double t_start, t_end;
    int i;

    /* Array declaration */
    DATA_TYPE **A = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NI,NK);
    DATA_TYPE **B = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NK,NJ);
    DATA_TYPE **C = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NI,NJ);
    DATA_TYPE **C_gpu = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NI,NJ);
    DATA_TYPE **D = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NJ,NL);
    DATA_TYPE **E = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NI,NL);
    DATA_TYPE **E_outputFromGpu = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NI,NL);

    /* Initialize array. */
    init_array(A, B, C, C_gpu, D, E, E_outputFromGpu);

    ctsar_pre_init();

    int iters;
    char * env;
    if(env = getenv("TEST_ITERATIONS")){
        iters = atoi(env);
    }else{
        iters = 1;
    }

    t_start = rtclock();
    ctsar_init(&s,NI,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
    ctsar_init(&ct,NI,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
    for(i=0; i < iters; i++)
    {
        double lt_start, lt_end;
        lt_start = rtclock();
        /* grunCorr(symmat_outputFromGpu, symmat_outputFromGpu, stddev_Gpu, mean_Gpu, float_n, eps); */
        if((env = getenv("OMP_CTSAR_FALLBACK")) && atoi(env)){
            twoMMloopa(A, B, C_gpu);
            twoMMloopb(C_gpu, D, E_outputFromGpu);
        }else{
            ctsar_init(&s,NI,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
            gtwoMMloopa(A, B, C_gpu);

            ctsar_init(&ct,NI,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
            gtwoMMloopb(C_gpu, D, E_outputFromGpu);
        }
        lt_end = rtclock();
        fprintf(stderr, "GPU Runtime:it%d: %0.6lfs\n", i, lt_end - lt_start);
    }

    t_end = rtclock();
    fprintf(stderr, "GPU Runtime: %0.6lfs\n", t_end - t_start);



    if((env = getenv("TEST_VERIFY")) && atoi(env)){
        t_start = rtclock();

        twoMMloopa(A, B, C);
        twoMMloopb(C, D, E);

        t_end = rtclock();
        fprintf(stderr, "CPU Runtime: %0.6lfs\n", t_end - t_start);

        compareResults(E, E_outputFromGpu);

    }
    return 0;
}
