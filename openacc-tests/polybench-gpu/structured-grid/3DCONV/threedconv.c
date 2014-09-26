/**
 * threedconv.c: This file is part of the PolyBench/GPU 1.0 test suite.
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
#include <errno.h>

#include "polybenchUtilFuncts.h"
#include "ctsar.h"
ctsar * s = NULL;

//define the error threshold for the results "not matching"
#define PERCENT_DIFF_ERROR_THRESHOLD 0.05


/* Problem dimensions */
#define NI 768 
#define NJ 768 
#define NK 768 

/* Can switch DATA_TYPE between float and double */
typedef float DATA_TYPE;



void conv3D(DATA_TYPE ***A, DATA_TYPE ***B)
{
    int i, j, k;
    DATA_TYPE c11, c12, c13, c21, c22, c23, c31, c32, c33;

    c11 = +2;  c21 = +5;  c31 = -8;
    c12 = -3;  c22 = +6;  c32 = -9;
    c13 = +4;  c23 = +7;  c33 = +10;

#pragma omp parallel for private(i,j,k) schedule(runtime)
    for (i = 1; i < NI - 1; ++i) // 0
    {
        for (j = 1; j < NJ - 1; ++j) // 1
        {
            for (k = 1; k < NK - 1; ++k) // 2
            {
                B[i][j][k] = 0
                    +   c11 * A[i - 1][j - 1][k - 1]  +  c13 * A[i + 1][j - 1][k - 1]
                    +   c21 * A[i - 1][j - 1][k - 1]  +  c23 * A[i + 1][j - 1][k - 1]
                    +   c31 * A[i - 1][j - 1][k - 1]  +  c33 * A[i + 1][j - 1][k - 1]
                    +   c12 * A[i + 0][j - 1][k + 0]
                    +   c22 * A[i + 0][j + 0][k + 0]
                    +   c32 * A[i + 0][j + 1][k + 0]
                    +   c11 * A[i - 1][j - 1][k + 1]  +  c13 * A[i + 1][j - 1][k + 1]
                    +   c21 * A[i - 1][j + 0][k + 1]  +  c23 * A[i + 1][j + 0][k + 1]
                    +   c31 * A[i - 1][j + 1][k + 1]  +  c33 * A[i + 1][j + 1][k + 1];
            }
        }
    }
}

void conv3D_gpu(DATA_TYPE *a_base, DATA_TYPE *b_base)
{
    int i, j, k;
    DATA_TYPE c11, c12, c13, c21, c22, c23, c31, c32, c33;

    c11 = +2;  c21 = +5;  c31 = -8;
    c12 = -3;  c22 = +6;  c32 = -9;
    c13 = +4;  c23 = +7;  c33 = +10;


#pragma omp parallel default(shared) private(j,i,k)
    {
        int tid = omp_get_thread_num();

        do{
            ctsar_next(s,NI);
            int gts = CSTART(s,tid) ? CSTART(s,tid) : 1;
            int gte = CEND(s,tid) == NI ? NI-1 : CEND(s,tid);
            DATA_TYPE *la,*lb;
            la = ctsar_reg_mem_2d(s, a_base, sizeof(DATA_TYPE), NI, NJ*NK,
                    CTSAR_MEM_INPUT | CTSAR_MEM_PARTIAL, 0, 1,NULL);
            lb = ctsar_reg_mem(s, b_base, sizeof(DATA_TYPE)*NJ*NK, NI,
                    CTSAR_MEM_PARTIAL | CTSAR_MEM_OUTPUT);

            ctsar_start(s);
#pragma acc data region deviceptr(la,lb)  if(ctsar_get_type(s) == CTSAR_DEV_GPU)
            {
#pragma acc region for independent if(ctsar_get_type(s) == CTSAR_DEV_GPU) \
                private(i,j,k)
                // copyout(b_base[0:NI*NJ*NK])
                for (i = gts; i < gte; ++i) // 0
                {
                    DATA_TYPE * ilb = lb+NJ*NK*i;
#pragma acc for gang independent
                    for (j = 1; j < NJ - 1; ++j) // 1
                    {
                        DATA_TYPE * jlb = ilb+NK*j;
#pragma acc for independent
                        for (k = 1; k < NK - 1; ++k) // 2
                        {
                            jlb[k] = 0
                                +   c11 * la[((NJ*NK) * (i - 1)) + (NK * (j - 1)) + (k - 1)]  +  c13 * la[((NJ*NK) * (i + 1)) + (NK * (j - 1)) + (k - 1)]
                                +   c21 * la[((NJ*NK) * (i - 1)) + (NK * (j - 1)) + (k - 1)]  +  c23 * la[((NJ*NK) * (i + 1)) + (NK * (j - 1)) + (k - 1)]
                                +   c31 * la[((NJ*NK) * (i - 1)) + (NK * (j - 1)) + (k - 1)]  +  c33 * la[((NJ*NK) * (i + 1)) + (NK * (j - 1)) + (k - 1)]
                                +   c12 * la[((NJ*NK) * (i + 0)) + (NK * (j - 1)) + (k + 0)]
                                +   c22 * la[((NJ*NK) * (i + 0)) + (NK * (j + 0)) + (k + 0)]
                                +   c32 * la[((NJ*NK) * (i + 0)) + (NK * (j + 1)) + (k + 0)]
                                +   c11 * la[((NJ*NK) * (i - 1)) + (NK * (j - 1)) + (k + 1)]  +  c13 * la[((NJ*NK) * (i + 1)) + (NK * (j - 1)) + (k + 1)]
                                +   c21 * la[((NJ*NK) * (i - 1)) + (NK * (j + 0)) + (k + 1)]  +  c23 * la[((NJ*NK) * (i + 1)) + (NK * (j + 0)) + (k + 1)]
                                +   c31 * la[((NJ*NK) * (i - 1)) + (NK * (j + 1)) + (k + 1)]  +  c33 * la[((NJ*NK) * (i + 1)) + (NK * (j + 1)) + (k + 1)];
                        }
                    }
                }
            }
            ctsar_end(s);
            /* printf("here: %d\n", omp_get_thread_num()); */
        }while(ctsar_loop(s));
    }
}

void init(DATA_TYPE ***A)
{
    int i, j, k;

    for (i = 0; i < NI; ++i)
    {
        for (j = 0; j < NJ; ++j)
        {
            for (k = 0; k < NK; ++k)
            {
                A[i][j][k] = i % 12 + 2 * (j % 7) + 3 * (k % 13);
            }
        }
    }
}


void compareResults(DATA_TYPE ***B, DATA_TYPE ***B_outputFromGpu)
{
    int i, j, k, fail;
    fail = 0;

    // Compare result from cpu and gpu...
#pragma omp parallel for reduction(+:fail)
    for (i = 1; i < NI - 1; ++i) // 0
    {
        for (j = 1; j < NJ - 1; ++j) // 1
        {
            for (k = 1; k < NK - 1; ++k) // 2
            {
                if (percentDiff(B[i][j][k], B_outputFromGpu[i][j][k]) > PERCENT_DIFF_ERROR_THRESHOLD)
                {
                    fail++;
                }
            }
        }
    }

    // Print results
    printf("Non-Matching CPU-GPU Outputs Beyond Error Threshold of %4.2f Percent: %d\n", PERCENT_DIFF_ERROR_THRESHOLD, fail);

}

void * calloc_test(size_t size, size_t count){
    static int cnt=0;
    void * ptr = NULL;
    ptr = calloc(size, count);
    if(ptr == NULL){
        perror("calloc failed");
        exit(errno);
    }
    return ptr;
}

int main(int argc, char *argv[])
{
    double t_start, t_end;

    /* DATA_TYPE A[NI][NJ][NK]; */
    DATA_TYPE ***A = (DATA_TYPE ***)ctsar_alloc_3d(sizeof(DATA_TYPE), NI, NJ, NK);
    DATA_TYPE *a_base = **A;

    DATA_TYPE ***B = (DATA_TYPE ***)ctsar_alloc_3d(sizeof(DATA_TYPE), NI, NJ, NK);
    DATA_TYPE *b_base = **B;

    DATA_TYPE ***B_outputFromGpu = (DATA_TYPE ***)ctsar_alloc_3d(sizeof(DATA_TYPE), NI, NJ, NK);  // GPU exec results
    DATA_TYPE *bg_base = **B_outputFromGpu;

    //initialize the arrays
    init(A);

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
    ctsar_init(&s,NI,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);

    for(int i=0; i < iters; i++)
    {
        double lt_start, lt_end;
        lt_start = rtclock();
        /* grunCorr(symmat_outputFromGpu, symmat_outputFromGpu, stddev_Gpu, mean_Gpu, float_n, eps); */
        if((env = getenv("OMP_CTSAR_FALLBACK")) && atoi(env)){
            conv3D(A, B_outputFromGpu);
        }else{
            ctsar_init(&s,NI,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
            conv3D_gpu(a_base, bg_base);
        }
        lt_end = rtclock();
        fprintf(stderr, "GPU Runtime:it%d: %0.6lfs\n", i, lt_end - lt_start);

    }
    t_end = rtclock();
    fprintf(stderr, "GPU Runtime: %0.6lfs\n", t_end - t_start);

    if((env = getenv("TEST_VERIFY")) && atoi(env)){
        t_start = rtclock();

        conv3D(A, B);

        t_end = rtclock();
        fprintf(stderr, "CPU Runtime: %0.6lfs\n", t_end - t_start);

        compareResults(B, B_outputFromGpu);
    }

    return 0;
}
