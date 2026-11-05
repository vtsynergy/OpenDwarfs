/**
 * gemm.c: This file is part of the PolyBench/GPU 1.0 test suite.
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
#include <omp.h>

#include "polybenchUtilFuncts.h"
#include "ctsar.h"
ctsar * s = NULL;

//define the error threshold for the results "not matching"
#define PERCENT_DIFF_ERROR_THRESHOLD 0.05

/* Problem size. */

/* #define NI 20240 */
/* #define NJ 20240 */
/* #define NK 20240 */
/* #define N  20240 */

static int  NI = 1024;
static int  NJ = 1024;
static int  NK = 1024;
static int  N  = 1024;


/* Can switch DATA_TYPE between float and double */
typedef double DATA_TYPE;

#define p_alpha 32412
#define p_beta 2123

void runGemm_gpu(DATA_TYPE **a_a, DATA_TYPE **b_a, DATA_TYPE **c_a)
{
    int i, j, k;
    DATA_TYPE *a = *a_a;
    DATA_TYPE *b = *b_a;
    DATA_TYPE *c = *c_a;

/*#pragma ctsar hetero(1) pcopyin(a[1:NJ][0:NI]) pcopy(c[1:NJ][0:NI]) copyin(b[1:NJ*NI]) | acc parallel for gang*/
    /*for (int i = 0; i < NI; ++i)*/
    /*{*/
/*#pragma acc for vector*/
        /*for (int j = 0; j < NJ; ++j)*/
        /*{*/
            /*c[(i * NJ) + j] *= p_beta;*/
/*#pragma acc for seq*/
            /*for (int k = 0; k < NK; ++k)*/
            /*{*/
                /*c[(i * NJ) + j] += p_alpha * a[(i*NJ) + k] * b[(k*NJ) + j];*/
            /*}*/
        /*}*/
    /*}*/
    static ctsar *CTSAR_INST_0 = NULL;
    ctsar_init(&CTSAR_INST_0, NI-0, CTSAR_STATIC, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
#pragma omp parallel
    {
        int __CTSAR_TID__ = omp_get_thread_num();
        do{
#pragma acc data region if(ctsar_get_type(CTSAR_INST_0) == CTSAR_DEV_GPU)
            {
                ctsar_next(CTSAR_INST_0,NI-0);
                int __CTSAR_START__ = CSTART(CTSAR_INST_0,__CTSAR_TID__);
                int __CTSAR_END__ = CEND(CTSAR_INST_0,__CTSAR_TID__);
                size_t a_width, b_width, c_width;

                DATA_TYPE * restrict __CTSAR__a_0__ = (DATA_TYPE *)ctsar_reg_mem_2d(CTSAR_INST_0, a, sizeof(DATA_TYPE), NJ, NI, CTSAR_MEM_RESHAPE | CTSAR_MEM_INPUT | CTSAR_MEM_PARTIAL | CTSAR_MEM_PACKED, 0, 0, &a_width);
                DATA_TYPE * restrict __CTSAR__c_1__ = (DATA_TYPE *)ctsar_reg_mem_2d(CTSAR_INST_0, c, sizeof(DATA_TYPE), NJ, NI, CTSAR_MEM_RESHAPE | CTSAR_MEM_INPUT | CTSAR_MEM_OUTPUT | CTSAR_MEM_PARTIAL | CTSAR_MEM_PACKED, 0, 0, &c_width);
                DATA_TYPE * restrict __CTSAR__b_2__ = (DATA_TYPE *)ctsar_reg_mem_2d(CTSAR_INST_0, b, sizeof(DATA_TYPE), NJ, NI, CTSAR_MEM_RESHAPE | CTSAR_MEM_INPUT /* | CTSAR_MEM_PERSIST */, 0, 0, &b_width);
                uint32_t a_i_width = a_width;
                uint32_t b_i_width = b_width;
                uint32_t c_i_width = c_width;
                if(ctsar_get_type(CTSAR_INST_0) == CTSAR_DEV_GPU){
                    ctsar_start(CTSAR_INST_0);
#pragma acc kernels for independent deviceptr(__CTSAR__a_0__) deviceptr(__CTSAR__c_1__) deviceptr(__CTSAR__b_2__) if(ctsar_get_type(CTSAR_INST_0) == CTSAR_DEV_GPU)
                    for( int i = __CTSAR_START__; i < __CTSAR_END__; ++ i )
                    {
                      int c_off = i * c_i_width;
                      int a_off = i * a_i_width;
#pragma acc for independent
                        for (int j = 0; j < NJ; ++j)
                        {
                          DATA_TYPE loc = __CTSAR__c_1__[c_off + j] * p_beta;
#pragma acc for reduction(+:loc)
                          for (int k = 0; k < NK; ++k)
                          {
                            loc += p_alpha * __CTSAR__a_0__[a_off + k] * __CTSAR__b_2__[(k*b_i_width) + j];
                          }
                          __CTSAR__c_1__[c_off + j] = loc;
                        }
                    }
                    ctsar_end(CTSAR_INST_0);
                }
                else{
                    ctsar_start(CTSAR_INST_0);
                    for( int i = __CTSAR_START__; i < __CTSAR_END__; ++ i )
                    {
                        for (int j = 0; j < NJ; ++j)
                        {
                            __CTSAR__c_1__[(i * NJ) + j] *= p_beta;
                            for (int k = 0; k < NK; ++k)
                            {
                                __CTSAR__c_1__[(i * NJ) + j] += p_alpha * __CTSAR__a_0__[(i*NJ) + k] * __CTSAR__b_2__[(k*NJ) + j];
                            }
                        }
                    }
                    ctsar_end(CTSAR_INST_0);
                }
            }
        }while(ctsar_loop(CTSAR_INST_0));
    }
}

void runGemm(DATA_TYPE **a_a, DATA_TYPE **b_a, DATA_TYPE **c_a)
{
    int i, j, k;
    DATA_TYPE *restrict a = *a_a;
    DATA_TYPE *restrict b = *b_a;
    DATA_TYPE *restrict c = *c_a;

    /* C := alpha*A*B + beta*C */
#pragma omp parallel for private(i,j,k) schedule(runtime)
/* #pragma acc kernels for independent copyin(a[0:NI*NJ],b[0:NI*NJ]) copy(c[0:NI*NJ]) */
    for (i = 0; i < NI; i++)
    {
/* #pragma acc for independent */
        for (int j = 0; j < NJ; ++j)
        {
          DATA_TYPE loc = c[(i * NJ) + j] * p_beta;
/* #pragma acc for reduction(+:loc) */
            for (int k = 0; k < NK; ++k)
            {
                loc += p_alpha * a[(i*NJ) + k] * b[(k*NJ) + j];
            }
            c[(i * NJ) + j] = loc;
        }
    }
}



void init(DATA_TYPE **A, DATA_TYPE **B, DATA_TYPE **C, DATA_TYPE **C_outputFromGpu)
{
    int i, j;

    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NK; j++)
        {
            A[i][j] = ((DATA_TYPE) i*j) / NI;
        }
    }

    for (i = 0; i < NK; i++)
    {
        for (j = 0; j < NJ; j++)
        {
            B[i][j] = ((DATA_TYPE) i*j + 1) / NJ;
        }
    }

    for (i = 0; i < NI; i++)
    {
        for (j = 0; j < NJ; j++)
        {
            C[i][j] = ((DATA_TYPE) i*j + 2) / NJ;
            C_outputFromGpu[i][j] = ((DATA_TYPE) i*j + 2) / NJ;
        }
    }
}


void compareResults(DATA_TYPE **C, DATA_TYPE **C_outputFromGpu)
{
    int i, j, fail;
    fail = 0;

    // Compare C1 and C2
    for (i=0; i < NI; i++) 
    {
        for (j=0; j < NJ; j++) 
        {
            if (percentDiff(C[i][j], C_outputFromGpu[i][j]) > PERCENT_DIFF_ERROR_THRESHOLD) 
            {
                if(fail < 10){
                    printf("CPU: %f GPU: %f\n", C[i][j], C_outputFromGpu[i][j]);
                }
                fail++;
            }
        }
    }

    // Print results
    printf("Non-Matching CPU-GPU Outputs Beyond Error Threshold of %4.2f Percent: %d\n", PERCENT_DIFF_ERROR_THRESHOLD, fail);

}



int main(int argc, char** argv)
{
    double t_start, t_end;
    char * env = NULL;
    if((env = getenv("TEST_SIZE")) && atoi(env) > 0){
      NI = atoi(env);
      NJ = atoi(env);
      NK = atoi(env);
      N  = atoi(env);
    }

    /* Array declaration */
    DATA_TYPE alpha;
    DATA_TYPE beta;
    DATA_TYPE **C = (DATA_TYPE**)ctsar_alloc_2d(sizeof(DATA_TYPE),NI,NJ);
    DATA_TYPE **C_outputFromGpu = (DATA_TYPE**)ctsar_alloc_2d(sizeof(DATA_TYPE),NI,NJ);

    DATA_TYPE **A = (DATA_TYPE**)ctsar_alloc_2d(sizeof(DATA_TYPE),NI,NK);
    DATA_TYPE **B = (DATA_TYPE**)ctsar_alloc_2d(sizeof(DATA_TYPE),NK,NJ);

    /* Initialize array. */
    init(A, B, C, C_outputFromGpu);
    ctsar_pre_init();

    int iters;
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
            runGemm(A, B, C_outputFromGpu);
	}else{
            ctsar_init(&s,NI,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
            runGemm_gpu(A, B, C_outputFromGpu);
        }
        lt_end = rtclock();
        fprintf(stderr, "GPU Runtime:it%d: %0.6lfs\n", i, lt_end - lt_start);

    }


    t_end = rtclock();
    fprintf(stderr, "GPU Runtime: %0.6lfs\n", t_end - t_start);


    if((env = getenv("TEST_VERIFY")) && atoi(env)){
        t_start = rtclock();

        runGemm(A, B, C);

        t_end = rtclock();
        fprintf(stderr, "CPU Runtime: %0.6lfs\n", t_end - t_start);

        compareResults(C, C_outputFromGpu);
    }

    return 0;
}
