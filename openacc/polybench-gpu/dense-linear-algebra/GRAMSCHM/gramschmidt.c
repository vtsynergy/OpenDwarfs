/**
 * gramschmidt.c: This file is part of the PolyBench/GPU 1.0 test suite.
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
#define M 1024
#define N 1024

/* Can switch DATA_TYPE between float and double */
typedef float DATA_TYPE;


void transpose(DATA_TYPE **in, int x, int y)
{
    DATA_TYPE tmp;
/* #pragma omp parallel for */
    for(int i=0; i < y; ++i){
        for(int j=0; j < i; ++j){
            tmp = in[i][j];
            in[i][j] = in[j][i];
            in[j][i] = tmp;
        }
    }
}
void grunGramSchmidt(DATA_TYPE **pAa, DATA_TYPE **pRa, DATA_TYPE **pQa)
{
    int i, j, k;
    int m = M;
    int n = N;



/* #pragma omp parallel default(shared) private(j,i,k)
    {
        int tid = omp_get_thread_num();

        do{
            ctsar_next(s,N);
            int gts = CSTART(s,tid);
            int gte = CEND(s,tid);
            DATA_TYPE *pA = *pAa,
                      *pR = *pRa,
                      *pQ = *pQa;
            pA = ctsar_reg_mem(s, *pAa, sizeof(DATA_TYPE)*M, N,
                    CTSAR_MEM_PARTIAL | CTSAR_MEM_INPUT | CTSAR_MEM_OUTPUT);
            pR = ctsar_reg_mem(s, *pRa, sizeof(DATA_TYPE)*M, N,
                    CTSAR_MEM_PERSIST);
            pQ = ctsar_reg_mem(s, *pQa, sizeof(DATA_TYPE)*M, N,
                    CTSAR_MEM_PERSIST);

            ctsar_start(s);
#pragma acc data region deviceptr(pA,pR,pQ)  if(ctsar_get_type(s) == CTSAR_DEV_GPU)
            {
#pragma acc region for private(i,j,k) if(ctsar_get_type(s) == CTSAR_DEV_GPU)
                for (k = gts; k < gte; k++)
                {
                    DATA_TYPE pnrm = 0;

                    for (i = 0; i < m; i++)
                    {
                        pnrm += pA[(i * N) + k] * pA[(i * N) + k];
                    }

                    pR[(k * N) + k] = sqrtf(pnrm);

                    for (i = 0; i < m; i++)
                    {
                        pQ[(i * N) + k] = pA[(i * N) + k] / pR[(k * N) + k];
                    }

                    for (j = k + 1; j < n; j++)
                    {
                        pR[(k * N) + j] = 0;

                        for (i = 0; i < m; i++)
                        {
                            pR[(k * N) + j] += pQ[(i * N) + k] * pA[(i * N) + j];
                        }

                        for (i = 0; i < m; i++)
                        {
                            pA[(i * N) + j] = pA[(i * N) + j] - pQ[(i * N) + k] * pR[(k * N) + j];
                        }
                    }
                }
            }
            ctsar_end(s);
            /\* printf("here: %d\n", omp_get_thread_num()); *\/
        }while(ctsar_loop(s));
    } */
    /* transpose(pAa, M, N);
    transpose(pRa, M, N);
    transpose(pQa, M, N); */
    DATA_TYPE *pA = *pAa,
              *pR = *pRa,
              *pQ = *pQa;
    /* for (k = 0; k < n; k++)
    { */

/* #pragma acc region for private(i,j,k) copyin(pA[0:M*N]) copyin(pR[0:M*N], pQ[0:M*N]) */
/* #pragma acc data region copy(pA[0:M*N]) local(pR[0:M*N], pQ[0:M*N])
        {
#pragma acc region
            {
                DATA_TYPE pnrm = 0;
#pragma acc for private(i)
                for (i = 0; i < m; i++)
                {
                    pnrm += pA[(i * N) + k] * pA[(i * N) + k];
                }

                pR[(k * N) + k] = sqrtf(pnrm);
            }

#pragma acc region for private(i) independent
            for (i = 0; i < m; i++)
            {
                pQ[(i * N) + k] = pA[(i * N) + k] / pR[(k * N) + k];
            }

#pragma acc region for private(i,j) independent
            for (j = k + 1; j < n; j++)
            {
                DATA_TYPE pRl = 0;

                for (i = 0; i < m; i++)
                {
                    pRl += pQ[(i * N) + k] * pA[(i * N) + j];
                }

                pR[(k * N) + j] = pRl;

                for (i = 0; i < m; i++)
                {
                    pA[(i * N) + j] = pA[(i * N) + j] - pQ[(i * N) + k] * pR[(k * N) + j];
                }
            }
        }
    } */
    transpose(pAa,M,N);
/* #pragma omp parallel for private(k,i,j) schedule(runtime) */
    for (k = 0; k < n; k++)
    {
        DATA_TYPE   pnrm = 0;

#pragma omp parallel for private(i,j) schedule(runtime) reduction(+:pnrm)
        for (i = 0; i < m; i++)
        {
            pnrm += pA[(k * N) + i] * pA[(k * N) + i];
        }

        pR[(k * N) + k] = sqrtf(pnrm);

#pragma omp parallel for private(i,j) schedule(runtime)
        for (i = 0; i < m; i++)
        {
            pQ[(k * N) + i] = pA[(k * N) + i] / pR[(k * N) + k];
        }

        ctsar_init(&s,N,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
#pragma omp parallel default(shared) private(j,i)
        {
            int tid = omp_get_thread_num();
            DATA_TYPE *plA = *pAa,
                      *plR = *pRa,
                      *plQ = *pQa;

/* #pragma acc data create(plA[0:M*N]) if(ctsar_get_type(s) == CTSAR_DEV_GPU)  */
            do{
                ctsar_next(s,N);
                int gts = CSTART(s,tid);
                int gte = CEND(s,tid);
                /* fprintf(stderr,"[%d] gts=%d, gte=%d\n",tid,gts,gte); */
                plA = ctsar_reg_mem(s, *pAa, sizeof(DATA_TYPE)*M, N,
                        CTSAR_MEM_PARTIAL | CTSAR_MEM_OUTPUT);
                if(k==0)
                    acc_memcpy_from_device(plA, *pAa, sizeof(DATA_TYPE)*M*N);
                /* plR = ctsar_reg_mem(s, *pRa, sizeof(DATA_TYPE)*M, N,
                        CTSAR_MEM_PERSIST); */
                if(k>0) ctsar_retarget_mem(s, pQa[k-1], pQa[k]);
                plQ = ctsar_reg_mem(s, pQa[k], sizeof(DATA_TYPE), M,
                        CTSAR_MEM_INPUT);

/* #pragma acc update device(plA[(gts)*M:(gte)*M])  if(ctsar_get_type(s) == CTSAR_DEV_GPU)  */
                ctsar_start(s);
#pragma acc data if(ctsar_get_type(s) == CTSAR_DEV_GPU) deviceptr(plA,plR,plQ)
                //deviceptr(plA,plR,plQ)
                {
#pragma acc region for independent private(i,j) if(ctsar_get_type(s) == CTSAR_DEV_GPU)
                    for (int j = gts; j < gte; j++)
                    {
                            DATA_TYPE lr = 0;
                        if(j >= k+1){
                            /* plR[(j * N) + k] = 0; */

#pragma acc for seq
                            for (i = 0; i < m; i++)
                            {
                                lr += plQ[i] * plA[(j * N) + i];
                            }
#pragma acc for seq
                            for (i = 0; i < m; i++)
                            {
                                plA[(j * N) + i] -= plQ[i] * lr;
                            }
                        }else{
#pragma acc for seq
                            for (i = 0; i < m; i++)
                            {
                                lr +=  plA[(j * N) + i];
                            }

                        }
                    }
                }
                ctsar_end(s);
            }while(ctsar_loop(s));
        }
    }
    transpose(pAa,M,N);
}

void runGramSchmidt(DATA_TYPE **pA, DATA_TYPE **pR, DATA_TYPE **pQ)
{
    int i, j, k;
    int m = M;
    int n = N;

  /* for (k = 0; k < _PB_NJ; k++)
    {
      nrm = 0;
      for (i = 0; i < _PB_NI; i++)
        nrm += A[i][k] * A[i][k];
      R[k][k] = sqrt(nrm);
      for (i = 0; i < _PB_NI; i++)
        Q[i][k] = A[i][k] / R[k][k];
      for (j = k + 1; j < _PB_NJ; j++)
        {
          R[k][j] = 0;
          for (i = 0; i < _PB_NI; i++)
            R[k][j] += Q[i][k] * A[i][j];
          for (i = 0; i < _PB_NI; i++)
            A[i][j] = A[i][j] - Q[i][k] * R[k][j];
        }
    } */


    /* {//working
#pragma omp parallel for private(k,i,j) schedule(runtime)
    for (k = 0; k < n; k++)
    {
        DATA_TYPE   pnrm = 0;

        for (i = 0; i < m; i++)
        {
            pnrm += pA[i][k] * pA[i][k];
        }

        pR[k][k] = sqrtf(pnrm);

        for (i = 0; i < m; i++)
        {
            pQ[i][k] = pA[i][k] / pR[k][k];
        }

        for (j = k + 1; j < n; j++)
        {
            pR[k][j] = 0;

            for (i = 0; i < m; i++)
            {
                pR[k][j] += pQ[i][k] * pA[i][j];
            }

            for (i = 0; i < m; i++)
            {
                pA[i][j] = pA[i][j] - pQ[i][k] * pR[k][j];
            }
        }
    }
    }*/

    transpose(pA,M,N);
#pragma omp parallel for private(k,i,j) schedule(runtime)
    for (k = 0; k < n; k++)
    {
        DATA_TYPE   pnrm = 0;

        for (i = 0; i < m; i++)
        {
            pnrm += pA[k][i] * pA[k][i];
        }

        pR[k][k] = sqrtf(pnrm);

        for (i = 0; i < m; i++)
        {
            pQ[k][i] = pA[k][i] / pR[k][k];
        }

        for (j = k + 1; j < n; j++)
        {
            pR[j][k] = 0;

            for (i = 0; i < m; i++)
            {
                pR[j][k] += pQ[k][i] * pA[j][i];
            }

            for (i = 0; i < m; i++)
            {
                pA[j][i] = pA[j][i] - pQ[k][i] * pR[j][k];
            }
        }
    }
    transpose(pA,M,N);
}


void init_array(DATA_TYPE **A, DATA_TYPE **A_Gpu)
{
    int i, j;

    for (i = 0; i < M; i++)
    {
        for (j = 0; j < N; j++)
        {
            A[i][j] = ((DATA_TYPE) i*j+M*M) / M;
            A_Gpu[i][j] = ((DATA_TYPE) i*j+M*M) / M;
            /* A[i][j] = ((DATA_TYPE) (i+1)*(j+1)) / (M+1); */
            /* A_Gpu[i][j] = ((DATA_TYPE) (i+1)*(j+1)) / (M+1); */
        }
    }
}


void compareResults(DATA_TYPE **A, DATA_TYPE **A_outputFromGpu)
{
    int i, j, fail;
    fail = 0;

    FILE * co = fopen("./c.out", "w");
    FILE * go = fopen("./g.out", "w");

    for (i=0; i < M; i++)
    {
        for (j=0; j < N; j++)
        {
            fprintf(co, "%f\t", A[i][j]);
            fprintf(go, "%f\t", A_outputFromGpu[i][j]);
            if (percentDiff(A[i][j], A_outputFromGpu[i][j]) > PERCENT_DIFF_ERROR_THRESHOLD)
            {
                if(fail < 10)
                    printf("non-matching at i=%d, j=%d cpu=%f gpu=%f\n", i, j, A[i][j], A_outputFromGpu[i][j]);
                fail++;
            }
        }
        fprintf(co, "\n");
        fprintf(go, "\n");
    }
    fclose(co);
    fclose(go);

    // Print results
    printf("Non-Matching CPU-GPU Outputs Beyond Error Threshold of %4.2f Percent: %d\n", PERCENT_DIFF_ERROR_THRESHOLD, fail);

}


int main(int argc, char** argv)
{
    double t_start, t_end;

    /* Array declaration */
    DATA_TYPE **A = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),M,N);
    DATA_TYPE **A_outputFromGpu = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),M,N);
    DATA_TYPE **R = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),M,N);
    DATA_TYPE **R_Gpu = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),M,N);
    DATA_TYPE **Q = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),M,N);
    DATA_TYPE **Q_Gpu = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),M,N);

    /* Initialize array. */
    init_array(A, A_outputFromGpu);

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
            runGramSchmidt(A_outputFromGpu, R_Gpu, Q_Gpu);
        }else{
            ctsar_init(&s,N,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
            grunGramSchmidt(A_outputFromGpu, R_Gpu, Q_Gpu);
        }
        lt_end = rtclock();
        fprintf(stderr, "GPU Runtime:it%d: %0.6lfs\n", i, lt_end - lt_start);

    }


    t_end = rtclock();
    fprintf(stderr, "GPU Runtime: %0.6lfs\n", t_end - t_start);

    if((env = getenv("TEST_VERIFY")) && atoi(env)){
        t_start = rtclock();

        runGramSchmidt(A, R, Q);

        t_end = rtclock();
        fprintf(stderr, "CPU Runtime: %0.6lfs\n", t_end - t_start);

        compareResults(A, A_outputFromGpu);
    }

    return 0;
}
