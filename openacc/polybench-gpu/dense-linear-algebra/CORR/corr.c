/**
 * corr.c: This file is part of the PolyBench/GPU 1.0 test suite.
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
#include <accel.h>

#include "polybenchUtilFuncts.h"
#include "ctsar.h"

//define the error threshold for the results "not matching"
#define PERCENT_DIFF_ERROR_THRESHOLD 0.05
#define sqrt_of_array_cell(x,j) sqrt(x[j])

/* Problem size. */
#define M 2048
#define N 2048
/* #define M 4096
#define N 4096 */

/* Can switch DATA_TYPE between float and double */
typedef float DATA_TYPE;


void prepareCorr(DATA_TYPE **pdata, DATA_TYPE **psymmat, DATA_TYPE *pstddev, DATA_TYPE *pmean, DATA_TYPE pfloat_n, DATA_TYPE peps){
    int i,j;
#pragma omp parallel for private(i,j)
    for (j = 1; j <= M; j++)
    {
        pmean[j] = 0.0;

        for (i = 1; i <= N; i++)
        {
            pmean[j] += pdata[i][j];
        }
        pmean[j] /= pfloat_n;
    }

    /* Determine standard deviations of column vectors of data matrix. */
#pragma omp parallel for private(i,j)
    for (j = 1; j <= M; j++)
    {
        pstddev[j] = 0.0;

        for (i = 1; i <= N; i++)
        {
            pstddev[j] += (pdata[i][j] - pmean[j]) * (pdata[i][j] - pmean[j]);
        }
        pstddev[j] /= pfloat_n;
        pstddev[j] = sqrt_of_array_cell(pstddev, j);

        /* The following in an inelegant but usual way to handle
           near-zero std. dev. values, which below would cause a zero-
           divide. */

        pstddev[j] = pstddev[j] <= peps ? 1.0 : pstddev[j];
    }



    /* Center and reduce the column vectors. */
#pragma omp parallel for private(i,j)
    for (i = 1; i <= N; i++)
    {
        for (j = 1; j <= M; j++)
        {
            pdata[i][j] -= pmean[j];
            pdata[i][j] /= sqrt(pfloat_n) * pstddev[j];
        }
    }

}

void grunCorr(ctsar * s, DATA_TYPE **pdata, DATA_TYPE **psymmat, DATA_TYPE *pstddev, DATA_TYPE *pmean, DATA_TYPE pfloat_n, DATA_TYPE peps)
{
    int i, j, j1, j2;




    /* Determine mean of column vectors of input data matrix */
#pragma omp parallel default(shared) private(j,i,j1,j2)
    {
        int tid = omp_get_thread_num();
        float * lpdata, *lpsymmat;
        lpdata = ctsar_reg_mem(s, pdata[0], sizeof(float)*(M+1), (N+1),
                CTSAR_MEM_INPUT);
        lpsymmat = ctsar_reg_mem(s, psymmat[0], sizeof(float)*(M+1), (N+1),
                CTSAR_MEM_PARTIAL | CTSAR_MEM_OUTPUT);

        do{
            ctsar_next(s,M+1);
            int gts = CSTART(s,tid) ? CSTART(s,tid) : 1;
            int gte = CEND(s,tid) == M+1 ? M : CEND(s,tid);

            ctsar_start(s);
#pragma acc data region deviceptr(lpdata,lpsymmat)  if(ctsar_get_type(s) == CTSAR_DEV_GPU)
            {
#pragma acc region for independent private(j1,j2,i) if(ctsar_get_type(s) == CTSAR_DEV_GPU)
                for (j1 = gts; j1 < gte; j1++)
                {
#pragma acc for independent
                    for (j2 = 2; j2 <= M; j2++)
                    {
                        if(j2 > j1){
                            DATA_TYPE ag = 0.0f;
#pragma acc for seq
                            for (i = 1; i <= N; i++)
                            {
                                ag += (lpdata[i*(M+1)+j1] * lpdata[i*(M+1)+j2]);
                            }
                            lpsymmat[j1*(M+1)+j2] = ag;
                            /* lpsymmat[j2*(M+1)+j1] = ag; */
                        }
                    }
                }
            }
            ctsar_end(s);
            /* printf("here: %d\n", omp_get_thread_num()); */
        }while(ctsar_loop(s));
    }
    psymmat[M][M] = 1.0;

#pragma omp parallel for private(j1,j2) schedule(runtime)
    for (j1 = 1; j1 <= M-1; j1++)
    {
        psymmat[j1][j1] = 1.0;
        for (j2 = j1+1; j2 <= M; j2++)
        {
            psymmat[j2][j1] = psymmat[j1][j2];
        }
    }
}

void runCorr(DATA_TYPE **pdata, DATA_TYPE **psymmat, DATA_TYPE *pstddev, DATA_TYPE *pmean, DATA_TYPE pfloat_n, DATA_TYPE peps)
{
    int i, j, j1, j2;

    /* Calculate the m * m correlation matrix. */
#pragma omp parallel for private(i,j1,j2) schedule(runtime)
    for (j1 = 1; j1 <= M-1; j1++)
    {
        psymmat[j1][j1] = 1.0;

        for (j2 = j1+1; j2 <= M; j2++)
        {
            DATA_TYPE ag = 0.0;

            for (i = 1; i <= N; i++)
            {
                ag += (pdata[i][j1] * pdata[i][j2]);
            }
            psymmat[j1][j2] = ag;
            psymmat[j2][j1] = psymmat[j1][j2];
        }
    }

    psymmat[M][M] = 1.0;
}


void init_arrays(DATA_TYPE **data, DATA_TYPE **data_Gpu)
{
    int i, j;

    for (i=0; i < (M+1); i++)
    {
        for (j=0; j< (N+1); j++)
        {
            data[i][j] = ((DATA_TYPE) i*j)/ (M+1);
            data_Gpu[i][j] = ((DATA_TYPE) i*j)/ (M+1);
        }
    }
}

void compareResults(DATA_TYPE **symmat, DATA_TYPE **symmat_outputFromGpu)
{
    int i,j,fail;
    fail = 0;

    for (i=1; i < (M+1); i++)
    {
        for (j=1; j < (N+1); j++)
        {
            if (percentDiff(symmat[i][j], symmat_outputFromGpu[i][j]) > PERCENT_DIFF_ERROR_THRESHOLD)
            {
                fail++;
                if(fail < 20)
                    printf("i: %d j: %d\n1: %f 2: %f\n", i, j, symmat[i][j], symmat_outputFromGpu[i][j]);
            }
        }
    }

    // print results
    printf("Non-Matching CPU-GPU Outputs Beyond Error Threshold of %4.2f Percent: %d\n", PERCENT_DIFF_ERROR_THRESHOLD, fail);
}


int main(int argc, char** argv)
{
    int m = M;
    int n = N;
    double t_start, t_end;

    /* Array declaration */
    DATA_TYPE float_n = 321414134.01;
    DATA_TYPE eps = 0.005;
    DATA_TYPE *data[M + 1];
    DATA_TYPE *data_Gpu[M + 1];
    DATA_TYPE mean[M + 1];
    DATA_TYPE mean_Gpu[M + 1];
    DATA_TYPE stddev[M + 1];
    DATA_TYPE stddev_Gpu[M + 1];
    DATA_TYPE *symmat[M + 1];
    DATA_TYPE *symmat_outputFromGpu[M + 1];
    /* DATA_TYPE **data_Gpu = (DATA_TYPE **)malloc(sizeof(DATA_TYPE*) * (M + 1)); */
    data[0] = (DATA_TYPE *)calloc(sizeof(DATA_TYPE) ,(M + 1)*(N + 1));
    data_Gpu[0] = (DATA_TYPE *)calloc(sizeof(DATA_TYPE) ,(M + 1)*(N + 1));
    symmat[0] = (DATA_TYPE *)calloc(sizeof(DATA_TYPE) ,(M + 1)*(N + 1));
    symmat_outputFromGpu[0] = (DATA_TYPE *)calloc(sizeof(DATA_TYPE) ,(M + 1)*(N + 1));
    int i;
    for(i=1; i < M+1; i++){
        data_Gpu[i] = data_Gpu[i-1]+N+1;
        data[i] = data[i-1]+N+1;
        symmat[i] = symmat[i-1]+N+1;
        symmat_outputFromGpu[i] = symmat_outputFromGpu[i-1]+N+1;
    }
    printf("test, printing %f, size %d, should be %d\n", data_Gpu[N][M], (sizeof(DATA_TYPE) *(M + 1)*(N + 1)), (M+1)*(N+1)*4);

    /* Initialize array. */
    init_arrays(data, data_Gpu);

    ctsar_pre_init();

    /* #pragma acc data region copyin(data_Gpu, symmat_outputFromGpu, stddev_Gpu, mean_Gpu, float_n, eps) */
    prepareCorr(data, symmat, stddev, mean, float_n, eps);
    prepareCorr(data_Gpu, symmat_outputFromGpu, stddev, mean, float_n, eps);

    int iters;
    char * env;
    if(env = getenv("TEST_ITERATIONS")){
        iters = atoi(env);
    }else{
        iters = 1;
    }
    t_start = rtclock();
    ctsar * s = NULL;
    ctsar_init(&s,M+1,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
    for(i=0; i < iters; i++)
    {
        double lt_start, lt_end;
        lt_start = rtclock();
        /* grunCorr(symmat_outputFromGpu, symmat_outputFromGpu, stddev_Gpu, mean_Gpu, float_n, eps); */
        if((env = getenv("OMP_CTSAR_FALLBACK")) && atoi(env)){
            runCorr(data, symmat, stddev, mean, float_n, eps);
        }else{
            ctsar_init(&s,M+1,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
            grunCorr(s,data_Gpu, symmat_outputFromGpu, stddev_Gpu, mean_Gpu, float_n, eps);
        }
        lt_end = rtclock();
        fprintf(stderr, "GPU Runtime:it%d: %0.6lfs\n", i, lt_end - lt_start);
    }
    t_end = rtclock();
    fprintf(stderr, "GPU Runtime:total: %0.6lfs\n", t_end - t_start);
    

    if((env = getenv("TEST_VERIFY")) && atoi(env)){
        t_start = rtclock();
        /* runCorr(symmat, symmat, stddev, mean, float_n, eps); */
        runCorr(data, symmat, stddev, mean, float_n, eps);

        t_end = rtclock();
        fprintf(stderr, "CPU Runtime: %0.6lfs\n", t_end - t_start);

        compareResults(symmat, symmat_outputFromGpu);
    }
    return 0;
}
