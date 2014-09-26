/**
 * covar.c: This file is part of the PolyBench/GPU 1.0 test suite.
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
#include <stdint.h>

#include "polybenchUtilFuncts.h"
#include "ctsar.h"
ctsar * s = NULL;


static struct timeval start_time;
static void init_time() {
    gettimeofday(&start_time, NULL);
}

static uint64_t get_time() {
    struct timeval t;
    gettimeofday(&t, NULL);
    return (uint64_t) (t.tv_sec - start_time.tv_sec) * 1000000
        + (t.tv_usec - start_time.tv_usec);
}

//define the error threshold for the results "not matching"
#define PERCENT_DIFF_ERROR_THRESHOLD 0.05

/* Problem size. */
/* #define M 2048
#define N 2048 */
int M=2048;
int N=2048;

/* Can switch DATA_TYPE between float and double */
typedef float DATA_TYPE;





void covarLoopa(DATA_TYPE *pmean, DATA_TYPE **pdata, DATA_TYPE pfloat_n)
{
    int i, j, j1, j2;

    /* Determine mean of column vectors of input data matrix */

    for (j = 1; j <= M; j++)
    {
        pmean[j] = 0.0;

        for (i = 1; i <= N; i++)
        {
            pmean[j] += pdata[i][j];
        }
        pmean[j] /= pfloat_n;
    }
}

void covarLoopb(DATA_TYPE **pdata, DATA_TYPE *pmean)
{
    int i, j;

    /* Center the column vectors. */

    for (i = 1; i <= N; i++)
    {
        for (j = 1; j <= M; j++)
        {
            pdata[i][j] -= pmean[j];
        }
    }
}

void covarLoopc_gpu(DATA_TYPE **psymmat, DATA_TYPE **pdata)
{
    int i, j1, j2;
    init_time();

    /* Calculate the m * m covariance matrix. */


#pragma omp parallel default(shared) private(j1,j2,i)
    {
        int tid = omp_get_thread_num();
        DATA_TYPE *d,*sy;
        d = ctsar_reg_mem(s, *pdata, sizeof(DATA_TYPE)*(M+1), M+1,
                CTSAR_MEM_INPUT);
        sy = ctsar_reg_mem(s, *psymmat, sizeof(DATA_TYPE)*(M+1), M+1,
                CTSAR_MEM_PARTIAL | CTSAR_MEM_OUTPUT);

        do{
            ctsar_next(s,M+1);
            int gts = CSTART(s,tid) ? CSTART(s,tid) : 1;
            int gte = CEND(s,tid);

            ctsar_start(s);
#pragma acc data region deviceptr(d,sy)  if(ctsar_get_type(s) == CTSAR_DEV_GPU)
            {
#pragma acc region for private(j1,j2,i) independent   if(ctsar_get_type(s) == CTSAR_DEV_GPU)
                for (j1 = gts; j1 < gte; j1++)
                {
#pragma acc for independent
                    for (j2 = 1; j2 <= M; j2++)
                    {
                        if( j2 >= j1 ){
                            DATA_TYPE l = 0.0;
#pragma acc for reduction(+:l)
                            for (i = 1; i <= N; i++) {
                                l += d[((N+1) *i) + j1] * d[((N+1) *i) + j2];
                            }
                            sy[((N+1)*j1)+j2] = l;
                        }
                    }
                }
            }
            ctsar_end(s);
            /* printf("here: %d\n", omp_get_thread_num()); */
        }while(ctsar_loop(s));
    }

    printf("covarc took %lu\n", get_time());

    for (j1 = 1; j1 <= M; j1++)
    {
        for (j2 = j1; j2 <= M; j2++)
        {
            psymmat[j2][j1] = psymmat[j1][j2];
        }
    }
}
void covarLoopc(DATA_TYPE **psymmat, DATA_TYPE **pdata)
{
    int i, j1, j2;
    init_time();

    /* Calculate the m * m covariance matrix. */


#pragma omp parallel for private(j1,j2,i) schedule(runtime)
    for (j1 = 1; j1 <= M; j1++)
    {
        for (j2 = j1; j2 <= M; j2++)
        {
            psymmat[j1][j2] = 0.0;

            for (i = 1; i <= N; i++)
            {
                psymmat[j1][j2] += pdata[i][j1] * pdata[i][j2];
            }

            psymmat[j2][j1] = psymmat[j1][j2];
        }
    }
    printf("covarc took %lu\n", get_time());
}


void init_arrays(DATA_TYPE **data, DATA_TYPE **data_Gpu)
{
    int i, j;

    for (i = 1; i < (M+1); i++)
    {
        for (j = 1; j < (N+1); j++)
        {
            data[i][j] = ((DATA_TYPE) i*j) / M;
            data_Gpu[i][j] = ((DATA_TYPE) i*j) / M;
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
                if(fail < 10)
                    printf("i=%d, j=%d, symmat:%f gpu:%f\n",i,j,symmat[i][j], symmat_outputFromGpu[i][j]);
                fail++;
            }			
        }
    }
    printf("Non-Matching CPU-GPU Outputs Beyond Error Threshold of %4.2f Percent: %d\n", PERCENT_DIFF_ERROR_THRESHOLD, fail);

}


int main(int argc, char** argv)
{
    char * env;
    if((env = getenv("TEST_SIZE")) && atoi(env) > 0){
      N  = atoi(env);
      M  = atoi(env);
    }


    int iters;
    if(env = getenv("TEST_ITERATIONS")){
        iters = atoi(env);
    }else{
        iters = 1;
    }

    int m = M;
    int n = N;
    double t_start, t_end;

    /* Array declaration */
    DATA_TYPE float_n = 321414134.01;
    DATA_TYPE **data = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),M + 1,N + 1);
    DATA_TYPE **data_Gpu = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),M + 1,N + 1);
    DATA_TYPE **symmat = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),M + 1,M + 1);
    DATA_TYPE **symmat_outputFromGpu = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),M + 1,M + 1);	
    DATA_TYPE *mean = (DATA_TYPE *)ctsar_calloc_test(sizeof(DATA_TYPE),M + 1);
    DATA_TYPE *mean_Gpu = (DATA_TYPE *)ctsar_calloc_test(sizeof(DATA_TYPE),M + 1);

    /* Initialize array. */
    init_arrays(data, data_Gpu);

    ctsar_pre_init();

    t_start = rtclock();
    ctsar_init(&s,M+1,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);

    covarLoopa(mean_Gpu, data_Gpu, float_n);
    covarLoopb(data_Gpu, mean_Gpu);
    for(int i=0; i < iters; i++)
    {
        double lt_start, lt_end;
        lt_start = rtclock();
        /* grunCorr(symmat_outputFromGpu, symmat_outputFromGpu, stddev_Gpu, mean_Gpu, float_n, eps); */
        if((env = getenv("OMP_CTSAR_FALLBACK")) && atoi(env)){
            covarLoopc(symmat_outputFromGpu, data_Gpu);
        }else{
            ctsar_init(&s,M+1,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
            covarLoopc_gpu(symmat_outputFromGpu, data_Gpu);
        }
        lt_end = rtclock();
        fprintf(stderr, "GPU Runtime:it%d: %0.6lfs\n", i, lt_end - lt_start);
    }

    t_end = rtclock();
    fprintf(stderr, "GPU Runtime: %0.6lfs\n", t_end - t_start);




    if((env = getenv("TEST_VERIFY")) && atoi(env)){
        t_start = rtclock();

        covarLoopa(mean, data, float_n);
        covarLoopb(data, mean);
        covarLoopc(symmat, data);

        t_end = rtclock();
        fprintf(stderr, "CPU Runtime: %0.6lfs\n", t_end - t_start);

        compareResults(symmat, symmat_outputFromGpu);
    }

    return 0;
}
