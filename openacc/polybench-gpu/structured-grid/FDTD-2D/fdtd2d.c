/**
 * fdtd2d.c: This file is part of the PolyBench/GPU 1.0 test suite.
 *
 *
 * Contact: Scott Grauer-Gray <sgrauerg@gmail.com>
 * Louis-Noel Pouchet <pouchet@cse.ohio-state.edu>
 * Web address: http://www.cse.ohio-state.edu/~pouchet/software/polybench/GPU
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <unistd.h>
#include <sys/time.h>

#include "polybenchUtilFuncts.h"
#include "ctsar.h"
ctsar * s = NULL;
ctsar * c = NULL;

//define the error threshold for the results "not matching"
#define PERCENT_DIFF_ERROR_THRESHOLD 0.05

/* Problem dimensions */
#define NX 2048 * 3
#define NY 2048 * 3
/* #define NX 1024
#define NY 1024 */
#define T_MAX 5

/* Can switch DATA_TYPE between float and double */
typedef float DATA_TYPE;

void grunFdtd(DATA_TYPE *fict, DATA_TYPE **exa, DATA_TYPE **eya, DATA_TYPE **hza)
{
    int t, i, j;
    for(t=0; t < T_MAX; t++)
    {
#pragma omp for
        for (j=0; j < NY; j++)
        {
            eya[0][j] = fict[t];
        }


#pragma omp parallel default(shared) private(j,i)
        {
            int tid = omp_get_thread_num();
            do{
                DATA_TYPE *ex = *exa,
                          *ey = *eya,
                          *hz = *hza;
                ctsar_next(s,NX);
                int gts = CSTART(s,tid);
                int gte = CEND(s,tid);
                ex = ctsar_reg_mem_2d(s, exa[0], sizeof(DATA_TYPE),(NX+1), NX,
                        CTSAR_MEM_PARTIAL | CTSAR_MEM_INPUT | CTSAR_MEM_OUTPUT,0,0,NULL);
                ey = ctsar_reg_mem_2d(s, eya[0], sizeof(DATA_TYPE),NX, NX+1,
                        CTSAR_MEM_PARTIAL | CTSAR_MEM_INPUT | CTSAR_MEM_OUTPUT,0,0,NULL);
                hz = ctsar_reg_mem_2d(s, hza[0], sizeof(DATA_TYPE),NX, NX,
                        CTSAR_MEM_PARTIAL | CTSAR_MEM_INPUT,0,1,NULL);

                /* if(ctsar_get_type(s) == CTSAR_DEV_GPU && gts < gte)
                   cmcpyhd(hz, hza[0], sizeof(DATA_TYPE)*NX*NY); */
                /* cmcpyhd(ex, exa[0], sizeof(DATA_TYPE)*NX*NY);
                   cmcpyhd(ey, eya[0], sizeof(DATA_TYPE)*NX*NY); */
                /* fprintf(stderr, "before 1\n"); */
                ctsar_start(s);
#pragma acc data if(ctsar_get_type(s) == CTSAR_DEV_GPU)
                {
                    /* #pragma acc region for independent copy(ey[0:NX*NX], fict[0:T_MAX]) */
#pragma acc region for independent deviceptr(ey,ex,hz) if(ctsar_get_type(s) == CTSAR_DEV_GPU)
                    for (i = gts; i < gte; i++)
                    {
                        if(i >= 1){
#pragma acc for independent
                            for (j = 0; j < NY; j++)
                            {
                                /* eya[i][j] = eya[i][j] - 0.5*(hza[i][j] - hza[i-1][j]); */
                                ey[((i) * (NX)) + (j)] -= 0.5*(hz[((i) * NX) + (j)] - hz[((i-1) * NX) + (j)]);
                                ex[((i) * (NX+1)) + (j)] -= 0.5*(hz[((i) * NX) + (j)] - hz[((i) * NX) + (j-1)]);
                            }
                        }
                    }
                }
                ctsar_end(s);
                /* printf("here: %d\n", omp_get_thread_num()); */
            }while(ctsar_loop(s));

            do{
                DATA_TYPE *ex = *exa,
                          *ey = *eya,
                          *hz = *hza;
                ctsar_next(c,NX);
                int gts = CSTART(c,tid);
                int gte = CEND(c,tid);
                /* fprintf(stderr,"%d: gts=%d\n", tid, gts); */
                ex = ctsar_reg_mem_2d(c, exa[0], sizeof(DATA_TYPE),(NX+1), NX,
                        CTSAR_MEM_PARTIAL | CTSAR_MEM_INPUT,0,1,NULL);
                ey = ctsar_reg_mem_2d(c, eya[0], sizeof(DATA_TYPE),NX, NX+1,
                        CTSAR_MEM_PARTIAL | CTSAR_MEM_INPUT,0,1,NULL);
                hz = ctsar_reg_mem_2d(c, hza[0], sizeof(DATA_TYPE),NX, NX,
                        CTSAR_MEM_PARTIAL | CTSAR_MEM_INPUT | CTSAR_MEM_OUTPUT,0,1,NULL);
                /* hz = ctsar_reg_mem(c, hza[0], sizeof(DATA_TYPE)*NX, NX,
                   CTSAR_MEM_PARTIAL | CTSAR_MEM_OUTPUT); */
                /* if(ctsar_get_type(c) == CTSAR_DEV_GPU && gts < gte){
                   /\* cmcpydh(exa[0] + (gts*(NX+1)), ex + (gts*(NX+1)), sizeof(DATA_TYPE)*(NX+1)); *\/
                   cmcpydh(eya[0] + (gts*(NX)), ey + (gts*(NX)), sizeof(DATA_TYPE)*(NX));
                   csync();
                   } */
                /* #pragma omp barrier */
                /* if(ctsar_get_type(c) == CTSAR_DEV_GPU && gts < gte){
                   /\* cmcpyhd(ex, exa[0], sizeof(DATA_TYPE)*NX*NY); *\/
                   /\* cmcpyhd(ey, eya[0], sizeof(DATA_TYPE)*NX*NY); *\/
                   cmcpyhd(ey + ((gte)*(NX)), eya[0] + ((gte)*(NX)), sizeof(DATA_TYPE)*(NX));
                   } */
                ctsar_start(c);
#pragma acc data if(ctsar_get_type(c) == CTSAR_DEV_GPU)
                {
#pragma acc region for independent deviceptr(ex,ey,hz) if(ctsar_get_type(c) == CTSAR_DEV_GPU)
                    for (i = gts; i < gte; i++)
                    {
#pragma acc for independent
                        for (j = 0; j < NY; j++)
                        {
                            hz[((i) * NX) + (j)] -= 0.7*(ex[((i) * NX) + (j+1)] - ex[((i) * NX) + (j)] + ey[((i+1) * NX) + (j)] - ey[((i) * NX) + (j)]);
                            /* hza[i][j] = hza[i][j] - 0.7*(exa[i][j+1] - exa[i][j] + eya[i+1][j] - eya[i][j]); */
                        }
                    }
                }
                ctsar_end(c);
                /* printf("here: %d\n", omp_get_thread_num()); */
            }while(ctsar_loop(c));
        }
    }
}
            /* #pragma acc data if(ctsar_get_type(s) == CTSAR_DEV_GPU)
               {
#pragma acc region for independent copy(ex[gts*(NX+1):(gte+1)*(NX+1)]) copyin(hz[0:NX*NX]) if(ctsar_get_type(s) == CTSAR_DEV_GPU)
for (i = gts; i < gte; i++)
{
for (j = 1; j < NY; j++)
{
ex[((i) * (NX+1)) + (j)] -= 0.5*(hz[((i) * NX) + (j)] - hz[((i) * NX) + (j-1)]);
/\* exa[i][j] = exa[i][j] - 0.5*(hza[i][j] - hza[i][j-1]); *\/
}
}
} */

void runFdtd(DATA_TYPE fict[T_MAX], DATA_TYPE **ex, DATA_TYPE **ey, DATA_TYPE **hz)
{
    int t, i, j;

    for(t=0; t< T_MAX; t++)
    {
#pragma omp parallel
        {
#pragma omp  for private(i,j) schedule(runtime)
            for (j=0; j < NY; j++)
            {
                ey[0][j] = fict[t];
            }

#pragma omp  for private(i,j) schedule(runtime)
            for (i = 1; i < NX; i++)
            {
                for (j = 0; j < NY; j++)
                {
                    ey[i][j] = ey[i][j] - 0.5*(hz[i][j] - hz[i-1][j]);
                    ex[i][j] = ex[i][j] - 0.5*(hz[i][j] - hz[i][j-1]);
                }
            }

            /* #pragma omp parallel for private(i,j)
               for (i = 0; i < NX; i++)
               {
               for (j = 1; j < NY; j++)
               {
               ex[i][j] = ex[i][j] - 0.5*(hz[i][j] - hz[i][j-1]);
               }
               } */

#pragma omp  for private(i,j) schedule(runtime)
            for (i = 0; i < NX; i++)
            {
                for (j = 0; j < NY; j++)
                {
                    hz[i][j] = hz[i][j] - 0.7*(ex[i][j+1] - ex[i][j] + ey[i+1][j] - ey[i][j]);
                }
            }
        }
    }
}


void init_arrays(DATA_TYPE *_fict_, DATA_TYPE **ex, DATA_TYPE **ex_Gpu, DATA_TYPE **ey, DATA_TYPE **ey_Gpu, DATA_TYPE **hz, DATA_TYPE **hz_Gpu)
{
    int i, j;

    for (i = 0; i < T_MAX; i++)
    {
        _fict_[i] = (DATA_TYPE) i;
    }

    for (i = 0; i < NX; i++)
    {
        for (j = 0; j < NY; j++)
        {
            ex[i][j] = ((DATA_TYPE) i*(j+1) + 1) / NX;
            ex_Gpu[i][j] = ((DATA_TYPE) i*(j+1) + 1) / NX;
            ey[i][j] = ((DATA_TYPE) (i-1)*(j+2) + 2) / NX;
            ey_Gpu[i][j] = ((DATA_TYPE) (i-1)*(j+2) + 2) / NX;
            hz[i][j] = ((DATA_TYPE) (i-9)*(j+4) + 3) / NX;
            hz_Gpu[i][j] = ((DATA_TYPE) (i-9)*(j+4) + 3) / NX;
        }
    }
}



void compareResults(DATA_TYPE **hz, DATA_TYPE **hz_outputFromGpu)
{
    int i, j, fail;
    fail = 0;

    for (i=0; i < NX; i++)
    {
        for (j=0; j < NY; j++)
        {
            if (percentDiff(hz[i][j], hz_outputFromGpu[i][j]) > PERCENT_DIFF_ERROR_THRESHOLD)
            {
                fail++;
            }
        }
    }

    // Print results
    printf("Non-Matching CPU-GPU Outputs Beyond Error Threshold of %4.2f Percent: %d\n", PERCENT_DIFF_ERROR_THRESHOLD, fail);

}


int main()
{
    double t_start, t_end;
    ctsar_pre_init();

    DATA_TYPE **ex = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NX,NY+1);
    DATA_TYPE **ex_Gpu = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NX,NY+1);
    DATA_TYPE **ey = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NX+1,NY);
    DATA_TYPE **ey_Gpu = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NX+1,NY);
    DATA_TYPE *fict = (DATA_TYPE *)ctsar_calloc_test(sizeof(DATA_TYPE),T_MAX);
    DATA_TYPE **hz = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NX,NY);
    DATA_TYPE **hz_outputFromGpu = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NX,NY);



    int iters;
    char * env;
    if(env = getenv("TEST_ITERATIONS")){
        iters = atoi(env);
    }else{
        iters = 1;
    }
    t_start = rtclock();

    ctsar_init(&s,NX,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
    ctsar_init(&c,NX,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
    for(int i=0; i < iters; i++)
    {
        double lt_start, lt_end;
        lt_start = rtclock();
        /* grunCorr(symmat_outputFromGpu, symmat_outputFromGpu, stddev_Gpu, mean_Gpu, float_n, eps); */
        if((env = getenv("OMP_CTSAR_FALLBACK")) && atoi(env)){
            runFdtd(fict, ex, ey, hz);
        }else{
            ctsar_init(&s,NX,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
            ctsar_init(&c,NX,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
            grunFdtd(fict, ex_Gpu, ey_Gpu, hz_outputFromGpu);
        }
        lt_end = rtclock();
        fprintf(stderr, "GPU Runtime:it%d: %0.6lfs\n", i, lt_end - lt_start);

    }


    t_end = rtclock();
    fprintf(stderr, "GPU Runtime: %0.6lfs\n", t_end - t_start);


    if((env = getenv("TEST_VERIFY")) && atoi(env)){
        t_start = rtclock();

        runFdtd(fict, ex, ey, hz);

        t_end = rtclock();
        fprintf(stderr, "CPU Runtime: %0.6lf\n", t_end - t_start);

        compareResults(hz, hz_outputFromGpu);
    }

    return 0;
}
