/**
 * mvt.c: This file is part of the PolyBench/GPU 1.0 test suite.
 *
 *
 * Contact: Scott Grauer-Gray <sgrauerg@gmail.com>
 * Louis-Noel Pouchet <pouchet@cse.ohio-state.edu>
 * Web address: http://www.cse.ohio-state.edu/~pouchet/software/polybench/GPU
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <sys/time.h>

#include "polybenchUtilFuncts.h"
#include "ctsar.h"
ctsar * s = NULL;

//define the error threshold for the results "not matching"
#define PERCENT_DIFF_ERROR_THRESHOLD 0.05

/* Problem size */
#define N 6144 * 2

/* Can switch DATA_TYPE between float and double */
typedef float DATA_TYPE;


void grunMvt(DATA_TYPE **aa, DATA_TYPE *x1a, DATA_TYPE *x2a, DATA_TYPE *y1a, DATA_TYPE *y2a)
{
#pragma omp parallel default(shared)
    {
        int tid = omp_get_thread_num();
        size_t a_pitch = N;

        do{
            ctsar_next(s,N);
            int gts = CSTART(s,tid);
            int gte = CEND(s,tid);
            int i;
            DATA_TYPE *a, *x1,*x2,*y1,*y2;
            /* a = ctsar_reg_mem(s, *aa, sizeof(DATA_TYPE)*N, N,
                    CTSAR_MEM_INPUT); */
            a = ctsar_reg_mem_2d(s, *aa, sizeof(DATA_TYPE), N, N,
                    CTSAR_MEM_INPUT | CTSAR_MEM_PARTIAL | CTSAR_MEM_COLUMN , 0, 0, &a_pitch);
            x1 = ctsar_reg_mem(s, x1a, sizeof(DATA_TYPE)*1, N,
                    CTSAR_MEM_PARTIAL | CTSAR_MEM_INPUT | CTSAR_MEM_OUTPUT);
            x2 = ctsar_reg_mem(s, x2a, sizeof(DATA_TYPE)*1, N,
                    CTSAR_MEM_PARTIAL | CTSAR_MEM_INPUT | CTSAR_MEM_OUTPUT);
            y1 = ctsar_reg_mem(s, y1a, sizeof(DATA_TYPE)*1, N,
                    CTSAR_MEM_INPUT);
            y2 = ctsar_reg_mem(s, y2a, sizeof(DATA_TYPE)*1, N,
                    CTSAR_MEM_INPUT);

            ctsar_start(s);
#pragma acc data region deviceptr(a,x1,x2,y1,y2)  if(ctsar_get_type(s) == CTSAR_DEV_GPU)
            {
#pragma acc region for private(i) independent if(ctsar_get_type(s) == CTSAR_DEV_GPU)
                //copy(x1[0:N],x2[0:N]) copyin(y1[0:N],y2[0:N],a[0:N*N])
                for (i=gts; i<gte; i++)
                {
                    int j;
                    DATA_TYPE lx1 = x1[i], lx2 = x2[i];
                #pragma acc for reduction(+: lx1, lx2)
                    for (j=0; j<N; j++)
                    {
                        lx1 += a[(i * a_pitch) + j] * y1[j];
                        lx2 += a[(j * a_pitch) + i] * y2[j];
                    }
                    x1[i] = lx1;
                    x2[i] = lx2;
                }
            }
            ctsar_end(s);
            /* printf("here: %d\n", omp_get_thread_num()); */
        }while(ctsar_loop(s));
    }
}

void runMvt(DATA_TYPE **a, DATA_TYPE *x1, DATA_TYPE *x2, DATA_TYPE *y1, DATA_TYPE *y2)
{
    int i, j;


#pragma omp parallel for schedule(runtime) private(i,j)
    for (i=0; i<N; i++)
    {
        for (j=0; j<N; j++)
        {
            x1[i] = x1[i] + a[i][j] * y1[j];
        }
        for (j=0; j<N; j++)
        {
            x2[i] = x2[i] + a[j][i] * y2[j];
        }
    }
}

void init_array(DATA_TYPE **A, DATA_TYPE *x1, DATA_TYPE *x1_Gpu, DATA_TYPE *x2, DATA_TYPE *x2_Gpu, DATA_TYPE *y1, DATA_TYPE *y2)
{
    int i, j;

    for (i = 0; i < N; i++)
    {
        x1[i] = ((DATA_TYPE) i) / N;
        x1_Gpu[i] = ((DATA_TYPE) i) / N;
        x2[i] = ((DATA_TYPE) i + 1) / N;
        x2_Gpu[i] = ((DATA_TYPE) i + 1) / N;
        y1[i] = ((DATA_TYPE) i + 3) / N;
        y2[i] = ((DATA_TYPE) i + 4) / N;
        for (j = 0; j < N; j++)
        {
            A[i][j] = ((DATA_TYPE) i*j) / N;
        }
    }
}


void compareResults(DATA_TYPE *x1, DATA_TYPE *x1_outputFromGpu, DATA_TYPE *x2, DATA_TYPE *x2_outputFromGpu)
{
    int i, fail;
    fail = 0;

    for (i=0; i<N; i++)
    {
        if (percentDiff(x1[i], x1_outputFromGpu[i]) > PERCENT_DIFF_ERROR_THRESHOLD)
        {
            fail++;
            if(fail < 10)
                printf("x1:non-matching at i=%d cpu=%f gpu=%f\n", i, x1[i], x1_outputFromGpu[i]);
        }

        if (percentDiff(x2[i], x2_outputFromGpu[i]) > PERCENT_DIFF_ERROR_THRESHOLD)
        {
            fail++;
            if(fail < 10)
                printf("x2:non-matching at i=%d cpu=%f gpu=%f\n", i, x2[i], x2_outputFromGpu[i]);
        }
    }

    // Print results
    printf("Non-Matching CPU-GPU Outputs Beyond Error Threshold of %4.2f Percent: %d\n", PERCENT_DIFF_ERROR_THRESHOLD, fail);

}



int main()
{
    double t_start, t_end;

    DATA_TYPE **a = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),N,N);
    DATA_TYPE *x1 = (DATA_TYPE *)calloc(sizeof(DATA_TYPE),N);
    DATA_TYPE *x1_outputFromGpu = (DATA_TYPE *)calloc(sizeof(DATA_TYPE),N);
    DATA_TYPE *x2 = (DATA_TYPE *)calloc(sizeof(DATA_TYPE),N);
    DATA_TYPE *x2_outputFromGpu = (DATA_TYPE *)calloc(sizeof(DATA_TYPE),N);
    DATA_TYPE *y1 = (DATA_TYPE *)calloc(sizeof(DATA_TYPE),N);
    DATA_TYPE *y2 = (DATA_TYPE *)calloc(sizeof(DATA_TYPE),N);


    //initialize the arrays for running on the CPU and GPU
    init_array(a, x1, x1_outputFromGpu, x2, x2_outputFromGpu, y1, y2);

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

    //run the algorithm on the GPU
    /* grunMvt(a, x1_outputFromGpu, x2_outputFromGpu, y1, y2); // parameters are initialized in decls.h and are initialized with init_array() */

    for(int i=0; i < iters; i++)
    {
        double lt_start, lt_end;
        lt_start = rtclock();
        /* grunCorr(symmat_outputFromGpu, symmat_outputFromGpu, stddev_Gpu, mean_Gpu, float_n, eps); */
        if((env = getenv("OMP_CTSAR_FALLBACK")) && atoi(env)){
            runMvt(a, x1_outputFromGpu, x2_outputFromGpu, y1, y2); // parameters are initialized in decls.h and are initialized with init_array()
        }else{
            ctsar_init(&s,N,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
            grunMvt(a, x1_outputFromGpu, x2_outputFromGpu, y1, y2); // parameters are initialized in decls.h and are initialized with init_array()
        }
        lt_end = rtclock();
        fprintf(stderr, "GPU Runtime:it%d: %0.6lfs\n", i, lt_end - lt_start);

    }

    t_end = rtclock();
    fprintf(stderr, "GPU Runtime: %0.6lf\n", t_end - t_start);

    if((env = getenv("TEST_VERIFY")) && atoi(env)){
        t_start = rtclock();

        //run the algorithm on the CPU
        runMvt(a, x1, x2, y1, y2);

        t_end = rtclock();
        fprintf(stderr, "CPU Runtime: %0.6lf\n", t_end - t_start);

        compareResults(x1, x1_outputFromGpu, x2, x2_outputFromGpu);
    }

    return 0;
}
