/**
 * bicg.c: This file is part of the PolyBench/GPU 1.0 test suite.
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
ctsar * ct = NULL;

//define the error threshold for the results "not matching"
#define PERCENT_DIFF_ERROR_THRESHOLD 0.05

/* Problem size */
#define NX 4096 * 5
#define NY 4096 * 5

/* Can switch DATA_TYPE between float and double */
typedef float DATA_TYPE;


void grunBicg(DATA_TYPE **aa, DATA_TYPE *pa, DATA_TYPE *qa, DATA_TYPE *ra, DATA_TYPE *sa)
{
    int i, j, k, l;
    static DATA_TYPE * a_far = NULL;
#pragma omp parallel default(shared) private(j,i,k)
    {
        int tid = omp_get_thread_num();
        size_t a_pitch = NX;
        if(a_far == NULL){
#pragma omp barrier
            if(tid == 4)
                a_far = (DATA_TYPE *) malloc(NX*NY*sizeof(DATA_TYPE));
#pragma omp barrier
        }
        DATA_TYPE * restrict a, * restrict p, * restrict q, * restrict r, * restrict s;
        a = ctsar_reg_mem_2d(ct, aa[0], sizeof(DATA_TYPE), NX, NY,
                CTSAR_MEM_INPUT | CTSAR_MEM_PARTIAL | CTSAR_MEM_COLUMN, 0, 0, &a_pitch);
        p = ctsar_reg_mem(ct, pa, sizeof(DATA_TYPE), NX,
                CTSAR_MEM_INPUT);
        r = ctsar_reg_mem(ct, ra, sizeof(DATA_TYPE), NX,
                CTSAR_MEM_INPUT);

        q = ctsar_reg_mem(ct, qa, sizeof(DATA_TYPE), NX,
                CTSAR_MEM_INPUT | CTSAR_MEM_OUTPUT | CTSAR_MEM_PARTIAL);
        s = ctsar_reg_mem(ct, sa, sizeof(DATA_TYPE), NX,
                CTSAR_MEM_INPUT | CTSAR_MEM_OUTPUT | CTSAR_MEM_PARTIAL);

        do{
            ctsar_next(ct,NX);
            int gts = CSTART(ct,tid);
            int gte = CEND(ct,tid);

            ctsar_start(ct);
            /* if(tid >3){
                fprintf(stderr,"tid:%d start:%d, end:%d\n", tid, gts, gte);
                for(i=0; i<gts; i++){
                    memcpy(a_far+(NX*i)+gts, (*aa)+(NX*i)+gts, sizeof(DATA_TYPE) * (gte-gts));
                }
                memcpy(a_far+(NX*gts), (*aa)+(NX*gts), NX * sizeof(DATA_TYPE)* (gte-gts));
                for(i=gte; i<NX; i++){
                    memcpy(a_far+(NX*i)+gts, (*aa)+(NX*i)+gts, sizeof(DATA_TYPE) * (gte-gts));
                }
                a = a_far;
            } */
#pragma acc data deviceptr(s,q,r,a,p) if(ctsar_get_type(ct) == CTSAR_DEV_GPU)
            {
#pragma acc region for independent private(i,j) if(ctsar_get_type(ct) == CTSAR_DEV_GPU)
                for (i = gts; i < gte; i++)
                {
                    DATA_TYPE ls = s[i], lq = q[i];
                    #pragma acc for reduction(+:ls, lq)
                    for (j = 0; j < NY; j++)
                    {
                        ls += r[j] * a[((j) * a_pitch) + (i)];
                        lq += a[((i) * a_pitch) + (j)] * p[j];
                    }
                    s[i] = ls;
                    q[i] = lq;
                }
            }
            ctsar_end(ct);
            /* printf("here: %d\n", omp_get_thread_num()); */
        }while(ctsar_loop(ct));
    }
    /* free(a_far); */
}

void runBicg(DATA_TYPE **a, DATA_TYPE *p, DATA_TYPE *q, DATA_TYPE *r, DATA_TYPE *s)
{
    int i, j, k, l;

#pragma omp parallel for schedule(runtime) private(i,j)
    for (i = 0; i < NX; i++)
    {
        for (j = 0; j < NY; j++)
        {
            s[i] = s[i] + r[j] * a[j][i];
            q[i] = q[i] + a[i][j] * p[j];
        }
    }
}


void init_array(DATA_TYPE **A, DATA_TYPE *p, DATA_TYPE *r)
{
    int i, j;

    for (i = 0; i < NX; i++)
    {
        r[i] = i * M_PI;

        for (j = 0; j < NY; j++)
        {
            A[i][j] = ((DATA_TYPE) i*j) / NX;
        }
    }

    for (i = 0; i < NY; i++)
    {
        p[i] = i * M_PI;
    }
}


void compareResults(DATA_TYPE *s, DATA_TYPE *s_outputFromGpu, DATA_TYPE *q, DATA_TYPE *q_outputFromGpu)
{
    int i,fail;
    fail = 0;

    // Compare s with s_cuda
    for (i=0; i<NX; i++)
    {
        if (percentDiff(q[i], q_outputFromGpu[i]) > PERCENT_DIFF_ERROR_THRESHOLD)
        {
            fail++;
        }
    }

    for (i=0; i<NY; i++)
    {
        if (percentDiff(s[i], s_outputFromGpu[i]) > PERCENT_DIFF_ERROR_THRESHOLD)
        {
            fail++;
        }		
    }

    // print results
    printf("Non-Matching CPU-GPU Outputs Beyond Error Threshold of %4.2f Percent: %d\n", PERCENT_DIFF_ERROR_THRESHOLD, fail);

}


int main(int argc, char** argv)
{
    double t_start, t_end;

    /* Array declaration */
    DATA_TYPE **A = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NX,NY);
    DATA_TYPE *p = (DATA_TYPE *)ctsar_calloc_test(sizeof(DATA_TYPE),NY);
    DATA_TYPE *q = (DATA_TYPE *)ctsar_calloc_test(sizeof(DATA_TYPE),NX);
    DATA_TYPE *q_outputFromGpu = (DATA_TYPE *)ctsar_calloc_test(sizeof(DATA_TYPE),NX);
    DATA_TYPE *r = (DATA_TYPE *)ctsar_calloc_test(sizeof(DATA_TYPE),NX);
    DATA_TYPE *s = (DATA_TYPE *)ctsar_calloc_test(sizeof(DATA_TYPE),NY);
    DATA_TYPE *s_outputFromGpu = (DATA_TYPE *)ctsar_calloc_test(sizeof(DATA_TYPE),NY);

    /* Initialize array. */
    init_array(A, p, r);

    ctsar_pre_init();


    int iters;
    char * env;
    if(env = getenv("TEST_ITERATIONS")){
        iters = atoi(env);
    }else{
        iters = 1;
    }

    t_start = rtclock();
    ctsar_init(&ct,NX,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);

    for(int i=0; i < iters; i++)
    {
        double lt_start, lt_end;
        lt_start = rtclock();
        /* grunCorr(symmat_outputFromGpu, symmat_outputFromGpu, stddev_Gpu, mean_Gpu, float_n, eps); */
        if((env = getenv("OMP_CTSAR_FALLBACK")) && atoi(env)){
            runBicg(A, p, q_outputFromGpu, r, s_outputFromGpu);
        }else{
            ctsar_init(&ct, NX, CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);

            grunBicg(A, p, q_outputFromGpu, r, s_outputFromGpu);
        }
        lt_end = rtclock();
        fprintf(stderr, "GPU Runtime:it%d: %0.6lfs\n", i, lt_end - lt_start);

    }

    t_end = rtclock();
    fprintf(stderr, "GPU Runtime: %0.6lfs\n", t_end - t_start);



    if((env = getenv("TEST_VERIFY")) && atoi(env)){
        t_start = rtclock();

        runBicg(A, p, q, r, s);

        t_end = rtclock();
        fprintf(stderr, "CPU Runtime: %0.6lfs\n", t_end - t_start);

        compareResults(s, s_outputFromGpu, q, q_outputFromGpu);
    }
    return 0;
}
