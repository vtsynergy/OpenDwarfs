/**
 * threemm.c: This file is part of the PolyBench/GPU 1.0 test suite.
 *
 *
 * Contact: Scott Grauer-Gray <sgrauerg@gmail.com>
 * Louis-Noel Pouchet <pouchet@cse.ohio-state.edu>
 * Web address: http://www.cse.ohio-state.edu/~pouchet/software/polybench/GPU
 */

#include <unistd.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#include "polybenchUtilFuncts.h"
#include "ctsar.h"
ctsar *s = NULL;

//define the error threshold for the results "not matching"
#define PERCENT_DIFF_ERROR_THRESHOLD 0.05

#define MUL 3
/* Problem size. */
#define NI (512 * MUL)
#define NJ (512 * MUL)
#define NK (512 * MUL)
#define NL (512 * MUL)
#define NM (512 * MUL)

/* Can switch DATA_TYPE between float and double */
typedef float DATA_TYPE;

void threeMMloopg(DATA_TYPE **aa, DATA_TYPE **ba, DATA_TYPE **ca)
{
    ctsar_init(&s,NI,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL);
#pragma omp parallel default(shared)
    {
        int i, j, k;

        /* E := A*B */

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
                for (i = gts; i < gte; i++)
                {
#pragma acc for independent
                    for (j = 0; j < NJ; j++)
                    {
                        DATA_TYPE lc=0.0;

                        for (k = 0; k < NK; ++k)
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

void threeMMloopa(DATA_TYPE **a, DATA_TYPE **b, DATA_TYPE **e)
{
	int i, j, k;

	/* E := A*B */

#pragma omp parallel for private(i,j,k)
	for (i = 0; i < NI; i++)
	{
		for (j = 0; j < NJ; j++)
		{
			e[i][j] = 0;
         
			for (k = 0; k < NK; ++k)
			{
				e[i][j] += a[i][k] * b[k][j];
			}
		}
	}
}

void threeMMloopb(DATA_TYPE **c, DATA_TYPE **d, DATA_TYPE **f)
{
	int i, j, k; 

	/* F := C*D */
#pragma omp parallel for private(i,j,k)
	for (i = 0; i < NJ; i++)
	{
		for (j = 0; j < NL; j++)
		{
			f[i][j] = 0;

			for (k = 0; k < NM; ++k)
			{
				f[i][j] += c[i][k] * d[k][j];
			}
		}
	}
}

void threeMMloopc(DATA_TYPE **e, DATA_TYPE **f, DATA_TYPE **g)
{
	int i, j, k;

	/* G := E*F */
#pragma omp parallel for private(i,j,k)
	for (i = 0; i < NI; i++)
	{      

		for (j = 0; j < NL; j++)
		{
			g[i][j] = 0;
          
			for (k = 0; k < NJ; ++k)
			{
				g[i][j] += e[i][k] * f[k][j];
			}
		}
	}
}


void compareResults(DATA_TYPE **G, DATA_TYPE **G_outputFromGpu)
{
	int i,j,fail;
	fail = 0;

	for (i=0; i < NI; i++)
	{
		for (j=0; j < NL; j++)
		{
			if (percentDiff(G[i][j], G_outputFromGpu[i][j]) > PERCENT_DIFF_ERROR_THRESHOLD)
			{
				fail++;				
			}
		}
	}
	
	// print results
	printf("Non-Matching CPU-GPU Outputs Beyond Error Threshold of %4.2f Percent: %d\n", PERCENT_DIFF_ERROR_THRESHOLD, fail);

}


void init_array(DATA_TYPE **A, DATA_TYPE **B, DATA_TYPE **C, DATA_TYPE **D)
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
			B[i][j] = ((DATA_TYPE) i*(j+1)) / NJ;
		}
	}
  
	for (i = 0; i < NJ; i++)
	{
		for (j = 0; j < NM; j++)
		{
			C[i][j] = ((DATA_TYPE) i*(j+3)) / NL;
		}
	}
  
	for (i = 0; i < NM; i++)
	{
		for (j = 0; j < NL; j++)
		{
			D[i][j] = ((DATA_TYPE) i*(j+2)) / NK;
		}
	}
}

int main(int argc, char** argv)
{
    double t_start, t_end;

    /* Array declaration */
    DATA_TYPE **A = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NI,NK);
    DATA_TYPE **B = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NK,NJ);
    DATA_TYPE **C = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NJ,NM);
    DATA_TYPE **D = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NM,NL);
    DATA_TYPE **E = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NI,NJ);
    DATA_TYPE **E_gpu = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NI,NJ);	
    DATA_TYPE **F = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NJ,NL);
    DATA_TYPE **F_gpu = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NJ,NL);
    DATA_TYPE **G = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NI,NL);
    DATA_TYPE **G_outputFromGpu = (DATA_TYPE **)ctsar_alloc_2d(sizeof(DATA_TYPE),NI,NL);

    /* Initialize array. */
    init_array(A, B, C, D);

    ctsar_pre_init();

    int iters;
    char * env;
    if(env = getenv("TEST_ITERATIONS")){
        iters = atoi(env);
    }else{
        iters = 1;
    }

    t_start = rtclock();
    /* ctsar_init(&ct,NI,CTSAR_START, CTSAR_DEV_CPU | CTSAR_DEV_GPU, NULL,NULL,NULL); */
    for(int i=0; i < iters; i++)
    {
        double lt_start, lt_end;
        lt_start = rtclock();
        /* grunCorr(symmat_outputFromGpu, symmat_outputFromGpu, stddev_Gpu, mean_Gpu, float_n, eps); */
        if((env = getenv("OMP_CTSAR_FALLBACK")) && atoi(env)){
            threeMMloopa(A, B, E_gpu);
            threeMMloopb(C, D, F_gpu);
            threeMMloopc(E_gpu, F_gpu, G_outputFromGpu);
        }else{
            threeMMloopg(A, B, E_gpu);
            threeMMloopg(C, D, F_gpu);
            threeMMloopg(E_gpu, F_gpu, G_outputFromGpu);
        }
        lt_end = rtclock();
        fprintf(stderr, "GPU Runtime:it%d: %0.6lfs\n", i, lt_end - lt_start);
    }
    t_end = rtclock();
    fprintf(stderr, "GPU Runtime: %0.6lfs\n", t_end - t_start);

    if((env = getenv("TEST_VERIFY")) && atoi(env)){
        t_start = rtclock();

        threeMMloopa(A, B, E);
        threeMMloopb(C, D, F);
        threeMMloopc(E, F, G);

        t_end = rtclock();
        fprintf(stderr, "CPU Runtime: %0.6lfs\n", t_end - t_start);

        compareResults(G, G_outputFromGpu);
    }

    return 0;
}
