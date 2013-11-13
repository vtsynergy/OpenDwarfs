#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>
#include <getopt.h>

#include "../../include/rdtsc.h"
#include "../../include/common_args.h"

/*
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif
*/

/*
#define CHKERR(err, str) \
    if ( err != CL_SUCCESS) \
    {\
        fprintf(stderr, "Error in executing \"%s\", %d\n", str, err);  \
        exit(EXIT_FAILURE);\
    }
*/
#define CHECK_NULL_ERROR(err, str) \
    if ( err == NULL) \
    { \
        fprintf(stderr, "Error creating objects in \"%s\"\n", str); \
        exit(EXIT_FAILURE);\
    }

#define T           1000        /* Number of static observations */
#define S           2           /* Number of static symbols */
#define N           60          /* Number of static states */
#define ITERATIONS  1           /* Number of iterations */

#define MAX_THREADS_PER_BLOCK 256
#define BLOCK_DIM 16
#define EXIT_ERROR 1

#define SDOT_BLOCK_SIZE (128)
#define SDOT_BLOCK_NUM (80)

#define MVMUL_BLOCK_SIZE (128)
#define MVMUL_BLOCK_NUM  (64)

cl_platform_id           _cl_firstPlatform;
cl_device_id             _cl_firstDevice;
//cl_context               _cl_context;
//cl_command_queue         _cl_commandQueue;
cl_program               _cl_program;

cl_kernel _cl_kernel_init_ones_dev;
cl_kernel _cl_kernel_init_alpha_dev;
cl_kernel _cl_kernel_calc_alpha_dev;
cl_kernel _cl_kernel_scale_alpha_dev;
cl_kernel _cl_kernel_init_beta_dev;
cl_kernel _cl_kernel_calc_beta_dev;
cl_kernel _cl_kernel_calc_gamma_dev;
cl_kernel _cl_kernel_calc_xi_dev;
cl_kernel _cl_kernel_est_a_dev;
cl_kernel _cl_kernel_scale_a_dev;
cl_kernel _cl_kernel_acc_b_dev;
cl_kernel _cl_kernel_est_b_dev;
cl_kernel _cl_kernel_scale_b_dev;
cl_kernel _cl_kernel_est_pi_dev;
cl_kernel _cl_kernel_s_dot_kernel_naive;
cl_kernel _cl_kernel_sgemvt_kernel_naive;
cl_kernel _cl_kernel_sgemvn_kernel_naive;

size_t LoadProgramSource(const char *filename, const char **progSrc) 
{
    FILE *f = fopen(filename, "r");
    fseek(f, 0, SEEK_END);
    size_t len = (size_t) ftell(f);
    *progSrc = (const char *) malloc(sizeof(char)*len);
    rewind(f);
    fread((void *) *progSrc, len, 1, f);
    fclose(f);
    return len;
}

/* Simple init the openCL platform */
void init_cl() 
{
    cl_int errNum;
    size_t progLen;
    const char *progSrc;
    
    ocd_initCL();//KK
    
    /* identify the first platform */
    //errNum = clGetPlatformIDs(1, &_cl_firstPlatform, NULL);
    //CHKERR(errNum , "finding any available OpenCL platform!");
    
    /* identify the first accessible device associated with the platform */
    //errNum = clGetDeviceIDs(_cl_firstPlatform, CL_DEVICE_TYPE_GPU, 1, &_cl_firstDevice, NULL);
    //CHKERR(errNum, "finding any available OpenCL device!");
    
    /* create a context containing only one device which is created earlier */
    //context = clCreateContext(NULL, 1, &_cl_firstDevice, NULL, NULL, &errNum);
    //CHKERR(errNum, "creating GPU context!");

    /* create command queue to target the device, which can support profiling */
    //_cl_commandQueue = clCreateCommandQueue(context, _cl_firstDevice, CL_QUEUE_PROFILING_ENABLE, NULL);
    //CHECK_NULL_ERROR(_cl_commandQueue, "_cl_commandQueue");
    
    /* read the source code into the char array */
    progLen = LoadProgramSource("bwa_hmm_opencl.cl", &progSrc);

    /* create program from the source code */
    _cl_program = clCreateProgramWithSource(context, 1, &progSrc, &progLen, NULL);
    CHECK_NULL_ERROR(_cl_program, "_cl_program");
    free((void *)progSrc);

    /* compile the source code for the device in the context */
    errNum = clBuildProgram(_cl_program, 1, &device_id, "-I .", NULL, NULL);
    if(errNum != CL_SUCCESS)
    {
        char buildLog[16384];
        clGetProgramBuildInfo(_cl_program, device_id, CL_PROGRAM_BUILD_LOG, sizeof(buildLog), buildLog, NULL);
        fprintf(stderr, "Error in compiling the kernels: \n");
        fprintf(stderr, "%s\n", buildLog);
    }
	
    /* create the kernels from the program */
	_cl_kernel_init_ones_dev   = clCreateKernel(_cl_program, "init_ones_dev", NULL);
	_cl_kernel_init_alpha_dev  = clCreateKernel(_cl_program, "init_alpha_dev", NULL);
	_cl_kernel_calc_alpha_dev  = clCreateKernel(_cl_program, "calc_alpha_dev", NULL);
	_cl_kernel_scale_alpha_dev = clCreateKernel(_cl_program, "scale_alpha_dev", NULL);
	_cl_kernel_init_beta_dev   = clCreateKernel(_cl_program, "init_beta_dev", NULL);
	_cl_kernel_calc_beta_dev   = clCreateKernel(_cl_program, "calc_beta_dev", NULL);
	_cl_kernel_calc_gamma_dev  = clCreateKernel(_cl_program, "calc_gamma_dev", NULL);
	_cl_kernel_calc_xi_dev     = clCreateKernel(_cl_program, "calc_xi_dev", NULL);
	_cl_kernel_est_a_dev       = clCreateKernel(_cl_program, "est_a_dev", NULL);
	_cl_kernel_scale_a_dev     = clCreateKernel(_cl_program, "scale_a_dev", NULL);
	_cl_kernel_acc_b_dev       = clCreateKernel(_cl_program, "acc_b_dev", NULL);
	_cl_kernel_est_b_dev       = clCreateKernel(_cl_program, "est_b_dev", NULL);
	_cl_kernel_scale_b_dev     = clCreateKernel(_cl_program, "scale_b_dev", NULL);
	_cl_kernel_est_pi_dev      = clCreateKernel(_cl_program, "est_pi_dev", NULL);
	_cl_kernel_s_dot_kernel_naive    = clCreateKernel(_cl_program, "s_dot_kernel_naive", NULL);
	_cl_kernel_sgemvt_kernel_naive   = clCreateKernel(_cl_program, "mvm_trans_kernel_naive", NULL);
	_cl_kernel_sgemvn_kernel_naive   = clCreateKernel(_cl_program, "mvm_non_kernel_naive", NULL);
	
    CHECK_NULL_ERROR( _cl_kernel_init_ones_dev, "_cl_kernel_init_ones_dev");
    CHECK_NULL_ERROR( _cl_kernel_init_alpha_dev, "_cl_kernel_init_alpha_dev");
	CHECK_NULL_ERROR( _cl_kernel_calc_alpha_dev, "_cl_kernel_calc_alpha_dev");
	CHECK_NULL_ERROR( _cl_kernel_scale_alpha_dev, "_cl_kernel_scale_alpha_dev");
	CHECK_NULL_ERROR( _cl_kernel_init_beta_dev, "_cl_kernel_init_beta_dev");
	CHECK_NULL_ERROR( _cl_kernel_calc_beta_dev, "_cl_kernel_calc_beta_dev");
    CHECK_NULL_ERROR( _cl_kernel_calc_gamma_dev, "_cl_kernel_calc_gamma_dev");
	CHECK_NULL_ERROR( _cl_kernel_calc_xi_dev, "_cl_kernel_calc_xi_dev");
	CHECK_NULL_ERROR( _cl_kernel_est_a_dev, "_cl_kernel_est_a_dev");
	CHECK_NULL_ERROR( _cl_kernel_scale_a_dev, "_cl_kernel_scale_a_dev");
	CHECK_NULL_ERROR( _cl_kernel_acc_b_dev, "_cl_kernel_acc_b_dev");
    CHECK_NULL_ERROR( _cl_kernel_est_b_dev, "_cl_kernel_est_b_dev");
	CHECK_NULL_ERROR( _cl_kernel_scale_b_dev, "_cl_kernel_scale_b_dev");
	CHECK_NULL_ERROR( _cl_kernel_est_pi_dev, "_cl_kernel_est_pi_dev");
	CHECK_NULL_ERROR( _cl_kernel_s_dot_kernel_naive, "_cl_kernel_s_dot_kernel_naive");
	CHECK_NULL_ERROR( _cl_kernel_sgemvt_kernel_naive, "_cl_kernel_sgemvt_kernel_naive");
    CHECK_NULL_ERROR( _cl_kernel_sgemvn_kernel_naive, "_cl_kernel_sgemvn_kernel_naive");

}

static int imax(int x, int y)
{
    return (x > y) ? x: y;
}

/* Subtracts time values to determine run time */
int timeval_subtract(struct timeval *result, struct timeval *t2, 
                                            struct timeval *t1)
{
    long int diff = (t2->tv_usec + 1000000 * t2->tv_sec) - 
                    (t1->tv_usec + 1000000 * t1->tv_sec);
    result->tv_sec = diff / 1000000;
    result->tv_usec = diff % 1000000;

    return (diff<0);
}


/* Starts timer */
void tic(struct timeval *timer)
{
    gettimeofday(timer, NULL);
}

/* Stops timer and prints difference to the screen */
void toc(struct timeval *timer)
{
    struct timeval tv_end, tv_diff;
    
    gettimeofday(&tv_end, NULL);
    timeval_subtract(&tv_diff, &tv_end, timer);
    printf("%ld.%06ld\t", tv_diff.tv_sec, tv_diff.tv_usec);
}

/* individual paramters of a Hidden Markov Model */
typedef struct {
    int nstates;            /**< number of states in the HMM */
    int nsymbols;           /**< number of possible symbols */
    float *a;               /**< A matrix - state transition probabilities */
    float *b;               /**< B matrix - symbol output probabilities */
    float *pi;              /**< Pi matrix - initial state probabilities */
} Hmm;


/* the observation sequence and length */
typedef struct {
    int length;             /**< the length of the observation sequence */
    int *data;              /**< the observation sequence */
} Obs;

/* Free the memory associated with the HMM and observation sequence */
void free_vars(Hmm *hmm, Obs *obs)
{
    if (hmm != NULL) {
        if (hmm->a != NULL) {
            free(hmm->a);
        }
        if (hmm->b != NULL) {
            free(hmm->b);
        }
        if (hmm->pi != NULL) {
            free(hmm->pi);
        }
        free(hmm);
    }
    if (obs != NULL) {
        if (obs->data != NULL) {
            free(obs->data);
        }
        free(obs);
    }
}

/* Global variables for device */
int nstates;                        /* The number of states in the HMM */
int nsymbols;                       /* The number of possible symbols */
int *obs;                           /* The observation sequence */
int length;                         /* The length of the observation sequence */
float *scale;                       /* Scaling factor as determined by alpha */
cl_mem a_d;                         /* A matrix on GPU */
cl_mem b_d;                         /* B matrix on GPU */
cl_mem pi_d;                        /* Pi matrix on GPU */
cl_mem alpha_d;                     /* Forward variables (alpha) on GPU */
cl_mem beta_d;                      /* Backward variables (beta) on GPU */
cl_mem gamma_sum_d;                 /* Sum of gamma variables on GPU */
cl_mem xi_sum_d;                    /* Sum of xi variables on GPU */
cl_mem c_d;                         /* Temporary array on GPU */
cl_mem ones_n_d;                    /* Length <states> array of 1s on GPU */
cl_mem ones_s_d;                    /* Length <symbols> array of 1s on GPU */

/************************************************************************/
/* This function is used to compute the dot product of two single-      */
/* precision vectors. Add offset support to cl_mem type files           */
/* OUTPUTS:                                                             */
/*   returns single-precision dot product                               */
/*                                                                      */
/************************************************************************/
float dot_production(int n, cl_mem paramA, int offsetA, cl_mem paramB, int offsetB)
{
    float result = 0.0f;
    int i;
    cl_int errNum;
 
    int blocks, threads;
    if(n < SDOT_BLOCK_NUM)
    {
        blocks = n;
    } else 
    {
        blocks = SDOT_BLOCK_NUM ;
    }
    threads =SDOT_BLOCK_SIZE ;

    if(n <= 0)
    {
        return result;
    }

    /* for the sum reduction on CPU */
    float *partialSum = (float*)calloc(n, sizeof(float));
    
	cl_mem partialSum_d;
	partialSum_d = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * n, NULL, NULL);
    CHECK_NULL_ERROR(commands, "partialSum_d");

    /* set partialSum_d to all zeros */
  	errNum = clEnqueueWriteBuffer(commands, partialSum_d, CL_TRUE, 0,  sizeof(float) * n, partialSum, 0, NULL, &ocdTempEvent);
  	clFinish(commands);
    START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "partialSum_d Data Copy", ocdTempTimer)
    END_TIMER(ocdTempTimer)
	CHKERR(errNum, "writing buffer partialSum_d");
	
	errNum  = clSetKernelArg(_cl_kernel_s_dot_kernel_naive, 0, sizeof(int), &n);
	errNum |= clSetKernelArg(_cl_kernel_s_dot_kernel_naive, 1, sizeof(cl_mem), &paramA);
	errNum |= clSetKernelArg(_cl_kernel_s_dot_kernel_naive, 2, sizeof(int), &offsetA);
	errNum |= clSetKernelArg(_cl_kernel_s_dot_kernel_naive, 3, sizeof(cl_mem), &paramB);
	errNum |= clSetKernelArg(_cl_kernel_s_dot_kernel_naive, 4, sizeof(int), &offsetB);
	errNum |= clSetKernelArg(_cl_kernel_s_dot_kernel_naive, 5, sizeof(cl_mem), &partialSum_d);
	CHKERR(errNum, "setting kernel _cl_kernel_s_dot_kernel_naive arguments");

	size_t globalWorkSize[1] = { blocks * threads };
	size_t localWorkSize[1]  = { threads };
	
	errNum = clEnqueueNDRangeKernel(commands, _cl_kernel_s_dot_kernel_naive, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
	clFinish(commands);
    START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "_cl_kernel_s_dot_kernel_naive Kernel", ocdTempTimer)
    END_TIMER(ocdTempTimer)
	CHKERR(errNum, "queuing kernel _cl_kernel_s_dot_kernel_naive");
	
	errNum = clEnqueueReadBuffer(commands, partialSum_d, CL_TRUE, 0,  sizeof(float) * n, partialSum, 0, NULL, &ocdTempEvent);
	clFinish(commands);
    START_TIMER(ocdTempEvent, OCD_TIMER_D2H, "partialSum_d Data Copy", ocdTempTimer)
    END_TIMER(ocdTempTimer)
	CHKERR(errNum, "reading buffer partialSum_d");
    
    for(i = 0; i < n; i++)
    {
        result += partialSum[i];
    }

    return result;
}

/************************************************************************/
/* This function is used for matrix and vector multiplication.          */
/* Add transposed matrix support                                        */
/************************************************************************/
void mat_vec_mul(char trans, int m, int n, cl_mem A, int lda, cl_mem x, int offsetX, cl_mem y, int offsetY)
{
    cl_int errNum;

    if((trans != 'n') && (trans != 't'))
    {
        return;
    }
	
	if(trans == 't')
    {
        /*
        errNum  = clSetKernelArg(_cl_kernel_sgemvt_kernel, 0, sizeof(int), &m);
		errNum |= clSetKernelArg(_cl_kernel_sgemvt_kernel, 1, sizeof(int), &n);
		errNum |= clSetKernelArg(_cl_kernel_sgemvt_kernel, 2, sizeof(cl_mem), &A);
		errNum |= clSetKernelArg(_cl_kernel_sgemvt_kernel, 3, sizeof(int), &lda);
		errNum |= clSetKernelArg(_cl_kernel_sgemvt_kernel, 4, sizeof(cl_mem), &x);
		errNum |= clSetKernelArg(_cl_kernel_sgemvt_kernel, 5, sizeof(int), &offsetX);
		errNum |= clSetKernelArg(_cl_kernel_sgemvt_kernel, 6, sizeof(cl_mem), &y);
		errNum |= clSetKernelArg(_cl_kernel_sgemvt_kernel, 7, sizeof(int), &offsetY);
		CHKERR(errNum, "setting kernel _cl_kernel_sgemvt_kernel arguments");
		
		size_t globalWorkSize[1] = { MVMUL_BLOCK_NUM * MVMUL_BLOCK_SIZE };
		size_t localWorkSize[1]  = { MVMUL_BLOCK_SIZE };
		
		errNum = clEnqueueNDRangeKernel(_cl_commandQueue, _cl_kernel_sgemvt_kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
		CHKERR(errNum, "queuing kernel _cl_kernel_sgemvt_kernel for execution");
        */
        errNum  = clSetKernelArg(_cl_kernel_sgemvt_kernel_naive, 0, sizeof(int), &m);
		errNum |= clSetKernelArg(_cl_kernel_sgemvt_kernel_naive, 1, sizeof(int), &n);
		errNum |= clSetKernelArg(_cl_kernel_sgemvt_kernel_naive, 2, sizeof(cl_mem), &A);
		errNum |= clSetKernelArg(_cl_kernel_sgemvt_kernel_naive, 3, sizeof(int), &lda);
		errNum |= clSetKernelArg(_cl_kernel_sgemvt_kernel_naive, 4, sizeof(cl_mem), &x);
		errNum |= clSetKernelArg(_cl_kernel_sgemvt_kernel_naive, 5, sizeof(int), &offsetX);
		errNum |= clSetKernelArg(_cl_kernel_sgemvt_kernel_naive, 6, sizeof(cl_mem), &y);
		errNum |= clSetKernelArg(_cl_kernel_sgemvt_kernel_naive, 7, sizeof(int), &offsetY);
		CHKERR(errNum, "setting kernel _cl_kernel_sgemvt_kernel_naive arguments");
		
		size_t globalWorkSize[1] = { MVMUL_BLOCK_NUM * MVMUL_BLOCK_SIZE };
		size_t localWorkSize[1]  = { MVMUL_BLOCK_SIZE };
		
		errNum = clEnqueueNDRangeKernel(commands, _cl_kernel_sgemvt_kernel_naive, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
		clFinish(commands);
    	START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "_cl_kernel_sgemvt_kernel_naive Kernel", ocdTempTimer)
    	END_TIMER(ocdTempTimer)
		CHKERR(errNum, "queuing kernel _cl_kernel_sgemvt_kernel_naive for execution");
    } else
	{
        /*
		errNum  = clSetKernelArg(_cl_kernel_sgemvn_kernel, 0, sizeof(int), &m);
		errNum |= clSetKernelArg(_cl_kernel_sgemvn_kernel, 1, sizeof(int), &n);
		errNum |= clSetKernelArg(_cl_kernel_sgemvn_kernel, 2, sizeof(cl_mem), &A);
		errNum |= clSetKernelArg(_cl_kernel_sgemvn_kernel, 3, sizeof(int), &lda);
		errNum |= clSetKernelArg(_cl_kernel_sgemvn_kernel, 4, sizeof(cl_mem), &x);
		errNum |= clSetKernelArg(_cl_kernel_sgemvn_kernel, 5, sizeof(int), &offsetX);
		errNum |= clSetKernelArg(_cl_kernel_sgemvn_kernel, 6, sizeof(cl_mem), &y);
		errNum |= clSetKernelArg(_cl_kernel_sgemvn_kernel, 7, sizeof(int), &offsetY);
		CHKERR(errNum, "setting kernel _cl_kernel_sgemvn_kernel arguments");
		
		size_t globalWorkSize[1] = { MVMUL_BLOCK_NUM * MVMUL_BLOCK_SIZE };
		size_t localWorkSize[1]  = { MVMUL_BLOCK_SIZE };
		
		errNum = clEnqueueNDRangeKernel(_cl_commandQueue, _cl_kernel_sgemvn_kernel, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, NULL);
		CHKERR(errNum, "queuing kernel _cl_kernel_sgemvn_kernel for execution");
        */
        errNum  = clSetKernelArg(_cl_kernel_sgemvn_kernel_naive, 0, sizeof(int), &m);
		errNum |= clSetKernelArg(_cl_kernel_sgemvn_kernel_naive, 1, sizeof(int), &n);
		errNum |= clSetKernelArg(_cl_kernel_sgemvn_kernel_naive, 2, sizeof(cl_mem), &A);
		errNum |= clSetKernelArg(_cl_kernel_sgemvn_kernel_naive, 3, sizeof(int), &lda);
		errNum |= clSetKernelArg(_cl_kernel_sgemvn_kernel_naive, 4, sizeof(cl_mem), &x);
		errNum |= clSetKernelArg(_cl_kernel_sgemvn_kernel_naive, 5, sizeof(int), &offsetX);
		errNum |= clSetKernelArg(_cl_kernel_sgemvn_kernel_naive, 6, sizeof(cl_mem), &y);
		errNum |= clSetKernelArg(_cl_kernel_sgemvn_kernel_naive, 7, sizeof(int), &offsetY);
		CHKERR(errNum, "setting kernel _cl_kernel_sgemvn_kernel_naive arguments");
		
		size_t globalWorkSize[1] = { MVMUL_BLOCK_NUM * MVMUL_BLOCK_SIZE };
		size_t localWorkSize[1]  = { MVMUL_BLOCK_SIZE };
		
		errNum = clEnqueueNDRangeKernel(commands, _cl_kernel_sgemvn_kernel_naive, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
		clFinish(commands);
    	START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "_cl_kernel_sgemvn_kernel_naive Kernel", ocdTempTimer)
    	END_TIMER(ocdTempTimer)
		CHKERR(errNum, "queuing kernel _cl_kernel_sgemvn_kernel_naive for execution");
	}

}

/*******************************************************************************
 * Supporting functions
 */
 
/* Calculates the forward variables (alpha) for an HMM and obs. sequence */
float calc_alpha()
{

    int threads_per_block;
    int nblocks;
    int offset_cur;
    int offset_prev;
    float log_lik;
    int t;
    cl_int errNum;

    /* Initialize alpha variables */
    threads_per_block = MAX_THREADS_PER_BLOCK;
    nblocks = (nstates + threads_per_block - 1) / threads_per_block;
	size_t globalWorkSize[1] = { nblocks * threads_per_block };
	size_t localWorkSize[1]  = { threads_per_block };
	
	// init_alpha_dev<<<nblocks, threads_per_block>>>( b_d, 
    //                                                 pi_d, 
    //                                                 nstates, 
    //                                                 alpha_d,
    //                                                 ones_n_d,
    //                                                 obs[0]);
	errNum  = clSetKernelArg(_cl_kernel_init_alpha_dev, 0, sizeof(cl_mem), &b_d);
	errNum |= clSetKernelArg(_cl_kernel_init_alpha_dev, 1, sizeof(cl_mem), &pi_d);
	errNum |= clSetKernelArg(_cl_kernel_init_alpha_dev, 2, sizeof(int), &nstates);
	errNum |= clSetKernelArg(_cl_kernel_init_alpha_dev, 3, sizeof(cl_mem), &alpha_d);
	errNum |= clSetKernelArg(_cl_kernel_init_alpha_dev, 4, sizeof(cl_mem), &ones_n_d);
	errNum |= clSetKernelArg(_cl_kernel_init_alpha_dev, 5, sizeof(int), obs);
	CHKERR(errNum, "setting kernel _cl_kernel_init_alpha_dev arguments");
	
	errNum = clEnqueueNDRangeKernel(commands, _cl_kernel_init_alpha_dev, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
	clFinish(commands);
    START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "_cl_kernel_init_alpha_dev Kernel", ocdTempTimer)
    END_TIMER(ocdTempTimer)
	CHKERR(errNum, "queuing kernel _cl_kernel_init_alpha_dev for execution");
    
    /* Sum alpha values to get scaling factor */
    scale[0] = dot_production(nstates, alpha_d, 0, ones_n_d, 0);
    
    /* Scale alpha values */
    // scale_alpha_dev<<<nblocks, threads_per_block>>>(    nstates, 
    //                                                     alpha_d, 
    //                                                     scale[0]);
    int tmp = 0;
	errNum  = clSetKernelArg(_cl_kernel_scale_alpha_dev, 0, sizeof(int), &nstates);
	errNum |= clSetKernelArg(_cl_kernel_scale_alpha_dev, 1, sizeof(cl_mem), &alpha_d);
	errNum |= clSetKernelArg(_cl_kernel_scale_alpha_dev, 2, sizeof(int), &tmp);
	errNum |= clSetKernelArg(_cl_kernel_scale_alpha_dev, 3, sizeof(float), scale);
	CHKERR(errNum, "setting kernel _cl_kernel_scale_alpha_dev arguments");
	
	errNum = clEnqueueNDRangeKernel(commands, _cl_kernel_scale_alpha_dev, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
	clFinish(commands);
    START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "_cl_kernel_scale_alpha_dev Kernel", ocdTempTimer)
    END_TIMER(ocdTempTimer)
	CHKERR(errNum, "queuing kernel _cl_kernel_scale_alpha_dev for execution");
    
    /* Initilialize log likelihood */
    log_lik = log10(scale[0]);
    
    /* Calculate the rest of the alpha variables */
    for (t = 1; t < length; t++) {
    
        /* Calculate offsets */
        offset_prev = (t - 1) * nstates;
        offset_cur = t * nstates;
        
        /* Multiply transposed A matrix by alpha(t-1) */
        /* Note: the matrix is auto-transposed by cublas reading column-major */
        // mat_vec_mul( 'N', nstates, nstates, 1.0f, a_d, nstates, 
        //              alpha_d + offset_prev, 1, 0, alpha_d + offset_cur, 1 );
        mat_vec_mul( 'n', nstates, nstates, a_d, nstates, 
                      alpha_d, offset_prev, alpha_d, offset_cur);
        
        /* Calculate alpha(t) */
        // calc_alpha_dev<<<nblocks, threads_per_block>>>( nstates, 
        //                                                 alpha_d + offset_cur, 
        //                                                 b_d, 
        //                                                 obs[t]);
		
		errNum  = clSetKernelArg(_cl_kernel_calc_alpha_dev, 0, sizeof(int), &nstates);
		errNum |= clSetKernelArg(_cl_kernel_calc_alpha_dev, 1, sizeof(cl_mem), &alpha_d);
		errNum |= clSetKernelArg(_cl_kernel_calc_alpha_dev, 2, sizeof(int), &offset_cur);
		errNum |= clSetKernelArg(_cl_kernel_calc_alpha_dev, 3, sizeof(cl_mem), &b_d);
		errNum |= clSetKernelArg(_cl_kernel_calc_alpha_dev, 4, sizeof(int), obs + t);
		CHKERR(errNum, "setting kernel _cl_kernel_calc_alpha_dev arguments");
		
		errNum = clEnqueueNDRangeKernel(commands, _cl_kernel_calc_alpha_dev, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
		clFinish(commands);
    	START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "_cl_kernel_calc_alpha_dev Kernel", ocdTempTimer)
    	END_TIMER(ocdTempTimer)
		CHKERR(errNum, "queuing kernel _cl_kernel_calc_alpha_dev for execution");
                                                        
        /* Sum alpha values to get scaling factor */
        scale[t] = dot_production(nstates, alpha_d, offset_cur, ones_n_d, 0);
        
        /* Scale alpha values */
        // scale_alpha_dev<<<nblocks, threads_per_block>>>(nstates, 
        //                                                alpha_d + offset_cur, 
        //                                                 scale[t]);
		errNum  = clSetKernelArg(_cl_kernel_scale_alpha_dev, 0, sizeof(int), &nstates);
		errNum |= clSetKernelArg(_cl_kernel_scale_alpha_dev, 1, sizeof(cl_mem), &alpha_d);
		errNum |= clSetKernelArg(_cl_kernel_scale_alpha_dev, 2, sizeof(int), &offset_cur);
		errNum |= clSetKernelArg(_cl_kernel_scale_alpha_dev, 3, sizeof(float), scale + t);
        CHKERR(errNum, "setting kernel _cl_kernel_scale_alpha_dev arguments");
		
		errNum = clEnqueueNDRangeKernel(commands, _cl_kernel_scale_alpha_dev, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
		clFinish(commands);
    	START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "_cl_kernel_scale_alpha_dev Kernel", ocdTempTimer)
    	END_TIMER(ocdTempTimer)
		CHKERR(errNum, "queuing kernel _cl_kernel_scale_alpha_dev for execution");                                                
        /* Update log likelihood */
        log_lik += log10(scale[t]);
    }
    
    return log_lik;
}

/* Calculates the backward variables (beta) */
int calc_beta()
{

    int threads_per_block;
    int nblocks;
    int t;
	int offset;
    cl_int errNum;
    
    /* Initialize beta variables */
    threads_per_block = MAX_THREADS_PER_BLOCK;
    nblocks = (nstates + threads_per_block - 1) / threads_per_block;
	size_t globalWorkSize[1] = { nblocks * threads_per_block };
	size_t localWorkSize[1]  = { threads_per_block };
	
    // init_beta_dev<<<nblocks, threads_per_block>>>(  nstates, beta_d + 
    //                                                 ((length - 1) * nstates),
    //                                                 scale[length - 1]);
	offset  = ((length - 1) * nstates);
	errNum  = clSetKernelArg(_cl_kernel_init_beta_dev, 0, sizeof(int), &nstates);
	errNum |= clSetKernelArg(_cl_kernel_init_beta_dev, 1, sizeof(cl_mem), &beta_d);
	errNum |= clSetKernelArg(_cl_kernel_init_beta_dev, 2, sizeof(int), &offset);
	errNum |= clSetKernelArg(_cl_kernel_init_beta_dev, 3, sizeof(float), scale + length - 1);
	CHKERR(errNum, "setting kernel _cl_kernel_init_beta_dev arguments");
	
	errNum = clEnqueueNDRangeKernel(commands, _cl_kernel_init_beta_dev, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
	clFinish(commands);
    START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "_cl_kernel_init_beta_dev Kernel", ocdTempTimer)
    END_TIMER(ocdTempTimer)
	CHKERR(errNum, "queuing kernel _cl_kernel_init_beta_dev for execution");
	
                                                    
    /* Calculate the rest of the beta variables */
    for (t = length - 2; t >= 0; t--) {
        
        /* Calculate first step of beta: B.*beta/scale */
        // calc_beta_dev<<<nblocks, threads_per_block>>>(  beta_d, 
        //                                                 b_d, 
        //                                                 scale[t], 
        //                                                 nstates, 
        //                                                 obs[t+1], 
        //                                                 t);
		errNum  = clSetKernelArg(_cl_kernel_calc_beta_dev, 0, sizeof(cl_mem), &beta_d);
		errNum |= clSetKernelArg(_cl_kernel_calc_beta_dev, 1, sizeof(cl_mem), &b_d);
		errNum |= clSetKernelArg(_cl_kernel_calc_beta_dev, 2, sizeof(float), scale + t);
		errNum |= clSetKernelArg(_cl_kernel_calc_beta_dev, 3, sizeof(int), &nstates);
		errNum |= clSetKernelArg(_cl_kernel_calc_beta_dev, 4, sizeof(int), obs + t + 1);
		errNum |= clSetKernelArg(_cl_kernel_calc_beta_dev, 5, sizeof(int), &t);
		CHKERR(errNum, "setting kernel _cl_kernel_calc_beta_dev arguments");
		
		errNum = clEnqueueNDRangeKernel(commands, _cl_kernel_calc_beta_dev, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
		clFinish(commands);
    	START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "_cl_kernel_calc_beta_dev Kernel", ocdTempTimer)
    	END_TIMER(ocdTempTimer)
		CHKERR(errNum, "queuing kernel _cl_kernel_calc_beta_dev for execution");
                                                        
        /* Multiply transposed A matrix by beta(t) */
        // mat_vec_mul( 'T', nstates, nstates, 1.0f, a_d, nstates, 
        //             beta_d + (t * nstates), 1, 0, beta_d + (t * nstates), 1 );
		mat_vec_mul( 'n', nstates, nstates, a_d, nstates, 
                      beta_d, t * nstates, beta_d, t * nstates);
                                                        
    }
    
    return 0;
}

/* Calculates the gamma sum */
void calc_gamma_sum()
{
    int threads_per_block;
    int nblocks;
    int size;
    int t;
	cl_int errNum;
    
    threads_per_block = MAX_THREADS_PER_BLOCK;
    nblocks = (nstates + threads_per_block - 1) / threads_per_block;
	size_t globalWorkSize[1] = { nblocks * threads_per_block };
	size_t localWorkSize[1]  = { threads_per_block };
    
	// init to zeros
	int *gamma_sum_zeros = (int*)calloc(nstates, sizeof(float));
	errNum = clEnqueueWriteBuffer(commands, gamma_sum_d, CL_TRUE, 0, sizeof(float) * nstates, gamma_sum_zeros, 0, NULL, &ocdTempEvent);
	clFinish(commands);
    START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "gamma_sum_d Data Copy", ocdTempTimer)
    END_TIMER(ocdTempTimer)
	CHKERR(errNum, "writing buffer gamma_sum_d");
	free(gamma_sum_zeros);
	
    
    /* Find sum of gamma variables */
    for (t = 0; t < length; t++) {
        // calc_gamma_dev<<<nblocks, threads_per_block>>>( gamma_sum_d,
                                                        // alpha_d,
                                                        // beta_d,
                                                        // nstates,
                                                        // t);
		errNum  = clSetKernelArg(_cl_kernel_calc_gamma_dev, 0, sizeof(cl_mem), &gamma_sum_d);
		errNum |= clSetKernelArg(_cl_kernel_calc_gamma_dev, 1, sizeof(cl_mem), &alpha_d);
		errNum |= clSetKernelArg(_cl_kernel_calc_gamma_dev, 2, sizeof(cl_mem), &beta_d);
		errNum |= clSetKernelArg(_cl_kernel_calc_gamma_dev, 3, sizeof(int), &nstates);
		errNum |= clSetKernelArg(_cl_kernel_calc_gamma_dev, 4, sizeof(int), &t);
		CHKERR(errNum, "setting kernel _cl_kernel_calc_gamma_dev arguments");
		
		errNum = clEnqueueNDRangeKernel(commands, _cl_kernel_calc_gamma_dev, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
		clFinish(commands);
    	START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "_cl_kernel_calc_gamma_dev Kernel", ocdTempTimer)
    	END_TIMER(ocdTempTimer)
		CHKERR(errNum, "queuing kernel _cl_kernel_calc_gamma_dev for execution");
    }

}

/* Calculates the sum of xi variables */
int calc_xi_sum()
{
    float sum_ab;
    int nblocks;
    int size;
    int t;
    cl_int errNum;
    
    // init to zeros
	int *xi_sum_d_zeros = (int*)calloc(nstates * nstates, sizeof(float));
	errNum = clEnqueueWriteBuffer(commands, xi_sum_d, CL_TRUE, 0, sizeof(float) * nstates * nstates, xi_sum_d_zeros, 0, NULL, &ocdTempEvent);
	clFinish(commands);
    START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "xi_sum_d Data Copy", ocdTempTimer)
    END_TIMER(ocdTempTimer)
	CHKERR(errNum, "writing buffer xi_sum_d");
	free(xi_sum_d_zeros);
    
    /* Calculate running xi sum */
    // dim3 threads(BLOCK_DIM, BLOCK_DIM);
    nblocks = (nstates + BLOCK_DIM - 1) / BLOCK_DIM;
    // dim3 grid(nblocks, nblocks);
	
	size_t globalWorkSize[2] = { nblocks * BLOCK_DIM, nblocks * BLOCK_DIM };
	size_t localWorkSize[2]  = { BLOCK_DIM, BLOCK_DIM };
    
    /* Find the sum of xi variables */
    for (t = 0; t < length - 1; t++) {
    
        /* Calculate denominator */
        sum_ab = dot_production(nstates, alpha_d, t * nstates, 
                                        beta_d, t * nstates);
        
        /* Calculate xi sum */
        // calc_xi_dev<<<grid, threads>>>( xi_sum_d,
                                        // a_d,
                                        // b_d,
                                        // alpha_d,
                                        // beta_d,
                                        // sum_ab,
                                        // nstates,
                                        // obs[t+1],
                                        // t);
		errNum  = clSetKernelArg(_cl_kernel_calc_xi_dev, 0, sizeof(cl_mem), &xi_sum_d);
		errNum |= clSetKernelArg(_cl_kernel_calc_xi_dev, 1, sizeof(cl_mem), &a_d);
		errNum |= clSetKernelArg(_cl_kernel_calc_xi_dev, 2, sizeof(cl_mem), &b_d);
		errNum |= clSetKernelArg(_cl_kernel_calc_xi_dev, 3, sizeof(cl_mem), &alpha_d);
		errNum |= clSetKernelArg(_cl_kernel_calc_xi_dev, 4, sizeof(cl_mem), &beta_d);
		errNum |= clSetKernelArg(_cl_kernel_calc_xi_dev, 5, sizeof(float), &sum_ab);
		errNum |= clSetKernelArg(_cl_kernel_calc_xi_dev, 6, sizeof(int), &nstates);
		errNum |= clSetKernelArg(_cl_kernel_calc_xi_dev, 7, sizeof(int), obs + t + 1);
		errNum |= clSetKernelArg(_cl_kernel_calc_xi_dev, 8, sizeof(int), &t);
		CHKERR(errNum, "setting kernel _cl_kernel_calc_xi_dev arguments");
		
		errNum = clEnqueueNDRangeKernel(commands, _cl_kernel_calc_xi_dev, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
		clFinish(commands);
    	START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "_cl_kernel_calc_xi_dev Kernel", ocdTempTimer)
    	END_TIMER(ocdTempTimer)
		CHKERR(errNum, "queuing kernel _cl_kernel_calc_xi_dev for execution");
   }
                                  
    return 0;
}

/* Re-estimates the state transition probabilities (A) */
int estimate_a()
{
    float sum_ab;
    int nblocks;
    cl_int errNum;
    
    /* Calculate running xi sum */
    // dim3 threads(BLOCK_DIM, BLOCK_DIM);
    nblocks = (nstates + BLOCK_DIM - 1) / BLOCK_DIM;
    // dim3 grid(nblocks, nblocks);
	
	size_t globalWorkSize[2] = { nblocks * BLOCK_DIM, nblocks * BLOCK_DIM };
	size_t localWorkSize[2]  = { BLOCK_DIM, BLOCK_DIM };
    
    /* Calculate denominator */
    sum_ab = dot_production(nstates, alpha_d, (length - 1) * nstates, 
                                    beta_d, (length - 1) * nstates);
    
    /* Calculate new value of A */
    // est_a_dev<<<grid, threads>>>(   a_d,
                                    // alpha_d,
                                    // beta_d,
                                    // xi_sum_d,
                                    // gamma_sum_d,
                                    // sum_ab,
                                    // nstates,
                                    // length);
	errNum  = clSetKernelArg(_cl_kernel_est_a_dev, 0, sizeof(cl_mem), &a_d);
	errNum |= clSetKernelArg(_cl_kernel_est_a_dev, 1, sizeof(cl_mem), &alpha_d);
	errNum |= clSetKernelArg(_cl_kernel_est_a_dev, 2, sizeof(cl_mem), &beta_d);
	errNum |= clSetKernelArg(_cl_kernel_est_a_dev, 3, sizeof(cl_mem), &xi_sum_d);
	errNum |= clSetKernelArg(_cl_kernel_est_a_dev, 4, sizeof(cl_mem), &gamma_sum_d);
	errNum |= clSetKernelArg(_cl_kernel_est_a_dev, 5, sizeof(float), &sum_ab);
	errNum |= clSetKernelArg(_cl_kernel_est_a_dev, 6, sizeof(int), &nstates);
	errNum |= clSetKernelArg(_cl_kernel_est_a_dev, 7, sizeof(int), &length);
	CHKERR(errNum, "setting kernel _cl_kernel_est_a_dev arguments");
	
	errNum = clEnqueueNDRangeKernel(commands, _cl_kernel_est_a_dev, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
	clFinish(commands);
    START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "_cl_kernel_est_a_dev Kernel", ocdTempTimer)
    END_TIMER(ocdTempTimer)
	CHKERR(errNum, "queuing kernel _cl_kernel_est_a_dev for execution");
      
    /* Sum rows of A to get scaling values */
    // mat_vec_mul( 'T', nstates, nstates, 1.0f, a_d, nstates, 
                // ones_n_d, 1, 0, c_d, 1 );
	mat_vec_mul( 't', nstates, nstates, a_d, nstates, 
                 ones_n_d, 0, c_d, 0);
    
    
    /* Normalize A matrix */
    // scale_a_dev<<<grid, threads>>>( a_d,
                                    // c_d,
                                    // nstates);
	errNum  = clSetKernelArg(_cl_kernel_scale_a_dev, 0, sizeof(cl_mem), &a_d);
	errNum |= clSetKernelArg(_cl_kernel_scale_a_dev, 1, sizeof(cl_mem), &c_d);
	errNum |= clSetKernelArg(_cl_kernel_scale_a_dev, 2, sizeof(int), &nstates);
	CHKERR(errNum, "setting kernel _cl_kernel_scale_a_dev arguments");
	
	errNum = clEnqueueNDRangeKernel(commands, _cl_kernel_scale_a_dev, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
	clFinish(commands);
    START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "_cl_kernel_scale_a_dev Kernel", ocdTempTimer)
    END_TIMER(ocdTempTimer)
	CHKERR(errNum, "queuing kernel _cl_kernel_scale_a_dev for execution");

    return 0;
}

/* Re-estimates the output symbol probabilities (B) */
int estimate_b()
{

    float sum_ab;
    int size;
    int t;
    int grid_x, grid_y;
    cl_int errNum;

	// init to zeros
	int *b_d_zeros = (int*)calloc(nstates * nsymbols, sizeof(float));
	errNum = clEnqueueWriteBuffer(commands, b_d, CL_TRUE, 0, sizeof(float) * nstates * nsymbols, b_d_zeros, 0, NULL, &ocdTempEvent);
	clFinish(commands);
    START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "b_d Data Copy", ocdTempTimer)
    END_TIMER(ocdTempTimer)
	CHKERR(errNum, "writing buffer b_d");
	free(b_d_zeros);
    
    /* Calculate number of threads and blocks needed */
    grid_x = (nstates + BLOCK_DIM - 1) / BLOCK_DIM;
    grid_y = (nsymbols + BLOCK_DIM - 1) / BLOCK_DIM;

	size_t globalWorkSize[2] = { grid_x * BLOCK_DIM, grid_y * BLOCK_DIM };
	size_t localWorkSize[2]  = { BLOCK_DIM, BLOCK_DIM };
				
    
    for (t = 0; t < length; t++) {
        
        /* Calculate denominator */
        sum_ab = dot_production(nstates, alpha_d, t * nstates, 
                                        beta_d, t * nstates);
        
        /* Accumulate B values */
        // acc_b_dev<<<grid, threads>>>(   b_d, 
                                        // alpha_d, 
                                        // beta_d, 
                                        // sum_ab,
                                        // nstates, 
                                        // nsymbols, 
                                        // obs[t], 
                                        // t);
		errNum  = clSetKernelArg(_cl_kernel_acc_b_dev, 0, sizeof(cl_mem), &b_d);
		errNum |= clSetKernelArg(_cl_kernel_acc_b_dev, 1, sizeof(cl_mem), &alpha_d);
		errNum |= clSetKernelArg(_cl_kernel_acc_b_dev, 2, sizeof(cl_mem), &beta_d);
		errNum |= clSetKernelArg(_cl_kernel_acc_b_dev, 3, sizeof(float), &sum_ab);
		errNum |= clSetKernelArg(_cl_kernel_acc_b_dev, 4, sizeof(int), &nstates);
		errNum |= clSetKernelArg(_cl_kernel_acc_b_dev, 5, sizeof(int), &nsymbols);
		errNum |= clSetKernelArg(_cl_kernel_acc_b_dev, 6, sizeof(int), obs + t);
		errNum |= clSetKernelArg(_cl_kernel_acc_b_dev, 7, sizeof(int), &t);
		CHKERR(errNum, "setting kernel _cl_kernel_acc_b_dev arguments");
		
		errNum = clEnqueueNDRangeKernel(commands, _cl_kernel_acc_b_dev, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
		clFinish(commands);
    	START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "_cl_kernel_acc_b_dev Kernel", ocdTempTimer)
    	END_TIMER(ocdTempTimer)
		CHKERR(errNum, "queuing kernel _cl_kernel_acc_b_dev for execution");
                                        
    }
    
    /* Re-estimate B values */
    // est_b_dev<<<grid, threads>>>(b_d, gamma_sum_d, nstates, nsymbols);
	errNum  = clSetKernelArg(_cl_kernel_est_b_dev, 0, sizeof(cl_mem), &b_d);
	errNum |= clSetKernelArg(_cl_kernel_est_b_dev, 1, sizeof(cl_mem), &gamma_sum_d);
	errNum |= clSetKernelArg(_cl_kernel_est_b_dev, 2, sizeof(int), &nstates);
	errNum |= clSetKernelArg(_cl_kernel_est_b_dev, 3, sizeof(int), &nsymbols);
	CHKERR(errNum, "setting kernel _cl_kernel_est_b_dev arguments");
	
	errNum = clEnqueueNDRangeKernel(commands, _cl_kernel_est_b_dev, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
	clFinish(commands);
    START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "_cl_kernel_est_b_dev Kernel", ocdTempTimer)
    END_TIMER(ocdTempTimer)
	CHKERR(errNum, "queuing kernel _cl_kernel_est_b_dev for execution");
	
    
    /* Sum rows of B to get scaling values */
    // mat_vec_mul( 'N', nstates, nsymbols, 1.0f, b_d, nstates, 
                // ones_s_d, 1, 0, c_d, 1 );
	mat_vec_mul( 'N', nstates, nsymbols, b_d, nstates, 
                 ones_s_d, 0, c_d, 0);


    /* Normalize B matrix */
    // scale_b_dev<<<grid, threads>>>( b_d,
                                    // c_d,
                                    // nstates,
                                    // nsymbols);
	errNum  = clSetKernelArg(_cl_kernel_scale_b_dev, 0, sizeof(cl_mem), &b_d);
	errNum |= clSetKernelArg(_cl_kernel_scale_b_dev, 1, sizeof(cl_mem), &c_d);
	errNum |= clSetKernelArg(_cl_kernel_scale_b_dev, 2, sizeof(int), &nstates);
	errNum |= clSetKernelArg(_cl_kernel_scale_b_dev, 3, sizeof(int), &nsymbols);
	CHKERR(errNum, "setting kernel _cl_kernel_scale_b_dev arguments");
	
	errNum = clEnqueueNDRangeKernel(commands, _cl_kernel_scale_b_dev, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
	clFinish(commands);
    START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "_cl_kernel_scale_b_dev Kernel", ocdTempTimer)
    END_TIMER(ocdTempTimer)
	CHKERR(errNum, "queuing kernel _cl_kernel_scale_b_dev for execution");
        
    return 0;
}    

/* Re-estimates the initial state probabilities (Pi) */
int estimate_pi()
{

    float sum_ab;
    int threads_per_block;
    int nblocks;
    cl_int errNum;
    
    /* Calculate denominator */
    sum_ab = dot_production(nstates, alpha_d, 0, beta_d, 0);
    
    /* Estimate Pi values */
    threads_per_block = MAX_THREADS_PER_BLOCK;
    nblocks = (nstates + threads_per_block - 1) / threads_per_block;
	size_t globalWorkSize[1] = { nblocks * threads_per_block };
	size_t localWorkSize[1]  = { threads_per_block };
	
    // est_pi_dev<<<nblocks, threads_per_block>>>( pi_d, 
                                                // alpha_d, 
                                                // beta_d, 
                                                // sum_ab, 
                                                // nstates);
	errNum  = clSetKernelArg(_cl_kernel_est_pi_dev, 0, sizeof(cl_mem), &pi_d);
	errNum |= clSetKernelArg(_cl_kernel_est_pi_dev, 1, sizeof(cl_mem), &alpha_d);
	errNum |= clSetKernelArg(_cl_kernel_est_pi_dev, 2, sizeof(cl_mem), &beta_d);
	errNum |= clSetKernelArg(_cl_kernel_est_pi_dev, 3, sizeof(float), &sum_ab);
	errNum |= clSetKernelArg(_cl_kernel_est_pi_dev, 4, sizeof(int), &nstates);
	CHKERR(errNum, "setting kernel _cl_kernel_est_pi_dev arguments");
	
	errNum = clEnqueueNDRangeKernel(commands, _cl_kernel_est_pi_dev, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
	clFinish(commands);
    START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "_cl_kernel_est_pi_dev Kernel", ocdTempTimer)
    END_TIMER(ocdTempTimer)
	CHKERR(errNum, "queuing kernel _cl_kernel_est_pi_dev for execution");

    return 0;
}


/*******************************************************************************
 * BWA function
 */

/* Runs the Baum-Welch Algorithm on the supplied HMM and observation sequence */
float run_hmm_bwa(  Hmm *hmm, 
                    Obs *in_obs, 
                    int iterations, 
                    float threshold)
{

    cl_int errNum;
	
    /* Host-side variables */
    float *a;
    float *b;
    float *pi;
    int threads_per_block;
    int nblocks;
    int size;
    float new_log_lik;
    float old_log_lik = 0;
    int iter;
	/* init the opencl context and kernels */
	
	init_cl();
    
    /* Initialize HMM values */
    a = hmm->a;
    b = hmm->b;
    pi = hmm->pi;
    nsymbols = hmm->nsymbols;
    nstates = hmm->nstates;
    obs = in_obs->data;
    length = in_obs->length;
    
    /* Allocate host memory */
    scale = (float *) malloc(sizeof(float) * length);
    if (scale == 0) {
        fprintf (stderr, "ERROR: Host memory allocation error (scale)\n");
        return EXIT_ERROR;
    }
    
    /* Allocate device memory */
	a_d = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * nstates * nstates, NULL, NULL);
	b_d = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * nstates * nsymbols, NULL, NULL);
	pi_d = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * nstates, NULL, NULL);
	alpha_d = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * nstates * length, NULL, NULL);
	beta_d = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * nstates * length, NULL, NULL);
	gamma_sum_d = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * nstates, NULL, NULL);
	xi_sum_d = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * nstates * nstates, NULL, NULL);
	c_d = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * nstates, NULL, NULL);
	ones_n_d = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * nstates, NULL, NULL);
	ones_s_d = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float) * nsymbols, NULL, NULL);
	
	CHECK_NULL_ERROR( a_d, "Error creating buffer for a_d");
	CHECK_NULL_ERROR( b_d, "Error creating buffer for  b_d");
	CHECK_NULL_ERROR( pi_d, "Error creating buffer for  pi_d");
	CHECK_NULL_ERROR( alpha_d, "Error creating buffer for  alpha_d");
	CHECK_NULL_ERROR( beta_d, "Error creating buffer for  beta_d");
	CHECK_NULL_ERROR( gamma_sum_d, "Error creating buffer for  gamma_sum_d");
	CHECK_NULL_ERROR( xi_sum_d, "Error creating buffer for  xi_sum_d");
	CHECK_NULL_ERROR( c_d, "Error creating buffer for  c_d");
	CHECK_NULL_ERROR( ones_n_d, "Error creating buffer for  ones_n_d");
	CHECK_NULL_ERROR( ones_s_d, "Error creating buffer for  ones_s_d");
	
	/* Transfer device data */
	errNum = clEnqueueWriteBuffer(commands, a_d, CL_TRUE, 0, sizeof(float) * nstates * nstates, a, 0, NULL, &ocdTempEvent);
	clFinish(commands);
    START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "a_d Data Copy", ocdTempTimer)
    END_TIMER(ocdTempTimer)
	CHKERR(errNum, "Error writing buffer a_d");
	
	errNum = clEnqueueWriteBuffer(commands, b_d, CL_TRUE, 0, sizeof(float) * nstates * nsymbols, b, 0, NULL, &ocdTempEvent);
	clFinish(commands);
    START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "b_d Data Copy", ocdTempTimer)
    END_TIMER(ocdTempTimer)
	CHKERR(errNum, "Error writing buffer b_d");

	errNum = clEnqueueWriteBuffer(commands, pi_d, CL_TRUE, 0, sizeof(float) * nstates, pi, 0, NULL, &ocdTempEvent);
	clFinish(commands);
    START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "pi_d Data Copy", ocdTempTimer)
    END_TIMER(ocdTempTimer)
	CHKERR(errNum, "Error writing buffer pi_d");
    
    /* Initialize ones array */
    threads_per_block = MAX_THREADS_PER_BLOCK;
    nblocks = (nstates + threads_per_block - 1) / threads_per_block;
	size_t globalWorkSize[1] = { nblocks * threads_per_block };
	size_t localWorkSize[1]  = { threads_per_block };
	
	// init_ones_dev<<<nblocks, threads_per_block>>>(ones_s_d, nsymbols);
    errNum  = clSetKernelArg(_cl_kernel_init_ones_dev, 0, sizeof(cl_mem), &ones_s_d);
	errNum |= clSetKernelArg(_cl_kernel_init_ones_dev, 1, sizeof(int), &nsymbols);
	CHKERR(errNum, "Error setting kernel _cl_kernel_init_ones_dev arguments");
	
	errNum = clEnqueueNDRangeKernel(commands, _cl_kernel_init_ones_dev, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
	clFinish(commands);
    START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "_cl_kernel_init_ones_dev Kernel", ocdTempTimer)
    END_TIMER(ocdTempTimer)
	CHKERR(errNum, "Error queuing kernel _cl_kernel_init_ones_dev for execution");
	
    /* Run BWA for either max iterations or until threshold is reached */
    for (iter = 0; iter < iterations; iter++) {
    
        new_log_lik = calc_alpha();
        if (new_log_lik == EXIT_ERROR) {
            return EXIT_ERROR;
        }
        
        if (calc_beta() == EXIT_ERROR) {
            return EXIT_ERROR;
        }
        
        calc_gamma_sum();
        
        if (calc_xi_sum() == EXIT_ERROR) {
            return EXIT_ERROR;
        }
        
        if (estimate_a() == EXIT_ERROR) {
            return EXIT_ERROR;
        }
        
        if (estimate_b() == EXIT_ERROR) {
            return EXIT_ERROR;
        }
        
        if (estimate_pi() == EXIT_ERROR) {
            return EXIT_ERROR;
        }

        /* check log_lik vs. threshold */
        if (threshold > 0 && iter > 0) {
            if (fabs(pow(10,new_log_lik) - pow(10,old_log_lik)) < threshold) {
                break;
            }
        }
            
        old_log_lik = new_log_lik;   

    }
    
    /* Copy device variables back to host */
    errNum = clEnqueueReadBuffer(commands, a_d, CL_TRUE, 0, sizeof(float) * nstates * nstates, a, 0, NULL, &ocdTempEvent);
    clFinish(commands);
    START_TIMER(ocdTempEvent, OCD_TIMER_D2H, "a_d Data Copy", ocdTempTimer)
    END_TIMER(ocdTempTimer)
	CHKERR(errNum, "Error reading buffer a_d");
	
	errNum = clEnqueueReadBuffer(commands, b_d, CL_TRUE, 0, sizeof(float) * nstates * nsymbols, b, 0, NULL, &ocdTempEvent);
	clFinish(commands);
    START_TIMER(ocdTempEvent, OCD_TIMER_D2H, "b_d Data Copy", ocdTempTimer)
    END_TIMER(ocdTempTimer)
	CHKERR(errNum, "Error reading buffer b_d");
	
	errNum = clEnqueueReadBuffer(commands, pi_d, CL_TRUE, 0, sizeof(float) * nstates, pi, 0, NULL, &ocdTempEvent);
	clFinish(commands);
    START_TIMER(ocdTempEvent, OCD_TIMER_D2H, "pi_d Data Copy", ocdTempTimer)
    END_TIMER(ocdTempTimer)
	CHKERR(errNum, "Error reading buffer pi_d");
	
    clFinish(commands);
    /* Free memory */
    free(scale);

	clReleaseMemObject(a_d);
	clReleaseMemObject(b_d);
	clReleaseMemObject(pi_d);
	clReleaseMemObject(alpha_d);
	clReleaseMemObject(beta_d);
	clReleaseMemObject(gamma_sum_d);
	clReleaseMemObject(xi_sum_d);
	clReleaseMemObject(c_d);
	clReleaseMemObject(ones_n_d);
	clReleaseMemObject(ones_s_d);
	
	clReleaseKernel(_cl_kernel_init_ones_dev);
	clReleaseKernel(_cl_kernel_init_alpha_dev);
	clReleaseKernel(_cl_kernel_calc_alpha_dev);
	clReleaseKernel(_cl_kernel_scale_alpha_dev);
	clReleaseKernel(_cl_kernel_init_beta_dev);
	clReleaseKernel(_cl_kernel_calc_beta_dev);
	clReleaseKernel(_cl_kernel_calc_gamma_dev);
	clReleaseKernel(_cl_kernel_calc_xi_dev);
	clReleaseKernel(_cl_kernel_est_a_dev);
	clReleaseKernel(_cl_kernel_scale_a_dev);
	clReleaseKernel(_cl_kernel_acc_b_dev);
	clReleaseKernel(_cl_kernel_est_b_dev);
	clReleaseKernel(_cl_kernel_scale_b_dev);
	clReleaseKernel(_cl_kernel_est_pi_dev);
	clReleaseKernel(_cl_kernel_s_dot_kernel_naive);
	clReleaseKernel(_cl_kernel_sgemvt_kernel_naive);
	clReleaseKernel(_cl_kernel_sgemvn_kernel_naive);

	clReleaseProgram(_cl_program);

	clReleaseCommandQueue(commands);
	clReleaseContext(context);
  
    return new_log_lik;
}

static struct option size_opts[] = 
{
	/* name, has_tag, flag, val*/
	{"state number", 1, NULL, 'n'},
	{"symbol number", 1, NULL, 's'},
	{"observation number", 1, NULL, 't'},
	{"varying mode", 1, NULL, 'v'},
	{0, 0, 0, 0}	
};

/* Time the forward algorithm and vary the number of states */ 
int main(int argc, char *argv[]) 
{
    /* Initialize variables */
    Hmm *hmm;                /* Initial HMM */
    Obs *obs;                /* Observation sequence */
    float *a;
    float *b;
    float *pi;
    int *obs_seq;
    float log_lik;           /* Output likelihood of FO */
    int mul;
    int m;
	int s = S, t = T;
    int n = N;
    int i;
    struct timeval timer;
    
    printf("Starting bwa_hmm\n");
    ocd_init(&argc, &argv, NULL);
	
	int opt, opt_index = 0;
	char v_model;
	while((opt = getopt_long(argc, argv, "::n:s:t:v:", size_opts, &opt_index)) != -1)
	{
		//printf("here");
		switch(opt)
		{
			case 'v':
				v_model = optarg[0];
				break;
			case 'n':
				n = atoi(optarg);
				break;
			case 's':
				s = atoi(optarg);
				break;
			case 't':
				t = atoi(optarg);
				break;
			default:
				fprintf(stderr, "Usage %s [-n state number | -s symbol number | -t observation number] [-v varying model]\n", argv[0]);
				exit(EXIT_FAILURE);
		}	
	} 
  
    /* Initialize HMM and observation sequence */
    hmm = (Hmm *)malloc(sizeof(Hmm));
    obs = (Obs *)malloc(sizeof(Obs));
	
	if(v_model == 'n')
	{
		/* Create observation sequence */
		obs->length = T;
		obs_seq = (int *)malloc(sizeof(int) * T);
		for (i = 0; i < T; i++) {
			obs_seq[i] = 0;
		}
		obs->data = obs_seq;
		
		// printf("Vary states with S = %i, T = %i\n", S, T);
		// printf("States\tTime (s)\tLog_likelihood\n");
		
		/* Run timed tests from 1*mul to 9*mul states */
		if (n >= 8000) {
			return 0;
		}
		// n = 7000;           
		/* Assign HMM parameters */
		hmm->nstates = n;
		hmm->nsymbols = S;
		a = (float *)malloc(sizeof(float) * n * n);
		for (i = 0; i < (n * n); i++) {
			a[i] = 1.0f/(float)n;
		}
		hmm->a = a;
		b = (float *)malloc(sizeof(float) * n * S);
		for (i = 0; i < (n * S); i++) {
			b[i] = 1.0f/(float)S;
		}
		hmm->b = b;
		pi = (float *)malloc(sizeof(float) * n);
		for (i = 0; i < n; i++) {
			pi[i] = 1.0f/(float)n;
		}
		hmm->pi = pi;
		
		/* Run the BWA on the observation sequence */
		
		tic(&timer);
		log_lik = run_hmm_bwa(hmm, obs, ITERATIONS, 0);
		printf("Observations\tTime (s)\tLog_likelihood\n");
		printf("%i\t", n);
		toc(&timer);
		printf("%f\n", log_lik);
		
		/* Free memory */
		free(a);
		free(b);
		free(pi);

		hmm->a = NULL;
		hmm->b = NULL;
		hmm->pi = NULL;
		free_vars(hmm, obs);

		// printf("\n");
	} else if(v_model == 's')
	{	
		/* Create observation sequence */
		obs->length = T;
		obs_seq = (int *)malloc(sizeof(int) * T);
		for (i = 0; i < T; i++) {
			obs_seq[i] = 0;
		}
		obs->data = obs_seq;
		
		// printf("Vary symbols with N = %i, T = %i\n", N, T);
		// printf("Symbols\tTime (s)\tLog_likelihood\n");
			
		if (s >= 8000) {
			return 0;
		}
				
		/* Assign HMM parameters */
		hmm->nstates = N;
		hmm->nsymbols = s;
		a = (float *)malloc(sizeof(float) * N * N);
		for (i = 0; i < (N * N); i++) {
			a[i] = 1.0f/(float)N;
		}
		hmm->a = a;
		b = (float *)malloc(sizeof(float) * N * s);
		for (i = 0; i < (N * s); i++) {
			b[i] = 1.0f/(float)s;
		}
		hmm->b = b;
		pi = (float *)malloc(sizeof(float) * N);
		for (i = 0; i < N; i++) {
			pi[i] = 1.0f/(float)N;
		}
		hmm->pi = pi;
		
		/* Run the BWA on the observation sequence */
		
		tic(&timer);
		log_lik = run_hmm_bwa(hmm, obs, ITERATIONS, 0);
		printf("Observations\tTime (s)\tLog_likelihood\n");
		printf("%i\t", s);
		toc(&timer);
		printf("%f\n", log_lik);
		
		/* Free memory */
		free(a);
		free(b);
		free(pi);

		hmm->a = NULL;
		hmm->b = NULL;
		hmm->pi = NULL;
		free_vars(hmm, obs);
		
		// printf("\n");
    
	} else if(v_model == 't')
	{

		if (t >= 10000) {
			return 0;
		}
		/* Create HMM */
		hmm->nstates = N;
		hmm->nsymbols = S;
		a = (float *)malloc(sizeof(float) * N * N);
		for (i = 0; i < (N * N); i++) {
			a[i] = 1.0f/(float)N;
		}
		hmm->a = a;
		b = (float *)malloc(sizeof(float) * N * S);
		for (i = 0; i < (N * S); i++) {
			b[i] = 1.0f/(float)S;
		}
		hmm->b = b;
		pi = (float *)malloc(sizeof(float) * N);
		for (i = 0; i < N; i++) {
			pi[i] = 1.0f/(float)N;
		}
		hmm->pi = pi;
		
		// printf("Vary observations with N = %i, S = %i\n", N, S);
		// printf("Observations\tTime (s)\tLog_likelihood\n");
		
		/* Create observation sequence */
		obs->length = t;
		obs_seq = (int *)malloc(sizeof(int) * t);
		for (i = 0; i < t; i++) {
			obs_seq[i] = 0;
		}
		obs->data = obs_seq;
		
		/* Run the BWA on the observation sequence */
		
		tic(&timer);
		log_lik = run_hmm_bwa(hmm, obs, ITERATIONS, 0);
		printf("Observations\tTime (s)\tLog_likelihood\n");
		printf("%i\t", t);
		toc(&timer);
		printf("%f\n", log_lik);
		
		/* Free memory */
		free(obs_seq);
		
		obs->data = NULL;
		free_vars(hmm, obs);
		
		// printf("\n");
	}
	
	ocd_finalize();//KK
    return;
}
