#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <unistd.h>
#include <math.h>
#include <sys/time.h>
#include <getopt.h>
#include <malloc.h>
#include "dwarf_cuda.h"
#include "dwarf_timer.h"
static int g_rep = 0, hmm_valid = 1; /* common timer repetition index / validation result */

#define AOCL_ALIGNMENT                    0

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

#define MAX_THREADS_PER_BLOCK DWARF_LS  /* [FIX-OCL-34] pinned work-group size 256 (was literal 256) for the 1-D kernels */
#define BLOCK_DIM 16
#define EXIT_ERROR 1

#define SDOT_BLOCK_SIZE (128)
#define SDOT_BLOCK_NUM (80)

#define MVMUL_BLOCK_SIZE (128)
#define MVMUL_BLOCK_NUM  (64)

#include "hmm_kernel.cuh"

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
__attribute__((aligned(64)))
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
float *a_d;                         /* A matrix on GPU */
float *b_d;                         /* B matrix on GPU */
float *pi_d;                        /* Pi matrix on GPU */
float *alpha_d;                     /* Forward variables (alpha) on GPU */
float *beta_d;                      /* Backward variables (beta) on GPU */
float *gamma_sum_d;                 /* Sum of gamma variables on GPU */
float *xi_sum_d;                    /* Sum of xi variables on GPU */
float *c_d;                         /* Temporary array on GPU */
float *ones_n_d;                    /* Length <states> array of 1s on GPU */
float *ones_s_d;                    /* Length <symbols> array of 1s on GPU */

/************************************************************************/
/* This function is used to compute the dot product of two single-      */
/* precision vectors. Add offset support to cl_mem type files           */
/* OUTPUTS:                                                             */
/*   returns single-precision dot product                               */
/*                                                                      */
/************************************************************************/
float dot_production(int n, float *paramA, int offsetA, float *paramB, int offsetB)
{
	float result = 0.0f;
	int i;

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
	float *partialSum = (float*) memalign ( AOCL_ALIGNMENT,n * sizeof(float));

	float *partialSum_d;
	CUDA_CHECK(cudaMalloc((void**)&partialSum_d, sizeof(float) * n));

	/* set partialSum_d to all zeros */
	CUDA_CHECK(cudaMemcpy(partialSum_d, partialSum, sizeof(float) * n, cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaDeviceSynchronize());

	s_dot_kernel_naive<<<DWARF_NBLOCKS(n), DWARF_LS>>>(n, paramA, offsetA, paramB, offsetB, partialSum_d);
	CUDA_SYNC_LAUNCH();

	CUDA_CHECK(cudaMemcpy(partialSum, partialSum_d, sizeof(float) * n, cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaDeviceSynchronize());

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
void mat_vec_mul(char trans, int m, int n, float *A, int lda, float *x, int offsetX, float *y, int offsetY)
{

	if((trans != 'n') && (trans != 't'))
	{
		return;
	}

	if(trans == 't')
	{
		mvm_trans_kernel_naive<<<DWARF_NBLOCKS((lda == m) ? m : n), DWARF_LS>>>(m, n, A, lda, x, offsetX, y, offsetY);
		CUDA_SYNC_LAUNCH();
	} else
	{
		mvm_non_kernel_naive<<<DWARF_NBLOCKS((lda == m) ? m : n), DWARF_LS>>>(m, n, A, lda, x, offsetX, y, offsetY);
		CUDA_SYNC_LAUNCH();
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

	/* Initialize alpha variables */
	threads_per_block = MAX_THREADS_PER_BLOCK;
	nblocks = (nstates + threads_per_block - 1) / threads_per_block;

	// init_alpha_dev<<<nblocks, threads_per_block>>>( b_d, 
	//                                                 pi_d, 
	//                                                 nstates, 
	//                                                 alpha_d,
	//                                                 ones_n_d,
	//                                                 obs[0]);
	init_alpha_dev<<<nblocks, threads_per_block>>>(b_d, pi_d, nstates, alpha_d, ones_n_d, obs[0]);
	CUDA_SYNC_LAUNCH();

	/* Sum alpha values to get scaling factor */
	scale[0] = dot_production(nstates, alpha_d, 0, ones_n_d, 0);

	/* Scale alpha values */
	// scale_alpha_dev<<<nblocks, threads_per_block>>>(    nstates, 
	//                                                     alpha_d, 
	//                                                     scale[0]);
	int tmp = 0;
	scale_alpha_dev<<<nblocks, threads_per_block>>>(nstates, alpha_d, tmp, scale[0]);
	CUDA_SYNC_LAUNCH();

	/* Initilialize log likelihood */
	log_lik = log10((double)scale[0]);

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

		calc_alpha_dev<<<nblocks, threads_per_block>>>(nstates, alpha_d, offset_cur, b_d, obs[t]);
		CUDA_SYNC_LAUNCH();

		/* Sum alpha values to get scaling factor */
		scale[t] = dot_production(nstates, alpha_d, offset_cur, ones_n_d, 0);

		/* Scale alpha values */
		// scale_alpha_dev<<<nblocks, threads_per_block>>>(nstates, 
		//                                                alpha_d + offset_cur, 
		//                                                 scale[t]);
		scale_alpha_dev<<<nblocks, threads_per_block>>>(nstates, alpha_d, offset_cur, scale[t]);
		CUDA_SYNC_LAUNCH();
		/* Update log likelihood */
		log_lik += log10((double)scale[t]);
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

	/* Initialize beta variables */
	threads_per_block = MAX_THREADS_PER_BLOCK;
	nblocks = (nstates + threads_per_block - 1) / threads_per_block;

	// init_beta_dev<<<nblocks, threads_per_block>>>(  nstates, beta_d + 
	//                                                 ((length - 1) * nstates),
	//                                                 scale[length - 1]);
	offset  = ((length - 1) * nstates);
	init_beta_dev<<<nblocks, threads_per_block>>>(nstates, beta_d, offset, scale[length - 1]);
	CUDA_SYNC_LAUNCH();


	/* Calculate the rest of the beta variables */
	for (t = length - 2; t >= 0; t--) {

		/* Calculate first step of beta: B.*beta/scale */
		// calc_beta_dev<<<nblocks, threads_per_block>>>(  beta_d, 
		//                                                 b_d, 
		//                                                 scale[t], 
		//                                                 nstates, 
		//                                                 obs[t+1], 
		//                                                 t);
		calc_beta_dev<<<nblocks, threads_per_block>>>(beta_d, b_d, scale[t], nstates, obs[t+1], t);
		CUDA_SYNC_LAUNCH();

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

	threads_per_block = MAX_THREADS_PER_BLOCK;
	nblocks = (nstates + threads_per_block - 1) / threads_per_block;

	// init to zeros
	int *gamma_sum_zeros = (int*)calloc(nstates, sizeof(float));
	CUDA_CHECK(cudaMemcpy(gamma_sum_d, gamma_sum_zeros, sizeof(float) * nstates, cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaDeviceSynchronize());
	free(gamma_sum_zeros);


	/* Find sum of gamma variables */
	for (t = 0; t < length; t++) {
		// calc_gamma_dev<<<nblocks, threads_per_block>>>( gamma_sum_d,
		// alpha_d,
		// beta_d,
		// nstates,
		// t);
		calc_gamma_dev<<<nblocks, threads_per_block>>>(gamma_sum_d, alpha_d, beta_d, nstates, t);
		CUDA_SYNC_LAUNCH();
	}

}

/* Calculates the sum of xi variables */
int calc_xi_sum()
{
	float sum_ab;
	int nblocks;
	int size;
	int t;

	// init to zeros
	int *xi_sum_d_zeros = (int*)calloc(nstates * nstates, sizeof(float));
	CUDA_CHECK(cudaMemcpy(xi_sum_d, xi_sum_d_zeros, sizeof(float) * nstates * nstates, cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaDeviceSynchronize());
	free(xi_sum_d_zeros);

	/* Calculate running xi sum */
	// dim3 threads(BLOCK_DIM, BLOCK_DIM);
	nblocks = (nstates + BLOCK_DIM - 1) / BLOCK_DIM;
	// dim3 grid(nblocks, nblocks);


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
		calc_xi_dev<<<dim3(nblocks, nblocks), dim3(BLOCK_DIM, BLOCK_DIM)>>>(xi_sum_d, a_d, b_d, alpha_d, beta_d, sum_ab, nstates, obs[t+1], t);
		CUDA_SYNC_LAUNCH();
	}

	return 0;
}

/* Re-estimates the state transition probabilities (A) */
int estimate_a()
{
	float sum_ab;
	int nblocks;

	/* Calculate running xi sum */
	// dim3 threads(BLOCK_DIM, BLOCK_DIM);
	nblocks = (nstates + BLOCK_DIM - 1) / BLOCK_DIM;
	// dim3 grid(nblocks, nblocks);


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
	est_a_dev<<<dim3(nblocks, nblocks), dim3(BLOCK_DIM, BLOCK_DIM)>>>(a_d, alpha_d, beta_d, xi_sum_d, gamma_sum_d, sum_ab, nstates, length);
	CUDA_SYNC_LAUNCH();

	/* Sum rows of A to get scaling values */
	// mat_vec_mul( 'T', nstates, nstates, 1.0f, a_d, nstates, 
	// ones_n_d, 1, 0, c_d, 1 );
	mat_vec_mul( 't', nstates, nstates, a_d, nstates, 
			ones_n_d, 0, c_d, 0);


	/* Normalize A matrix */
	// scale_a_dev<<<grid, threads>>>( a_d,
	// c_d,
	// nstates);
	scale_a_dev<<<dim3(nblocks, nblocks), dim3(BLOCK_DIM, BLOCK_DIM)>>>(a_d, c_d, nstates);
	CUDA_SYNC_LAUNCH();

	return 0;
}

/* Re-estimates the output symbol probabilities (B) */
int estimate_b()
{

	float sum_ab;
	int size;
	int t;
	int grid_x, grid_y;

	// init to zeros
	int *b_d_zeros = (int*)calloc(nstates * nsymbols, sizeof(float));
	CUDA_CHECK(cudaMemcpy(b_d, b_d_zeros, sizeof(float) * nstates * nsymbols, cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaDeviceSynchronize());
	free(b_d_zeros);

	/* Calculate number of threads and blocks needed */
	grid_x = (nstates + BLOCK_DIM - 1) / BLOCK_DIM;
	grid_y = (nsymbols + BLOCK_DIM - 1) / BLOCK_DIM;



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
		acc_b_dev<<<dim3(grid_x, grid_y), dim3(BLOCK_DIM, BLOCK_DIM)>>>(b_d, alpha_d, beta_d, sum_ab, nstates, nsymbols, obs[t], t);
		CUDA_SYNC_LAUNCH();

	}

	/* Re-estimate B values */
	// est_b_dev<<<grid, threads>>>(b_d, gamma_sum_d, nstates, nsymbols);
	est_b_dev<<<dim3(grid_x, grid_y), dim3(BLOCK_DIM, BLOCK_DIM)>>>(b_d, gamma_sum_d, nstates, nsymbols);
	CUDA_SYNC_LAUNCH();


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
	scale_b_dev<<<dim3(grid_x, grid_y), dim3(BLOCK_DIM, BLOCK_DIM)>>>(b_d, c_d, nstates, nsymbols);
	CUDA_SYNC_LAUNCH();

	return 0;
}    

/* Re-estimates the initial state probabilities (Pi) */
int estimate_pi()
{

	float sum_ab;
	int threads_per_block;
	int nblocks;

	/* Calculate denominator */
	sum_ab = dot_production(nstates, alpha_d, 0, beta_d, 0);

	/* Estimate Pi values */
	threads_per_block = MAX_THREADS_PER_BLOCK;
	nblocks = (nstates + threads_per_block - 1) / threads_per_block;

	// est_pi_dev<<<nblocks, threads_per_block>>>( pi_d, 
	// alpha_d, 
	// beta_d, 
	// sum_ab, 
	// nstates);
	est_pi_dev<<<nblocks, threads_per_block>>>(pi_d, alpha_d, beta_d, sum_ab, nstates);
	CUDA_SYNC_LAUNCH();

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

	dwarf_region_begin();   /* common timer: allocation + uploads + all iterations + downloads (build/JIT is in init_cl) */

	/* Initialize HMM values */
	a = hmm->a;
	b = hmm->b;
	pi = hmm->pi;
	nsymbols = hmm->nsymbols;
	nstates = hmm->nstates;
	obs = in_obs->data;
	length = in_obs->length;

	/* Allocate host memory */
	//scale = (float *) malloc(sizeof(float) * length);
        scale = (float *) memalign(AOCL_ALIGNMENT, sizeof(float) * length);
	if (scale == 0) {
		fprintf (stderr, "ERROR: Host memory allocation error (scale)\n");
		return EXIT_ERROR;
	}

	/* Allocate device memory */
	/* [FIX-OCL-16] alpha_d/beta_d get one extra row */
	CUDA_CHECK(cudaMalloc((void**)&a_d, sizeof(float) * nstates * nstates));
	CUDA_CHECK(cudaMalloc((void**)&b_d, sizeof(float) * nstates * nsymbols));
	CUDA_CHECK(cudaMalloc((void**)&pi_d, sizeof(float) * nstates));
	CUDA_CHECK(cudaMalloc((void**)&alpha_d, sizeof(float) * nstates * (length + 1)));
	CUDA_CHECK(cudaMalloc((void**)&beta_d, sizeof(float) * nstates * (length + 1)));
	CUDA_CHECK(cudaMalloc((void**)&gamma_sum_d, sizeof(float) * nstates));
	CUDA_CHECK(cudaMalloc((void**)&xi_sum_d, sizeof(float) * nstates * nstates));
	CUDA_CHECK(cudaMalloc((void**)&c_d, sizeof(float) * nstates));
	CUDA_CHECK(cudaMalloc((void**)&ones_n_d, sizeof(float) * nstates));
	CUDA_CHECK(cudaMalloc((void**)&ones_s_d, sizeof(float) * nsymbols));

	/* Transfer device data */
	CUDA_CHECK(cudaMemcpy(a_d, a, sizeof(float) * nstates * nstates, cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaDeviceSynchronize());

	CUDA_CHECK(cudaMemcpy(b_d, b, sizeof(float) * nstates * nsymbols, cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaDeviceSynchronize());

	CUDA_CHECK(cudaMemcpy(pi_d, pi, sizeof(float) * nstates, cudaMemcpyHostToDevice));
	CUDA_CHECK(cudaDeviceSynchronize());

	/* Initialize ones array */
	threads_per_block = MAX_THREADS_PER_BLOCK;
	nblocks = (nsymbols + threads_per_block - 1) / threads_per_block;  /* [FIX-OCL-34] kernel guards idx < nsymbols (was sized by nstates: left ones_s_d tail unset if nsymbols > nstates; OpenMP loops over nsymbols) */

	// init_ones_dev<<<nblocks, threads_per_block>>>(ones_s_d, nsymbols);
	init_ones_dev<<<nblocks, threads_per_block>>>(ones_s_d, nsymbols);
	CUDA_SYNC_LAUNCH();

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
			if (fabs(pow(10.0,(double)new_log_lik) - pow(10.0,(double)old_log_lik)) < threshold) {
				break;
			}
		}

		old_log_lik = new_log_lik;   

	}

	/* Copy device variables back to host */
	CUDA_CHECK(cudaMemcpy(a, a_d, sizeof(float) * nstates * nstates, cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaDeviceSynchronize());

	CUDA_CHECK(cudaMemcpy(b, b_d, sizeof(float) * nstates * nsymbols, cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaDeviceSynchronize());

	CUDA_CHECK(cudaMemcpy(pi, pi_d, sizeof(float) * nstates, cudaMemcpyDeviceToHost));
	CUDA_CHECK(cudaDeviceSynchronize());

	CUDA_CHECK(cudaDeviceSynchronize());
	dwarf_region_end(g_rep);
	/* Free memory */
	free(scale);

	cudaFree(a_d);
	cudaFree(b_d);
	cudaFree(pi_d);
	cudaFree(alpha_d);
	cudaFree(beta_d);
	cudaFree(gamma_sum_d);
	cudaFree(xi_sum_d);
	cudaFree(c_d);
	cudaFree(ones_n_d);
	cudaFree(ones_s_d);



	return new_log_lik;
}

/* Repetitions for the common timer. run_hmm_bwa() overwrites a, b and pi with the estimated model, so each repetition
 * restarts from a pristine copy. Returns the log-likelihood of the last repetition. */
static float run_hmm_bwa_reps(Hmm *hmm, Obs *in_obs, int iterations, float threshold)
{
	int ns = hmm->nstates, nsym = hmm->nsymbols, rep;
	float *a0 = (float*) malloc(sizeof(float)*ns*ns), *b0 = (float*) malloc(sizeof(float)*ns*nsym), *pi0 = (float*) malloc(sizeof(float)*ns);
	float ll = 0;
	memcpy(a0, hmm->a, sizeof(float)*ns*ns); memcpy(b0, hmm->b, sizeof(float)*ns*nsym); memcpy(pi0, hmm->pi, sizeof(float)*ns);
	dwarf_init("hmm", "graphical-models", "cuda");
	dwarf_set_size(ns, in_obs->length);
	for(rep = dwarf_first_rep(); rep < dwarf_reps(); rep++)
	{
		memcpy(hmm->a, a0, sizeof(float)*ns*ns); memcpy(hmm->b, b0, sizeof(float)*ns*nsym); memcpy(hmm->pi, pi0, sizeof(float)*ns);
		g_rep = rep;
		ll = run_hmm_bwa(hmm, in_obs, iterations, threshold);
	}
	free(a0); free(b0); free(pi0);
	/* [FIX-OCL-17] validity check (none before): the uniform HMM with an all-zero observation sequence has the closed-form
	 * log10-likelihood  length * log10(1/nsymbols). */
	{
		float expect = (float)in_obs->length * log10f(1.0f / (float)nsym);
		hmm_valid = fabsf(ll - expect) <= 1e-3f * fabsf(expect) + 1e-3f;
		printf("%s: log-likelihood %f vs analytic %f\n", hmm_valid ? "SUCCESS" : "FAILURE", ll, expect);
	}
	return ll;
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
	dwarf_cuda_init(&argc, &argv);

	int opt, opt_index = 0;
	char v_model = 0; /* [FIX-OCL-18] was uninitialised: running without -v read an indeterminate value */
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
	//hmm = (Hmm *)malloc(sizeof(Hmm));
	//obs = (Obs *)malloc(sizeof(Obs));
	hmm = (Hmm *)memalign(AOCL_ALIGNMENT,sizeof(Hmm));
	obs = (Obs *)memalign(AOCL_ALIGNMENT,sizeof(Obs));
	if(v_model == 'n')
	{
		/* Create observation sequence */
		obs->length = T;
		//obs_seq = (int *)malloc(sizeof(int) * T);
                obs_seq = (int *)memalign(AOCL_ALIGNMENT,sizeof(int) * T);
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
		//a = (float *)malloc(sizeof(float) * n * n);
                a = (float *) memalign(AOCL_ALIGNMENT,sizeof(float) * n * n);
		for (i = 0; i < (n * n); i++) {
			a[i] = 1.0f/(float)n;
		}
		hmm->a = a;
		//b = (float *)malloc(sizeof(float) * n * S);
		b = (float *) memalign(AOCL_ALIGNMENT,sizeof(float) * n * S);

		for (i = 0; i < (n * S); i++) {
			b[i] = 1.0f/(float)S;
		}
		hmm->b = b;
		//pi = (float *)malloc(sizeof(float) * n);
		pi = (float *) memalign(AOCL_ALIGNMENT,sizeof(float) * n);
		for (i = 0; i < n; i++) {
			pi[i] = 1.0f/(float)n;
		}
		hmm->pi = pi;

		/* Run the BWA on the observation sequence */

		tic(&timer);
		log_lik = run_hmm_bwa_reps(hmm, obs, ITERATIONS, 0);
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
		//obs_seq = (int *)malloc(sizeof(int) * T);
		obs_seq = (int *)  memalign(AOCL_ALIGNMENT,sizeof(int) * T);
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
		//a = (float *)malloc(sizeof(float) * N * N);
		a = (float *) memalign(AOCL_ALIGNMENT,sizeof(float) * N * N);
		for (i = 0; i < (N * N); i++) {
			a[i] = 1.0f/(float)N;
		}
		hmm->a = a;
		//b = (float *)malloc(sizeof(float) * N * s);
		b = (float *) memalign(AOCL_ALIGNMENT,sizeof(float) * N * s);
		for (i = 0; i < (N * s); i++) {
			b[i] = 1.0f/(float)s;
		}
		hmm->b = b;
		//pi = (float *)malloc(sizeof(float) * N);
		pi = (float *) memalign(AOCL_ALIGNMENT,sizeof(float) * N);
		for (i = 0; i < N; i++) {
			pi[i] = 1.0f/(float)N;
		}
		hmm->pi = pi;

		/* Run the BWA on the observation sequence */

		tic(&timer);
		log_lik = run_hmm_bwa_reps(hmm, obs, ITERATIONS, 0);
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
		//a = (float *)malloc(sizeof(float) * N * N);
		a = (float *) memalign(AOCL_ALIGNMENT,sizeof(float) * N * N);
		for (i = 0; i < (N * N); i++) {
			a[i] = 1.0f/(float)N;
		}
		hmm->a = a;
		//b = (float *)malloc(sizeof(float) * N * S);
		b = (float *) memalign(AOCL_ALIGNMENT,sizeof(float) * N * S);
		for (i = 0; i < (N * S); i++) {
			b[i] = 1.0f/(float)S;
		}
		hmm->b = b;
		//pi = (float *)malloc(sizeof(float) * N);
		pi = (float *) memalign(AOCL_ALIGNMENT,sizeof(float) * N);
		for (i = 0; i < N; i++) {
			pi[i] = 1.0f/(float)N;
		}
		hmm->pi = pi;

		// printf("Vary observations with N = %i, S = %i\n", N, S);
		// printf("Observations\tTime (s)\tLog_likelihood\n");

		/* Create observation sequence */
		obs->length = t;
		//obs_seq = (int *)malloc(sizeof(int) * t);
		obs_seq = (int *) memalign(AOCL_ALIGNMENT,sizeof(int) * t);
		for (i = 0; i < t; i++) {
			obs_seq[i] = 0;
		}
		obs->data = obs_seq;

		/* Run the BWA on the observation sequence */

		tic(&timer);
		log_lik = run_hmm_bwa_reps(hmm, obs, ITERATIONS, 0);
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

	if(v_model != 'n' && v_model != 's' && v_model != 't') { fprintf(stderr, "-v n|s|t is required\n"); return EXIT_FAILURE; }
	/* [FIX-OCL-18] `return;` in int main() returned an indeterminate status */
	return dwarf_report(hmm_valid);
}
