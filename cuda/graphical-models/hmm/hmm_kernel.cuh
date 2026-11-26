/* CUDA translation of bwa_hmm_opencl.cl (mechanical: 1-D kernels use DWARF_LS threads/block, 2-D kernels 16x16 blocks). */
#define SDOT_BLOCK_SIZE        (128)
#define SDOT_BLOCK_NUM         (80)

#define MVMUL_BLOCK_SIZE       (128)
#define MVMUL_BLOCK_NUM        (64)
#define TILE_SIZE              (32)

/* Initialize ones vector */
__global__ void init_ones_dev(  float *ones_s_d,  /* lth = nsymbols */    
		int nsymbols)              /* number of symbols */   
{
	//unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < nsymbols)
	{
		ones_s_d[idx] = 1.0f;
	}
}

/* Initialize alpha variables */
__global__ void init_alpha_dev( float *b_d,      /* dim = nsymbols x nstates */   
		float *pi_d,    /* lth = nstates */         
		int nstates,              /* number of states */      
		float *alpha_d,  /* dim = length x nstates */     
		float *ones_n_d, /* lth = nstates */       
		int obs_t)                /* current observation */     
{
	// unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < nstates)
	{
		alpha_d[idx] = pi_d[idx] * b_d[(obs_t * nstates) + idx];
		ones_n_d[idx] = 1.0f;
	}
}

/* Calculate alpha variables */
__global__ void calc_alpha_dev( int nstates,             /* number of states */     
		float *alpha_d, /* dim = length x nstates */        
		int offset,              /* offset for alpha_d */  
		float *b_d,     /* dim = nsymbols x nstates */  
		int obs_t)               /* current observation */     
{
	// unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < nstates)
	{
		alpha_d[offset + idx] = alpha_d[offset + idx] * b_d[(obs_t * nstates) + idx];
	}
}

/* Scale alpha values */
__global__ void scale_alpha_dev( int nstates,             /* number of states */  
		float *alpha_d, /* dim = length x nstates */  
		int offset,              /* offset for alpha_d */   
		float scale)             /* scaling value */    
{
	// unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < nstates)
	{
		alpha_d[offset + idx] = alpha_d[offset + idx] / scale;
	}
}

/* Initialize beta values */
__global__ void init_beta_dev(  int nstates,            /* number of states */  
		float *beta_d, /* dim = length x nstates */   
		int offset,             /* offset for beta_d */   
		float scale)            /* scaling value */    
{
	// unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < nstates)
	{
		beta_d[offset + idx] = 1.0f / scale;
	}
}

/* Calculate beta variables */
__global__ void calc_beta_dev( float *beta_d, /* dim = length x nstates */ 
		float *b_d,    /* dim = nsymbols x nstates */ 
		float scale_t,          /* current scaling value */  
		int nstates,            /* number of states */   
		int obs_t,              /* current observation */    
		int t)                  /* current time */    
{
	// unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < nstates)
	{
		beta_d[(t * nstates) + idx] = beta_d[((t + 1) * nstates) + idx] *
			b_d[(obs_t * nstates) + idx] / scale_t;
	}
}

/* Sum next iteration of gamma variables */
__global__ void calc_gamma_dev(float *gamma_sum_d, /* lth = nstates */  
		float *alpha_d,     /* dim = length x nstates */  
		float *beta_d,      /* dim = length x nstates */    
		int nstates,                 /* number of states */   
		int t)                       /* current time */  
{
	// unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < nstates)
	{
		gamma_sum_d[idx] += alpha_d[(t * nstates) + idx] *
			beta_d[(t * nstates) + idx];
	}
}

/* Sum next iteration of xi variables */
__global__ void calc_xi_dev(   float *xi_sum_d,  /* dim = nstates x nstates */   
		float *a_d,       /* dim = nstates x nstates */  
		float *b_d,       /* dim = nsymbols x nstates */  
		float *alpha_d,   /* dim = length x nstates */    
		float *beta_d,    /* dim = length x nstates */   
		float sum_ab,              /* sum dot production of (alpha_d+t*nstates) and (beta_d+t*nstates) */  
		int nstates,               /* number of states */     
		int obs_t,                 /* next observation at t + 1 */    
		int t)                     /* current time */   
{
	// unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	// unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;
	unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;

	if (idx < nstates && idy < nstates)
	{
		xi_sum_d[(idy * nstates) + idx] += alpha_d[(t * nstates) + idy] *
			a_d[(idy * nstates) + idx] *
			b_d[(obs_t * nstates) + idx] *
			beta_d[((t+1) * nstates) + idx] /
			sum_ab;
	}
}

/* Re-estimate A matrix */
__global__ void est_a_dev(  float *a_d,          /* dim = nstates x nstates */ 
		float *alpha_d,      /* dim = length x nstates */      
		float *beta_d,       /* dim = length x nstates */   
		float *xi_sum_d,     /* dim = nstates x nstates */  
		float *gamma_sum_d,  /* lth = nstates */  
		float sum_ab,                 /* sum of dot production of (alpha_d+(length-1)*nstates) and  */ 
		/* (beta_d+(length-1)*nstates)  */ 
		int nstates,                  /* number of states */  
		int length)                   /* number of observations */    
{
	// unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	// unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;
	unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;

	if (idx < nstates && idy < nstates)
	{
		a_d[(idy * nstates) + idx] = xi_sum_d[(idy * nstates) + idx] /
			(gamma_sum_d[idy] -
			 alpha_d[(length * nstates) + idy] *
			 beta_d[(length * nstates) + idy] /
			 sum_ab);
	}
}

/* Normalize A matrix */
__global__ void scale_a_dev(  float *a_d, /* dim = nstates x nstates */  
		float *c_d, /* lth = nstates, matrix-vector product of a_d and ones_n_d */  
		int nstates)          /* number of states */    
{
	// unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	// unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;
	unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;

	if (idx < nstates && idy < nstates)
	{
		a_d[(idy * nstates) + idx] = a_d[(idy * nstates) + idx] / c_d[idy];
	}
}

/* Accumulate B values */
__global__ void acc_b_dev( float *b_d,      /* dim = nsymbols x nstates */  
		float *alpha_d,  /* dim = length x nstates */   
		float *beta_d,   /* dim = length x nstates */  
		float sum_ab,             /* sum dot production of (alpha_d+t*nstates) and (beta_d+t*nstates) */   
		int nstates,              /* number of states */    
		int nsymbols,             /* number of symbols */  
		int obs_t,                /* current observation */        
		int t)                    /* current time */   
{
	// unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	// unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;
	unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;

	if (idy < nsymbols && idx < nstates && obs_t == idy)
	{
		b_d[(idy * nstates) + idx] += alpha_d[(t * nstates) + idx] *
			beta_d[(t * nstates) + idx] / sum_ab;
	}
}

/* Re-estimate B values */
__global__ void est_b_dev( float *b_d,         /* dim = nsymbols x nstates */ 
		float *gamma_sum_d, /* lth = nstates */  
		int nstates,                 /* number of states */  
		int nsymbols)                /* number of symbols */  
{
	// unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	// unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;
	unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;

	if (idy < nsymbols && idx < nstates)
	{
		b_d[(idy * nstates) + idx] = b_d[(idy * nstates) + idx] /
			gamma_sum_d[idx];
	}
}

/* Normalize B matrix */
__global__ void scale_b_dev(  float *b_d, /* dim = nsymbols x nstates */   
		float *c_d, /* lth = nstates, matrix-vector product of b_d and ones_s_d */    
		int nstates,         /* number of states */     
		int nsymbols)        /* number of symbols */     
{
	// unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	// unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;
	unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;

	if (idx < nstates && idy < nsymbols)
	{
		if (b_d[(idy * nstates) + idx] == 0)
		{
			b_d[(idy * nstates) + idx] = 1e-10;
		}
		else
		{
			b_d[(idy * nstates) + idx] = b_d[(idy * nstates) + idx] / c_d[idx];
		}
	}
}

/* Re-estimate Pi values */
__global__ void est_pi_dev(float *pi_d,    /* lth = nstates */  
		float *alpha_d, /* dim = length x nstates */   
		float *beta_d,  /* dim = length x nstates */   
		float sum_ab,            /* sum dot production of alpha_d and beta_d */  
		int nstates)             /* number of states */  
{
	// unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
	unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;

	if (idx < nstates)
	{
		pi_d[idx] = alpha_d[idx] * beta_d[idx] / sum_ab;
	}
}

/* dot production naive version */
__global__ void s_dot_kernel_naive(int n,                        /* number of elements */ 
		float *paramA,       /* first vector A */  
		int offsetA,                  /* offset of vector A */ 
		float *paramB,       /* second vector B */ 
		int offsetB,                  /* offset of vector B */    
		float *partialSum_d) /* memory space for partial sum */  
{
	unsigned int i;
	unsigned int tid = threadIdx.x;
	unsigned int totalThreads = gridDim.x * SDOT_BLOCK_SIZE;
	unsigned int offset = SDOT_BLOCK_SIZE * blockIdx.x;

	/* [FIX-OCL-34] pinned launch (global = roundup(n,256), local 256): grid-stride over the global range (one element per work-item);
	   was offset/stride built from SDOT_BLOCK_SIZE=128 / 80 blocks */
	for(i = blockIdx.x * blockDim.x + threadIdx.x; i < n; i += gridDim.x * blockDim.x)
	{ 
		partialSum_d[i] = paramA[offsetA + i] * paramB[offsetB + i];
	}
}

/* mat x vec naive verion */
__global__ void mvm_non_kernel_naive(int m,
		int n,
		float *A,
		int lda,
		float *x,
		int offsetX,
		float *y,
		int offsetY)
{
	unsigned int i, j;
	unsigned int tid = threadIdx.x;
	unsigned int totalThreads = gridDim.x * MVMUL_BLOCK_SIZE;
	unsigned int offset = MVMUL_BLOCK_SIZE * blockIdx.x;
	int n_size, m_size;

	float sum;
	if(lda == m)
	{
		n_size = n;
		m_size = m;
	} else 
	{
		n_size = m;
		m_size = n;
	}

	/* [FIX-OCL-34] pinned launch (global = roundup(m_size,256), local 256): one row per work-item (was 128 x 64 grid-stride) */
	for(i = blockIdx.x * blockDim.x + threadIdx.x; i < m_size; i += gridDim.x * blockDim.x)
	{
		sum = 0.f;
		for(j = 0; j < n_size; j++)
		{
			sum += A[i * n_size + j] * x[j+offsetX];
		}
		y[i+offsetY] = sum;
	}

}

/* mat x vec transposed naive verion */
__global__ void mvm_trans_kernel_naive(int m,
		int n,
		float *A,
		int lda,
		float *x,
		int offsetX,
		float *y,
		int offsetY)
{
	unsigned int i, j;
	unsigned int tid = threadIdx.x;
	unsigned int totalThreads = gridDim.x * MVMUL_BLOCK_SIZE;
	unsigned int offset = MVMUL_BLOCK_SIZE * blockIdx.x;
	int n_size, m_size;

	float sum;
	if(lda == m)
	{
		n_size = n;
		m_size = m;
	} else 
	{
		n_size = m;
		m_size = n;
	}

	/* [FIX-OCL-34] pinned launch (global = roundup(m_size,256), local 256): one row per work-item (was 128 x 64 grid-stride) */
	for(i = blockIdx.x * blockDim.x + threadIdx.x; i < m_size; i += gridDim.x * blockDim.x)
	{
		sum = 0.f;
		for(j = 0; j < n_size; j++)
		{
			sum += A[j * n_size + i] * x[j+offsetX];
		}
		y[i+offsetY] = sum;
	}

}

