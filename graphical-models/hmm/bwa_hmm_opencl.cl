/************************************************************************/
/* OpenCL Kernels for BW Algorithm in Hidden Markov Model               */
/*                                                                      */
/************************************************************************/
#define SDOT_BLOCK_SIZE        (128)
#define SDOT_BLOCK_NUM         (80)

#define MVMUL_BLOCK_SIZE       (128)
#define MVMUL_BLOCK_NUM        (64)
#define TILE_SIZE              (32)

/* Initialize ones vector */
__kernel void init_ones_dev(  __global float *ones_s_d,  /* lth = nsymbols */    
                              int nsymbols)              /* number of symbols */   
{
    //unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int idx = get_group_id(0) * get_local_size(0) + get_local_id(0);

    if (idx < nsymbols)
    {
        ones_s_d[idx] = 1.0f;
    }
}

/* Initialize alpha variables */
__kernel void init_alpha_dev( __global float *b_d,      /* dim = nsymbols x nstates */   
                              __global  float *pi_d,    /* lth = nstates */         
                              int nstates,              /* number of states */      
                              __global float *alpha_d,  /* dim = length x nstates */     
                              __global float *ones_n_d, /* lth = nstates */       
                              int obs_t)                /* current observation */     
{
    // unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int idx = get_group_id(0) * get_local_size(0) + get_local_id(0);

    if (idx < nstates)
    {
        alpha_d[idx] = pi_d[idx] * b_d[(obs_t * nstates) + idx];
        ones_n_d[idx] = 1.0f;
    }
}

/* Calculate alpha variables */
__kernel void calc_alpha_dev( int nstates,             /* number of states */     
                              __global float *alpha_d, /* dim = length x nstates */        
                              int offset,              /* offset for alpha_d */  
                              __global float *b_d,     /* dim = nsymbols x nstates */  
                              int obs_t)               /* current observation */     
{
    // unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int idx = get_group_id(0) * get_local_size(0) + get_local_id(0);

    if (idx < nstates)
    {
        alpha_d[offset + idx] = alpha_d[offset + idx] * b_d[(obs_t * nstates) + idx];
    }
}

/* Scale alpha values */
__kernel void scale_alpha_dev( int nstates,             /* number of states */  
                               __global float *alpha_d, /* dim = length x nstates */  
                               int offset,              /* offset for alpha_d */   
                               float scale)             /* scaling value */    
{
    // unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int idx = get_group_id(0) * get_local_size(0) + get_local_id(0);

    if (idx < nstates)
    {
        alpha_d[offset + idx] = alpha_d[offset + idx] / scale;
    }
}

/* Initialize beta values */
__kernel void init_beta_dev(  int nstates,            /* number of states */  
                              __global float *beta_d, /* dim = length x nstates */   
                              int offset,             /* offset for beta_d */   
                              float scale)            /* scaling value */    
{
    // unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int idx = get_group_id(0) * get_local_size(0) + get_local_id(0);

    if (idx < nstates)
    {
        beta_d[offset + idx] = 1.0f / scale;
    }
}

/* Calculate beta variables */
__kernel void calc_beta_dev( __global float *beta_d, /* dim = length x nstates */ 
                             __global float *b_d,    /* dim = nsymbols x nstates */ 
                             float scale_t,          /* current scaling value */  
                             int nstates,            /* number of states */   
                             int obs_t,              /* current observation */    
                             int t)                  /* current time */    
{
    // unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int idx = get_group_id(0) * get_local_size(0) + get_local_id(0);

    if (idx < nstates)
    {
        beta_d[(t * nstates) + idx] = beta_d[((t + 1) * nstates) + idx] *
                                      b_d[(obs_t * nstates) + idx] / scale_t;
    }
}

/* Sum next iteration of gamma variables */
__kernel void calc_gamma_dev(__global float *gamma_sum_d, /* lth = nstates */  
                             __global float *alpha_d,     /* dim = length x nstates */  
                             __global float *beta_d,      /* dim = length x nstates */    
                             int nstates,                 /* number of states */   
                             int t)                       /* current time */  
{
    // unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int idx = get_group_id(0) * get_local_size(0) + get_local_id(0);

    if (idx < nstates)
    {
        gamma_sum_d[idx] += alpha_d[(t * nstates) + idx] *
                            beta_d[(t * nstates) + idx];
    }
}

/* Sum next iteration of xi variables */
__kernel void calc_xi_dev(   __global float *xi_sum_d,  /* dim = nstates x nstates */   
                             __global float *a_d,       /* dim = nstates x nstates */  
                             __global float *b_d,       /* dim = nsymbols x nstates */  
                             __global float *alpha_d,   /* dim = length x nstates */    
                             __global float *beta_d,    /* dim = length x nstates */   
                             float sum_ab,              /* sum dot production of (alpha_d+t*nstates) and (beta_d+t*nstates) */  
                             int nstates,               /* number of states */     
                             int obs_t,                 /* next observation at t + 1 */    
                             int t)                     /* current time */   
{
    // unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int idx = get_group_id(0) * get_local_size(0) + get_local_id(0);
    // unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned int idy = get_group_id(1) * get_local_size(1) + get_local_id(1);

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
__kernel void est_a_dev(  __global float *a_d,          /* dim = nstates x nstates */ 
                          __global float *alpha_d,      /* dim = length x nstates */      
                          __global float *beta_d,       /* dim = length x nstates */   
                          __global float *xi_sum_d,     /* dim = nstates x nstates */  
                          __global float *gamma_sum_d,  /* lth = nstates */  
                          float sum_ab,                 /* sum of dot production of (alpha_d+(length-1)*nstates) and  */ 
                                                        /* (beta_d+(length-1)*nstates)  */ 
                          int nstates,                  /* number of states */  
                          int length)                   /* number of observations */    
{
    // unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int idx = get_group_id(0) * get_local_size(0) + get_local_id(0);
    // unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned int idy = get_group_id(1) * get_local_size(1) + get_local_id(1);

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
__kernel void scale_a_dev(  __global  float *a_d, /* dim = nstates x nstates */  
                            __global  float *c_d, /* lth = nstates, matrix-vector product of a_d and ones_n_d */  
                            int nstates)          /* number of states */    
{
    // unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int idx = get_group_id(0) * get_local_size(0) + get_local_id(0);
    // unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned int idy = get_group_id(1) * get_local_size(1) + get_local_id(1);

    if (idx < nstates && idy < nstates)
    {
        a_d[(idy * nstates) + idx] = a_d[(idy * nstates) + idx] / c_d[idy];
    }
}

/* Accumulate B values */
__kernel void acc_b_dev( __global float *b_d,      /* dim = nsymbols x nstates */  
                         __global float *alpha_d,  /* dim = length x nstates */   
                         __global float *beta_d,   /* dim = length x nstates */  
                         float sum_ab,             /* sum dot production of (alpha_d+t*nstates) and (beta_d+t*nstates) */   
                         int nstates,              /* number of states */    
                         int nsymbols,             /* number of symbols */  
                         int obs_t,                /* current observation */        
                         int t)                    /* current time */   
{
    // unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int idx = get_group_id(0) * get_local_size(0) + get_local_id(0);
    // unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned int idy = get_group_id(1) * get_local_size(1) + get_local_id(1);

    if (idy < nsymbols && idx < nstates && obs_t == idy)
    {
        b_d[(idy * nstates) + idx] += alpha_d[(t * nstates) + idx] *
                                      beta_d[(t * nstates) + idx] / sum_ab;
    }
}

/* Re-estimate B values */
__kernel void est_b_dev( __global float *b_d,         /* dim = nsymbols x nstates */ 
                         __global float *gamma_sum_d, /* lth = nstates */  
                         int nstates,                 /* number of states */  
                         int nsymbols)                /* number of symbols */  
{
    // unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int idx = get_group_id(0) * get_local_size(0) + get_local_id(0);
    // unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned int idy = get_group_id(1) * get_local_size(1) + get_local_id(1);

    if (idy < nsymbols && idx < nstates)
    {
        b_d[(idy * nstates) + idx] = b_d[(idy * nstates) + idx] /
                                     gamma_sum_d[idx];
    }
}

/* Normalize B matrix */
__kernel void scale_b_dev(  __global float *b_d, /* dim = nsymbols x nstates */   
                            __global float *c_d, /* lth = nstates, matrix-vector product of b_d and ones_s_d */    
                            int nstates,         /* number of states */     
                            int nsymbols)        /* number of symbols */     
{
    // unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int idx = get_group_id(0) * get_local_size(0) + get_local_id(0);
    // unsigned int idy = blockIdx.y * blockDim.y + threadIdx.y;
    unsigned int idy = get_group_id(1) * get_local_size(1) + get_local_id(1);

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
__kernel void est_pi_dev(__global float *pi_d,    /* lth = nstates */  
                         __global float *alpha_d, /* dim = length x nstates */   
                         __global float *beta_d,  /* dim = length x nstates */   
                         float sum_ab,            /* sum dot production of alpha_d and beta_d */  
                         int nstates)             /* number of states */  
{
    // unsigned int idx = blockIdx.x * blockDim.x + threadIdx.x;
    unsigned int idx = get_group_id(0) * get_local_size(0) + get_local_id(0);

    if (idx < nstates)
    {
        pi_d[idx] = alpha_d[idx] * beta_d[idx] / sum_ab;
    }
}

/* dot production naive version */
__kernel void s_dot_kernel_naive(int n,                        /* number of elements */ 
                           __global float *paramA,       /* first vector A */  
                           int offsetA,                  /* offset of vector A */ 
                           __global float *paramB,       /* second vector B */ 
                           int offsetB,                  /* offset of vector B */    
                           __global float *partialSum_d) /* memory space for partial sum */  
{
    unsigned int i;
    unsigned int tid = get_local_id(0);
    unsigned int totalThreads = get_num_groups(0) * SDOT_BLOCK_SIZE;
    unsigned int offset = SDOT_BLOCK_SIZE * get_group_id(0);

    for(i = offset + tid; i < n; i += totalThreads)
    { 
        partialSum_d[i] = paramA[offsetA + i] * paramB[offsetB + i];
    }
}

/* mat x vec naive verion */
__kernel void mvm_non_kernel_naive(int m,
                            int n,
                            __global float *A,
                            int lda,
                            __global float *x,
                            int offsetX,
                            __global float *y,
                            int offsetY)
{
    unsigned int i, j;
    unsigned int tid = get_local_id(0);
    unsigned int totalThreads = get_num_groups(0) * MVMUL_BLOCK_SIZE;
    unsigned int offset = MVMUL_BLOCK_SIZE * get_group_id(0);
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

    for(i = offset + tid; i < m_size; i += totalThreads)
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
__kernel void mvm_trans_kernel_naive(int m,
                            int n,
                            __global float *A,
                            int lda,
                            __global float *x,
                            int offsetX,
                            __global float *y,
                            int offsetY)
{
    unsigned int i, j;
    unsigned int tid = get_local_id(0);
    unsigned int totalThreads = get_num_groups(0) * MVMUL_BLOCK_SIZE;
    unsigned int offset = MVMUL_BLOCK_SIZE * get_group_id(0);
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

    for(i = offset + tid; i < m_size; i += totalThreads)
    {
        sum = 0.f;
        for(j = 0; j < n_size; j++)
        {
            sum += A[j * n_size + i] * x[j+offsetX];
        }
        y[i+offsetY] = sum;
    }

}

