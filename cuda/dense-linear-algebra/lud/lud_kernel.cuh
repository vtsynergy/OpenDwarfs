/* CUDA translation of opencl/dense-linear-algebra/lud/lud_kernel.cl (block shapes = OpenCL work-groups: 16, 32, 16x16). */
#ifndef BLOCK_SIZE_K
#define BLOCK_SIZE_K 16
#endif

__global__ void lud_diagonal(float *m, int matrix_dim, int offset)
{
	int i,j;
	int array_offset = offset*matrix_dim+offset;
	int lid = threadIdx.x;

	for(i=0; i < BLOCK_SIZE_K; i++) {
		if (lid>i){
			for(j=0; j < i; j++)
				m[array_offset+lid*matrix_dim+i] -= m[array_offset+lid*matrix_dim+j]*m[array_offset+j*matrix_dim+i];
			m[array_offset+lid*matrix_dim+i] /= m[array_offset+i*matrix_dim+i];
		}
		__syncthreads();
		if (lid>i){
			for(j=0; j < i+1; j++)
				m[array_offset+(i+1)*matrix_dim+lid] -= m[array_offset+(i+1)*matrix_dim+j]*m[array_offset+j*matrix_dim+lid];
		}
		__syncthreads();
	}
}

__global__ void lud_perimeter(float *m, int matrix_dim, int offset)
{
	int i,j, array_offset;
	int array_offset_peri;
	int idx;
	int lid = threadIdx.x;
	int gid = blockIdx.x;

	if (lid < BLOCK_SIZE_K) { //peri-row
		idx=lid;
		array_offset = offset*matrix_dim+offset;
		for(i=1; i < BLOCK_SIZE_K; i++){
/* [CUDA-WA-01] nvcc 12.9 (V12.9.86) -O3 miscompiles this loop when it unrolls it: the stores of the unrolled copy are emitted as
 * [address already containing the +BLOCK_SIZE column offset + BLOCK_SIZE], i.e. 16 columns too far; in the last block they wrap
 * into the diagonal block. Seen in the PTX; -Xcicc -O0 or no unrolling gives correct results. Disable unrolling here only. */
#pragma unroll 1
			for (j=0; j < i; j++)
				m[array_offset+i*matrix_dim+(gid+1)*BLOCK_SIZE_K+idx] -= m[array_offset+i*matrix_dim+j]*m[array_offset+j*matrix_dim+(gid+1)*BLOCK_SIZE_K+idx];
		}
	} else { //peri-col
		idx=lid - BLOCK_SIZE_K;
		array_offset_peri = (offset+(gid+1)*BLOCK_SIZE_K)*matrix_dim+offset;
		array_offset = offset*matrix_dim+offset;
		for(i=0; i < BLOCK_SIZE_K; i++){
			for(j=0; j < i; j++)
				m[array_offset_peri+idx*matrix_dim+i] -= m[array_offset_peri+idx*matrix_dim+j] * m[array_offset+j*matrix_dim+i];
			m[array_offset_peri+idx*matrix_dim+i] /= m[array_offset+i*matrix_dim+i];
		}
	}
}

__global__ void lud_internal(float *m, int matrix_dim, int offset)
{
	int i;
	float sum;

	int global_row_id = offset + ((int)blockIdx.y+1)*BLOCK_SIZE_K;
	int global_col_id = offset + ((int)blockIdx.x+1)*BLOCK_SIZE_K;

	sum = 0;
	for (i=0; i < BLOCK_SIZE_K; i++)
		sum += m[(global_row_id+(int)threadIdx.y)*matrix_dim+offset+i] * m[(offset+i)*matrix_dim+global_col_id+(int)threadIdx.x];
	m[(global_row_id+(int)threadIdx.y)*matrix_dim+global_col_id+(int)threadIdx.x] -= sum;
}
