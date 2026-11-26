/* CUDA translation of needle_kernel.cl: one 16x16 tile per block of BLOCK_SIZE (16) threads, one barrier per anti-diagonal. */
__device__ static inline int maximum_dev(int a, int b, int c)
{
	int k;
	if( a <= b )
		k = b;
	else
		k = a;
	if( k <=c )
		return(c);
	else
		return(k);
}

__global__ void needle_opencl_shared_1(int* referrence, int* matrix_opencl, int cols, int penalty, int i, int block_width)
{
	int bx = blockIdx.x;
	int tx = threadIdx.x;

	int b_index_x = bx;
	int b_index_y = i - 1 - bx;
	int index   = cols * BLOCK_SIZE * b_index_y + BLOCK_SIZE * b_index_x + tx + ( cols + 1 );

	for( int m = 0 ; m < BLOCK_SIZE ; m++){
		if ( tx <= m ){
			int ref_x=index+(m-tx)*cols;
			matrix_opencl[ref_x] = maximum_dev( matrix_opencl[ref_x-(cols+1)] + referrence[ref_x],
					matrix_opencl[ref_x-1]  - penalty,
					matrix_opencl[ref_x-cols]  - penalty);
		}
		__syncthreads();
	}

	for( int m = BLOCK_SIZE - 2 ; m >=0 ; m--){
		if ( tx <= m){
			int ref_x=index+(m-tx)*cols+(cols+1)*(BLOCK_SIZE-1-m);
			matrix_opencl[ref_x] = maximum_dev( matrix_opencl[ref_x-(cols+1)] + referrence[ref_x],
					matrix_opencl[ref_x-1]  - penalty,
					matrix_opencl[ref_x-cols]  - penalty);
		}
		__syncthreads();
	}
}

__global__ void needle_opencl_shared_2(int* referrence, int* matrix_opencl, int cols, int penalty, int i, int block_width)
{
	int bx = blockIdx.x;
	int tx = threadIdx.x;

	int b_index_x = bx + block_width - i  ;
	int b_index_y = block_width - bx -1;
	int index   = cols * BLOCK_SIZE * b_index_y + BLOCK_SIZE * b_index_x + tx + ( cols + 1 );

	for( int m = 0 ; m < BLOCK_SIZE ; m++){
		if ( tx <= m ){
			int ref_x=index+(m-tx)*cols;
			matrix_opencl[ref_x] = maximum_dev( matrix_opencl[ref_x-(cols+1)] + referrence[ref_x],
					matrix_opencl[ref_x-1]  - penalty,
					matrix_opencl[ref_x-cols]  - penalty);
		}
		__syncthreads();
	}

	for( int m = BLOCK_SIZE - 2 ; m >=0 ; m--){
		if ( tx <= m){
			int ref_x=index+(m-tx)*cols+(cols+1)*(BLOCK_SIZE-1-m);
			matrix_opencl[ref_x] = maximum_dev( matrix_opencl[ref_x-(cols+1)] + referrence[ref_x],
					matrix_opencl[ref_x-1]  - penalty,
					matrix_opencl[ref_x-cols]  - penalty);
		}
		__syncthreads();
	}
}
