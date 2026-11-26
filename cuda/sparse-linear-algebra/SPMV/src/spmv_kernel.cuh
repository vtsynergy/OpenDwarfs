/* CUDA translation of src/spmv_kernel.cl: one matrix row per thread; DWARF_LS threads/block, guard row < num_rows. */
__global__ void csr_kernel(const unsigned int num_rows,
		const unsigned int * Ap,
		const unsigned int * Aj,
		const float * Ax,
		const float * x,
		float * y)
{
	unsigned int row = blockIdx.x * blockDim.x + threadIdx.x;

	if(row < num_rows)
	{
		float sum = y[row];

		const unsigned int row_start = Ap[row];
		const unsigned int row_end = Ap[row+1];

		unsigned int jj = 0;
		for (jj = row_start; jj < row_end; jj++)
			sum += Ax[jj] * x[Aj[jj]];

		y[row] = sum;
	}
}
