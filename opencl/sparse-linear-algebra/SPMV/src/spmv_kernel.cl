/*
 */

void __kernel csr(const unsigned int num_rows,
		__global unsigned int * Ap, 
		__global unsigned int * Aj, 
		__global float * Ax, 
		__global float * x, 
		__global float * y)
{
	unsigned int row = get_global_id(0);

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
