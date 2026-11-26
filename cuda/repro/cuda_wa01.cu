/* CUDA-WA-01 minimal reproducer: nvcc 12.9 (V12.9.86), -O3, sm_90.
 * lud_perimeter() is the kernel of ../dense-linear-algebra/lud/lud_kernel.cuh without the workaround (-DWORKAROUND restores it).
 * A correct kernel NEVER writes the 16x16 diagonal block: peri-row threads write rows 1..15 at columns >= 16, peri-col threads
 * write rows >= 16 at columns 0..15. The test also compares the whole result with a serial host version of the kernel. */
#include <cuda_runtime.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define BLOCK_SIZE_K 16
#define DIM 256

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
#ifdef WORKAROUND
#pragma unroll 1
#endif
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

int main()
{
	static float in[DIM*DIM], ref[DIM*DIM], out[DIM*DIM];
	const int B = BLOCK_SIZE_K, nb = DIM/B - 1;
	unsigned s = 12345;
	for (int i = 0; i < DIM*DIM; i++) { s = s*1103515245u + 12345u; in[i] = 0.1f + (s >> 16)/65536.0f*0.1f; ref[i] = 0; }
	for (int i = 0; i < DIM; i++) in[i*DIM+i] += 4.0f;
	for (int i = 0; i < DIM*DIM; i++) ref[i] = in[i];
	for (int i = 1; i < B; i++) for (int j = 0; j < i; j++) for (int c = B; c < DIM; c++) ref[i*DIM+c] -= ref[i*DIM+j]*ref[j*DIM+c];
	for (int r = B; r < DIM; r++) for (int i = 0; i < B; i++) { for (int j = 0; j < i; j++) ref[r*DIM+i] -= ref[r*DIM+j]*ref[j*DIM+i]; ref[r*DIM+i] /= ref[i*DIM+i]; }
	float *d;
	if (cudaMalloc((void**)&d, sizeof(in)) != cudaSuccess) { puts("no device"); return 2; }
	cudaMemcpy(d, in, sizeof(in), cudaMemcpyHostToDevice);
	lud_perimeter<<<nb, 2*BLOCK_SIZE_K>>>(d, DIM, 0);
	cudaError_t e = cudaDeviceSynchronize();
	cudaMemcpy(out, d, sizeof(out), cudaMemcpyDeviceToHost);
	int diag = 0, bad = 0;
	for (int r = 0; r < DIM; r++) for (int c = 0; c < DIM; c++) {
		if (r < B && c < B && out[r*DIM+c] != in[r*DIM+c]) { if (!diag) printf("first diagonal-block element written: row %d col %d\n", r, c); diag++; }
		if (!(fabsf(out[r*DIM+c]-ref[r*DIM+c]) <= 1e-3f*fmaxf(1.0f, fabsf(ref[r*DIM+c])))) bad++;
	}
	printf("%s: %d diagonal-block elements modified, %d elements differ from the serial reference (%s)\n",
	       diag || bad ? "FAIL" : "PASS", diag, bad, cudaGetErrorString(e));
	return diag || bad;
}
