// includes, system
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "srad.h"
#include <string.h>
#include <malloc.h>
//#define OUTPUT
#define AOCL_ALIGNMENT 64
#define BLOCK_SIZE 16   /* srad_kernel.cl */
#include "dwarf_sycl.hpp"
#include "dwarf_timer.h"
static int srad_valid = 1;
/* [FIX-OCL-09] The app had no correctness check, and its (disabled, #ifdef CPU) host path declared the diffusion
 * coefficient array as int (`int *c`), truncating every coefficient to 0 or 1. This independent serial float
 * implementation (same equations, clamped borders) replaces it as the reference for the device result. */
static void srad_reference(float *J, int rows, int cols, int r1, int r2, int c1, int c2, float lambda, int niter)
{
	int size_I = rows * cols, size_R = (r2 - r1 + 1) * (c2 - c1 + 1), i, j, k, iter;
	float *dN = (float*)malloc(sizeof(float)*size_I), *dS = (float*)malloc(sizeof(float)*size_I),
	      *dW = (float*)malloc(sizeof(float)*size_I), *dE = (float*)malloc(sizeof(float)*size_I), *cf = (float*)malloc(sizeof(float)*size_I);
	for (iter = 0; iter < niter; iter++) {
		float sum = 0, sum2 = 0, tmp, meanROI, varROI, q0sqr;
		for (i = r1; i <= r2; i++) for (j = c1; j <= c2; j++) { tmp = J[i*cols+j]; sum += tmp; sum2 += tmp*tmp; }
		meanROI = sum / size_R; varROI = (sum2 / size_R) - meanROI*meanROI; q0sqr = varROI / (meanROI*meanROI);
		for (i = 0; i < rows; i++) for (j = 0; j < cols; j++) {
			float Jc, G2, L, num, den, qsqr, c;
			int iN = i == 0 ? 0 : i-1, iS = i == rows-1 ? rows-1 : i+1, jW = j == 0 ? 0 : j-1, jE = j == cols-1 ? cols-1 : j+1;
			k = i*cols + j; Jc = J[k];
			dN[k] = J[iN*cols+j] - Jc; dS[k] = J[iS*cols+j] - Jc; dW[k] = J[i*cols+jW] - Jc; dE[k] = J[i*cols+jE] - Jc;
			G2 = (dN[k]*dN[k] + dS[k]*dS[k] + dW[k]*dW[k] + dE[k]*dE[k]) / (Jc*Jc);
			L = (dN[k] + dS[k] + dW[k] + dE[k]) / Jc;
			num = (0.5f*G2) - ((1.0f/16.0f)*(L*L)); den = 1 + (.25f*L); qsqr = num/(den*den);
			den = (qsqr-q0sqr) / (q0sqr * (1+q0sqr)); c = 1.0f / (1.0f+den);
			cf[k] = c < 0 ? 0 : (c > 1 ? 1 : c);
		}
		for (i = 0; i < rows; i++) for (j = 0; j < cols; j++) {
			float cN, cS, cW, cE, D; int iS = i == rows-1 ? rows-1 : i+1, jE = j == cols-1 ? cols-1 : j+1;
			k = i*cols + j; cN = cf[k]; cS = cf[iS*cols+j]; cW = cf[k]; cE = cf[i*cols+jE];
			D = cN*dN[k] + cS*dS[k] + cW*dW[k] + cE*dE[k];
			J[k] = J[k] + 0.25f*lambda*D;
		}
	}
	free(dN); free(dS); free(dW); free(dE); free(cf);
}

//If GPU is defined in srad.h the C implementation runs, otherwise the OpenCL one.

/* kernels (translated from srad_kernel.cl). OpenCL NDRange {cols, rows} local {16, 16}: OpenCL dim 0 (x) = SYCL dim 1. */
static void srad_cuda_1(sycl::queue &q, sycl::nd_range<2> ndr, float *E_C, float *W_C, float *N_C, float *S_C,
		float *J_cuda, float *C_cuda, int cols, int rows, float q0sqr)
{
	q.parallel_for<class srad_k1>(ndr, [=](sycl::nd_item<2> it) {

		//block id
		int bx = (int)it.get_group(1);
		int by = (int)it.get_group(0);

		//thread id
		int tx = (int)it.get_local_id(1);
		int ty = (int)it.get_local_id(0);

		//indices
		int index   = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + tx;
		int index_n = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + tx - cols;   //up
		int index_s = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * BLOCK_SIZE + tx; //down
		int index_w = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty - 1;  //left
		int index_e = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + BLOCK_SIZE;  //right

		float n, w, e, s, jc, g2, l, num, den, qsqr, c;


		if ( by == 0 ){
			index_n = BLOCK_SIZE * bx + tx;
		}
		else if ( by == (int)it.get_group_range(0) - 1 ){
			index_s = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ( BLOCK_SIZE - 1 ) + tx;
		}

		if ( bx == 0 ){
			index_w = cols * BLOCK_SIZE * by + cols * ty;
		}
		else if ( bx == (int)it.get_group_range(1) - 1 ){
			index_e = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + BLOCK_SIZE-1;
		}

		//Corresponding to [ty+1][x], [ty-1][x], etc.
		int index_yp1_x = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * (ty + 1) + tx;
		int index_ym1_x = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * (ty - 1) + tx;
		int index_y_xp1 = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols  * ty + (tx + 1);
		int index_y_xm1 = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + (tx - 1);


		jc = J_cuda[index];

		if ( ty == 0 && tx == 0 ){ //nw
			n  = J_cuda[index_n] - jc;
			s  = J_cuda[index_yp1_x] - jc;
			w  = J_cuda[index_w]  - jc; 
			e  = J_cuda[index_y_xp1] - jc;
		}	    
		else if ( ty == 0 && tx == BLOCK_SIZE-1 ){ //ne
			n  = J_cuda[index_n] - jc;
			s  = J_cuda[index_yp1_x] - jc;
			w  = J_cuda[index_y_xm1] - jc;
			e  = J_cuda[index_e] - jc;
		}
		else if ( ty == BLOCK_SIZE -1 && tx == BLOCK_SIZE - 1){ //se
			n  = J_cuda[index_ym1_x] - jc;
			s  = J_cuda[index_s] - jc;
			w  = J_cuda[index_y_xm1] - jc; 
			e  = J_cuda[index_e]  - jc;
		}
		else if ( ty == BLOCK_SIZE -1 && tx == 0 ){//sw
			n  = J_cuda[index_ym1_x] - jc;
			s  = J_cuda[index_s] - jc;
			w  = J_cuda[index_w]  - jc; 
			e  = J_cuda[index_y_xp1] - jc;
		}
		else if ( ty == 0 ){ //n
			n  = J_cuda[index_n] - jc;
			s  = J_cuda[index_yp1_x] - jc;
			w  = J_cuda[index_y_xm1] - jc; 
			e  = J_cuda[index_y_xp1] - jc;
		}
		else if ( tx == BLOCK_SIZE -1 ){ //e
			n  = J_cuda[index_ym1_x] - jc;
			s  = J_cuda[index_yp1_x] - jc;
			w  = J_cuda[index_y_xm1] - jc; 
			e  = J_cuda[index_e] - jc;
		}
		else if ( ty == BLOCK_SIZE -1){ //s
			n  = J_cuda[index_ym1_x] - jc;
			s  = J_cuda[index_s] - jc;
			w  = J_cuda[index_y_xm1] - jc; 
			e  = J_cuda[index_y_xp1] - jc;
		}
		else if ( tx == 0 ){ //w
			n  = J_cuda[index_ym1_x] - jc;
			s  = J_cuda[index_yp1_x] - jc;
			w  = J_cuda[index_w] - jc; 
			e  = J_cuda[index_y_xp1] - jc;
		}
		else{  //the data elements which are not on the borders 
			n  = J_cuda[index_ym1_x] - jc;
			s  = J_cuda[index_yp1_x] - jc;
			w  = J_cuda[index_y_xm1] - jc; 
			e  = J_cuda[index_y_xp1] - jc;
		}


		g2 = ( n * n + s * s + w * w + e * e ) / (jc * jc);

		l = ( n + s + w + e ) / jc;

		num  = (0.5f*g2) - ((1.0f/16.0f)*(l*l)) ;
		den  = 1 + (.25f*l);
		qsqr = num/(den*den);

		// diffusion coefficent (equ 33)
		den = (qsqr-q0sqr) / (q0sqr * (1+q0sqr)) ;
		c = 1.0f / (1.0f+den) ;

		// saturate diffusion coefficent
		if (c < 0){C_cuda[index] = 0;}
		else if (c > 1) {C_cuda[index] = 1;}
		else {C_cuda[index] = c;}

		//barrier(CLK_LOCAL_MEM_FENCE);

		E_C[index] = e;
		W_C[index] = w;
		S_C[index] = s;
		N_C[index] = n;

	}).wait_and_throw();
}

static void srad_cuda_2(sycl::queue &q, sycl::nd_range<2> ndr, float *E_C, float *W_C, float *N_C, float *S_C,
		float *J_cuda, float *C_cuda, int cols, int rows, float lambda, float q0sqr)
{
	q.parallel_for<class srad_k2>(ndr, [=](sycl::nd_item<2> it) {
		//block id
		int bx = (int)it.get_group(1);
		int by = (int)it.get_group(0);

		//thread id
		int tx = (int)it.get_local_id(1);
		int ty = (int)it.get_local_id(0);


		//indices
		int index   = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + tx;
		int index_s = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * BLOCK_SIZE + tx;
		int index_e = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + BLOCK_SIZE;
		float cc, cn, cs, ce, cw, d_sum;

		//The below correspond to the index that was in localmem [ty+1][tx] and [ty][tx+1], so we load from global memory for the corresponding conditional branch
		int index_yp1_x = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * (ty+1) + tx;
		int index_y_xp1 = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + (tx+1);

		cc = C_cuda[index];
		cn = cc;
		cw = cc;


	//This kind of branching degrades fpga performance.
		if ( ty == BLOCK_SIZE -1 && tx == BLOCK_SIZE - 1){ //se
			if ( by == (int)it.get_group_range(0) - 1 )
				cs = C_cuda[cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + tx];
			else
				cs = C_cuda[index_s];

			if ( bx == (int)it.get_group_range(1) - 1 )
				ce = C_cuda[cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + tx];
			else
				ce = C_cuda[index_e];
		} 
		else if ( tx == BLOCK_SIZE -1 ){ //e
			cs  = C_cuda[index_yp1_x];

			if ( bx == (int)it.get_group_range(1) - 1 )
				ce = C_cuda[cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + tx];
			else
				ce = C_cuda[index_e];
		}
		else if ( ty == BLOCK_SIZE -1){ //s
			if ( by == (int)it.get_group_range(0) - 1 )
				cs = C_cuda[cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + tx];
			else
				cs = C_cuda[index_s];
			ce  = C_cuda[index_y_xp1];
		}
		else{ //the data elements which are not on the borders 
			cs  = C_cuda[index_yp1_x];
			ce  = C_cuda[index_y_xp1];
		}

		// divergence (equ 58)
		d_sum = cn * N_C[index] + cs * S_C[index] + cw * W_C[index] + ce * E_C[index];

		// image update (equ 61)
		J_cuda[index] = J_cuda[index] + 0.25f * lambda * d_sum;

		//barrier(CLK_GLOBAL_MEM_FENCE);

	}).wait_and_throw();
}

void random_matrix(float *I, int rows, int cols);
void runTest( int argc, char** argv);
void usage(int argc, char **argv)
{
	fprintf(stderr, "Usage: %s <rows> <cols> <y1> <y2> <x1> <x2> <lamda> <no. of iter> [platform & device]\n", argv[0]);
	fprintf(stderr, "\t<rows>   - number of rows\n");
	fprintf(stderr, "\t<cols>    - number of cols\n");
	fprintf(stderr, "\t<y1> 	 - y1 value of the speckle\n");
	fprintf(stderr, "\t<y2>      - y2 value of the speckle\n");
	fprintf(stderr, "\t<x1>       - x1 value of the speckle\n");
	fprintf(stderr, "\t<x2>       - x2 value of the speckle\n");
	fprintf(stderr, "\t<lamda>   - lambda (0,1)\n");
	fprintf(stderr, "\t<no. of iter>   - number of iterations\n");
	exit(1);
}


	int
main( int argc, char** argv) 
{
	dsycl::parse_args(&argc, &argv);
	runTest( argc, argv);
	return dwarf_report(srad_valid) == 0 ? EXIT_SUCCESS : EXIT_FAILURE;
}


	void
runTest( int argc, char** argv) 
{
	int rows, cols, size_I, size_R, niter = 10, iter;
	float *I, *J, lambda, q0sqr, sum, sum2, tmp, meanROI,varROI ;
	int i, j, k;

#ifdef CPU
	float Jc, G2, L, num, den, qsqr;
	int *iN,*iS,*jE,*jW;// k;
	float *dN,*dS,*dW,*dE;
	float cN,cS,cW,cE,D;
#endif

	sycl::queue q = dsycl::make_queue();
	float *J_cuda, *C_cuda, *E_C, *W_C, *N_C, *S_C;

	unsigned int r1, r2, c1, c2;
	float *c; /* [FIX-OCL-09] was int *c: the diffusion coefficient array truncated to 0/1 in the CPU path */



	if (argc == 9)
	{
		rows = atoi(argv[1]);  //number of rows in the domain
		cols = atoi(argv[2]);  //number of cols in the domain
		if ((rows%16!=0) || (cols%16!=0)){
			fprintf(stderr, "rows and cols must be multiples of 16\n");
			exit(1);
		}
		r1   = atoi(argv[3]);  //y1 position of the speckle
		r2   = atoi(argv[4]);  //y2 position of the speckle
		c1   = atoi(argv[5]);  //x1 position of the speckle
		c2   = atoi(argv[6]);  //x2 position of the speckle
		lambda = atof(argv[7]); //Lambda value
		niter = atoi(argv[8]); //number of iterations
	}
	else{
		usage(argc, argv);
	}



	size_I = cols * rows;
	size_R = (r2-r1+1)*(c2-c1+1);   

	I = (float*) memalign(AOCL_ALIGNMENT, sizeof(float)*size_I );
	J = (float*) memalign(AOCL_ALIGNMENT, sizeof(float)*size_I );
	c  = (float*) memalign(AOCL_ALIGNMENT, sizeof(float)*size_I) ;

      /*  I = (float*) malloc(sizeof(float)*size_I );
	J = (float*) malloc(sizeof(float)*size_I );
	c  = (float*) malloc(sizeof(float)*size_I) ;
      */
     


#ifdef CPU

	iN = (int *)malloc(sizeof(unsigned int*) * rows) ;
	iS = (int *)malloc(sizeof(unsigned int*) * rows) ;
	jW = (int *)malloc(sizeof(unsigned int*) * cols) ;
	jE = (int *)malloc(sizeof(unsigned int*) * cols) ;    


	dN = (float *)malloc(sizeof(float)* size_I) ;
	dS = (float *)malloc(sizeof(float)* size_I) ;
	dW = (float *)malloc(sizeof(float)* size_I) ;
	dE = (float *)malloc(sizeof(float)* size_I) ;    


	for (i=0; i< rows; i++) {
		iN[i] = i-1;
		iS[i] = i+1;
	}    
	for (j=0; j< cols; j++) {
		jW[j] = j-1;
		jE[j] = j+1;
	}
	iN[0]    = 0;
	iS[rows-1] = rows-1;
	jW[0]    = 0;
	jE[cols-1] = cols-1;

#endif

	//Allocate device memory
	J_cuda = dsycl::alloc<float>(q, size_I);
	C_cuda = dsycl::alloc<float>(q, size_I);
	E_C = dsycl::alloc<float>(q, size_I);
	W_C = dsycl::alloc<float>(q, size_I);
	S_C = dsycl::alloc<float>(q, size_I);
	N_C = dsycl::alloc<float>(q, size_I);
	dsycl::prebuild(q);   /* kernel JIT, as clBuildProgram: outside the timed region */

	printf("Randomizing the input matrix\n");
	//Generate a random matrix
	random_matrix(I, rows, cols);

	for (k = 0;  k < size_I; k++ ) {
		J[k] = (float)exp((double)I[k]) ; /* double exp, as in the C (OpenCL tree) build */
	}
	printf("Start the SRAD main loop\n");
	/* common timer: repetitions around the whole iteration loop (host ROI statistics + J transfers + both kernels) */
	float *J0 = (float*) malloc(sizeof(float)*size_I);
	memcpy(J0, J, sizeof(float)*size_I);
	dwarf_init("srad", "structured-grids", "sycl");
	dwarf_set_size(rows, cols);
	int rep;
	for(rep = dwarf_first_rep(); rep < dwarf_reps(); rep++)
	{
	memcpy(J, J0, sizeof(float)*size_I);
	dwarf_region_begin();
	for (iter=0; iter< niter; iter++){     
		sum=0; sum2=0;
		for (i=r1; i<=r2; i++) {
			for (j=c1; j<=c2; j++) {
				tmp   = J[i * cols + j];
				sum  += tmp ;
				sum2 += tmp*tmp;
			}
		}
		meanROI = sum / size_R;
		varROI  = (sum2 / size_R) - meanROI*meanROI;
		q0sqr   = varROI / (meanROI*meanROI);

#ifdef CPU

		for (i = 0 ; i < rows ; i++) {
			for (j = 0; j < cols; j++) { 

				k = i * cols + j;
				Jc = J[k];

				// directional derivates
				dN[k] = J[iN[i] * cols + j] - Jc;
				dS[k] = J[iS[i] * cols + j] - Jc;
				dW[k] = J[i * cols + jW[j]] - Jc;
				dE[k] = J[i * cols + jE[j]] - Jc;

				G2 = (dN[k]*dN[k] + dS[k]*dS[k] 
						+ dW[k]*dW[k] + dE[k]*dE[k]) / (Jc*Jc);

				L = (dN[k] + dS[k] + dW[k] + dE[k]) / Jc;

				num  = (0.5*G2) - ((1.0/16.0)*(L*L)) ;
				den  = 1 + (.25*L);
				qsqr = num/(den*den);

				// diffusion coefficent (equ 33)
				den = (qsqr-q0sqr) / (q0sqr * (1+q0sqr)) ;
				c[k] = 1.0 / (1.0+den) ;

				// saturate diffusion coefficent
				if (c[k] < 0) {c[k] = 0;}
				else if (c[k] > 1) {c[k] = 1;}
			}
		}
		for (i = 0; i < rows; i++) {
			for (j = 0; j < cols; j++) {        

				// current index
				k = i * cols + j;

				// diffusion coefficent
				cN = c[k];
				cS = c[iS[i] * cols + j];
				cW = c[k];
				cE = c[i * cols + jE[j]];

				// divergence (equ 58)
				D = cN * dN[k] + cS * dS[k] + cW * dW[k] + cE * dE[k];

				// image update (equ 61)
				J[k] = J[k] + 0.25*lambda*D;
			}
		}

#endif // CPU


		//Currently the input size must be divided by 16 - the block size
		int block_x = cols/BLOCK_SIZE ;
		int block_y = rows/BLOCK_SIZE ;

		/* [FIX-OCL-34] algorithmic work-group: 16x16 tile, kept; rows/cols are multiples of 16, so no padding/guard is needed. */
		sycl::nd_range<2> ndr({(size_t)block_y*BLOCK_SIZE, (size_t)block_x*BLOCK_SIZE}, {BLOCK_SIZE, BLOCK_SIZE});

		//Copy data from main memory to device memory
		dsycl::write(q, J_cuda, J, sizeof(float)*size_I);
		//Run kernels
		srad_cuda_1(q, ndr, E_C, W_C, N_C, S_C, J_cuda, C_cuda, cols, rows, q0sqr);
		srad_cuda_2(q, ndr, E_C, W_C, N_C, S_C, J_cuda, C_cuda, cols, rows, lambda, q0sqr);
		//Copy data from device memory to main memory
		dsycl::read(q, J, J_cuda, sizeof(float)*size_I);
	}

	q.wait_and_throw();
	dwarf_region_end(rep);
	}

#ifdef GPU
	{
		float *Jr = (float*) malloc(sizeof(float)*size_I); long bad = 0;
		memcpy(Jr, J0, sizeof(float)*size_I);
		srad_reference(Jr, rows, cols, r1, r2, c1, c2, lambda, niter);
		for (k = 0; k < size_I; k++) { float d = fabsf(Jr[k]-J[k]), m = fmaxf(1.0f, fmaxf(fabsf(Jr[k]), fabsf(J[k]))); if (!(d/m <= 1e-3f)) bad++; }
		printf("%s: %ld pixels differ from the serial host reference (rel. tol 1e-3)\n", bad ? "FAILURE" : "SUCCESS", bad);
		srad_valid = (bad == 0);
		free(Jr);
	}
#endif

#ifdef OUTPUT
	//Printing output	
	printf("Printing Output:\n"); 
	for( i = 0 ; i < rows ; i++){
		for ( j = 0 ; j < cols ; j++){
			printf("%.5f ", J[i * cols + j]); 
		}	
		printf("\n"); 
	}
#endif 

	printf("Computation Done\n");

	free(I);
	free(J);
#ifdef CPU
	free(iN); free(iS); free(jW); free(jE);
	free(dN); free(dS); free(dW); free(dE);
#endif
	dsycl::release(q, C_cuda); dsycl::release(q, J_cuda); dsycl::release(q, E_C);
	dsycl::release(q, W_C); dsycl::release(q, N_C); dsycl::release(q, S_C);
	free(c);

}


void random_matrix(float *I, int rows, int cols){

	int i, j;
	srand(7);

	for( i = 0 ; i < rows ; i++){
		for ( j = 0 ; j < cols ; j++){
			I[i * cols + j] = rand()/(float)RAND_MAX ;
		}
	}

}
