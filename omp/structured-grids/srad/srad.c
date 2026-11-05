/* SRAD (Speckle Reducing Anisotropic Diffusion) -- OpenMP target-offload translation of
 * opencl/structured-grids/srad/srad.c + srad_kernel.cl (Rodinia SRAD, OpenDwarfs variant).
 *
 * Same host flow: for every iteration { ROI mean/variance on the host -> q0sqr; copy J to the
 * device; srad_cuda_1 (derivatives + diffusion coefficient); srad_cuda_2 (divergence + image
 * update); copy J back }. Same kernels (index arithmetic and border handling copied from
 * srad_kernel.cl), same random input (srand(7)), same CLI:
 *     srad <rows> <cols> <y1> <y2> <x1> <x2> <lambda> <niter>      (rows, cols multiples of 16)
 * The two kernels have no barriers, so the BLOCK_SIZE x BLOCK_SIZE work-groups are expressed as a
 * flat 4-D parallel loop over (group_y, group_x, local_y, local_x).
 * Timed region: the whole iteration loop (including the per-iteration host ROI reduction, which is
 * part of the algorithm), device allocation, per-iteration J transfers.
 * Validation (added): serial host evaluation of the same float kernels, compared with a relative
 * tolerance of 1e-3 (GPU and host may contract multiply-adds differently).
 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "srad.h"
#include "dwarf_timer.h"

#define BS 16  /* BLOCK_SIZE of the kernels */

void random_matrix(float *I, int rows, int cols);

static void usage(int argc, char **argv) {
  fprintf(stderr, "Usage: %s <rows> <cols> <y1> <y2> <x1> <x2> <lamda> <no. of iter>\n", argv[0]);
  exit(1);
}

#pragma omp declare target
/* ---- srad_cuda_1 for one work-item ---- */
static inline void srad1_item(float *E_C, float *W_C, float *N_C, float *S_C, const float *J, float *C, int cols,
                              int rows, float q0sqr, int bx, int by, int tx, int ty, int nbx, int nby) {
  const int BLOCK_SIZE = BS;


	//block id
	
	

	//thread id
	
	

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
	else if ( by == nby - 1 ){
		index_s = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ( BLOCK_SIZE - 1 ) + tx;
	}

	if ( bx == 0 ){
		index_w = cols * BLOCK_SIZE * by + cols * ty;
	}
	else if ( bx == nbx - 1 ){
		index_e = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + BLOCK_SIZE-1;
	}

	//Corresponding to [ty+1][x], [ty-1][x], etc.
	int index_yp1_x = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * (ty + 1) + tx;
	int index_ym1_x = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * (ty - 1) + tx;
	int index_y_xp1 = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols  * ty + (tx + 1);
	int index_y_xm1 = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + (tx - 1);


	jc = J[index];

	if ( ty == 0 && tx == 0 ){ //nw
		n  = J[index_n] - jc;
		s  = J[index_yp1_x] - jc;
		w  = J[index_w]  - jc; 
		e  = J[index_y_xp1] - jc;
	}	    
	else if ( ty == 0 && tx == BLOCK_SIZE-1 ){ //ne
		n  = J[index_n] - jc;
		s  = J[index_yp1_x] - jc;
		w  = J[index_y_xm1] - jc;
		e  = J[index_e] - jc;
	}
	else if ( ty == BLOCK_SIZE -1 && tx == BLOCK_SIZE - 1){ //se
		n  = J[index_ym1_x] - jc;
		s  = J[index_s] - jc;
		w  = J[index_y_xm1] - jc; 
		e  = J[index_e]  - jc;
	}
	else if ( ty == BLOCK_SIZE -1 && tx == 0 ){//sw
		n  = J[index_ym1_x] - jc;
		s  = J[index_s] - jc;
		w  = J[index_w]  - jc; 
		e  = J[index_y_xp1] - jc;
	}
	else if ( ty == 0 ){ //n
		n  = J[index_n] - jc;
		s  = J[index_yp1_x] - jc;
		w  = J[index_y_xm1] - jc; 
		e  = J[index_y_xp1] - jc;
	}
	else if ( tx == BLOCK_SIZE -1 ){ //e
		n  = J[index_ym1_x] - jc;
		s  = J[index_yp1_x] - jc;
		w  = J[index_y_xm1] - jc; 
		e  = J[index_e] - jc;
	}
	else if ( ty == BLOCK_SIZE -1){ //s
		n  = J[index_ym1_x] - jc;
		s  = J[index_s] - jc;
		w  = J[index_y_xm1] - jc; 
		e  = J[index_y_xp1] - jc;
	}
	else if ( tx == 0 ){ //w
		n  = J[index_ym1_x] - jc;
		s  = J[index_yp1_x] - jc;
		w  = J[index_w] - jc; 
		e  = J[index_y_xp1] - jc;
	}
	else{  //the data elements which are not on the borders 
		n  = J[index_ym1_x] - jc;
		s  = J[index_yp1_x] - jc;
		w  = J[index_y_xm1] - jc; 
		e  = J[index_y_xp1] - jc;
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
	if (c < 0){C[index] = 0;}
	else if (c > 1) {C[index] = 1;}
	else {C[index] = c;}

	//barrier(CLK_LOCAL_MEM_FENCE);

	E_C[index] = e;
	W_C[index] = w;
	S_C[index] = s;
	N_C[index] = n;

}

/* ---- srad_cuda_2 for one work-item ---- */
static inline void srad2_item(const float *E_C, const float *W_C, const float *N_C, const float *S_C, float *J,
                              const float *C, int cols, int rows, float lambda, float q0sqr, int bx, int by, int tx,
                              int ty, int nbx, int nby) {
  const int BLOCK_SIZE = BS;

	//block id
	
	

	//thread id
	
	


	//indices
	int index   = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + tx;
	int index_s = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * BLOCK_SIZE + tx;
	int index_e = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + BLOCK_SIZE;
	float cc, cn, cs, ce, cw, d_sum;

	//The below correspond to the index that was in localmem [ty+1][tx] and [ty][tx+1], so we load from global memory for the corresponding conditional branch
	int index_yp1_x = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * (ty+1) + tx;
	int index_y_xp1 = cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + (tx+1);

	cc = C[index];
	cn = cc;
	cw = cc;


//This kind of branching degrades fpga performance.
	if ( ty == BLOCK_SIZE -1 && tx == BLOCK_SIZE - 1){ //se
		if ( by == nby - 1 )
			cs = C[cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + tx];
		else
			cs = C[index_s];

		if ( bx == nbx - 1 )
			ce = C[cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + tx];
		else
			ce = C[index_e];
	} 
	else if ( tx == BLOCK_SIZE -1 ){ //e
		cs  = C[index_yp1_x];

		if ( bx == nbx - 1 )
			ce = C[cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + tx];
		else
			ce = C[index_e];
	}
	else if ( ty == BLOCK_SIZE -1){ //s
		if ( by == nby - 1 )
			cs = C[cols * BLOCK_SIZE * by + BLOCK_SIZE * bx + cols * ty + tx];
		else
			cs = C[index_s];
		ce  = C[index_y_xp1];
	}
	else{ //the data elements which are not on the borders 
		cs  = C[index_yp1_x];
		ce  = C[index_y_xp1];
	}

	// divergence (equ 58)
	d_sum = cn * N_C[index] + cs * S_C[index] + cw * W_C[index] + ce * E_C[index];

	// image update (equ 61)
	J[index] = J[index] + 0.25f * lambda * d_sum;

	//barrier(CLK_GLOBAL_MEM_FENCE);

}
#pragma omp end declare target

static void run_kernels(float *E_C, float *W_C, float *N_C, float *S_C, float *J, float *C, int cols, int rows,
                        float lambda, float q0sqr, size_t n) {
  int nbx = cols / BS, nby = rows / BS;
#pragma omp target teams distribute parallel for collapse(4) map(present, alloc : E_C[0 : n], W_C[0 : n], N_C[0 : n], S_C[0 : n], J[0 : n], C[0 : n])
  for (int by = 0; by < nby; by++)
    for (int bx = 0; bx < nbx; bx++)
      for (int ty = 0; ty < BS; ty++)
        for (int tx = 0; tx < BS; tx++) srad1_item(E_C, W_C, N_C, S_C, J, C, cols, rows, q0sqr, bx, by, tx, ty, nbx, nby);
#pragma omp target teams distribute parallel for collapse(4) map(present, alloc : E_C[0 : n], W_C[0 : n], N_C[0 : n], S_C[0 : n], J[0 : n], C[0 : n])
  for (int by = 0; by < nby; by++)
    for (int bx = 0; bx < nbx; bx++)
      for (int ty = 0; ty < BS; ty++)
        for (int tx = 0; tx < BS; tx++) srad2_item(E_C, W_C, N_C, S_C, J, C, cols, rows, lambda, q0sqr, bx, by, tx, ty, nbx, nby);
}

static void roi_q0(const float *J, int cols, unsigned r1, unsigned r2, unsigned c1, unsigned c2, int size_R, float *q0sqr) {
  float sum = 0, sum2 = 0, tmp;
  for (unsigned i = r1; i <= r2; i++)
    for (unsigned j = c1; j <= c2; j++) { tmp = J[i * cols + j]; sum += tmp; sum2 += tmp * tmp; }
  float meanROI = sum / size_R;
  float varROI = (sum2 / size_R) - meanROI * meanROI;
  *q0sqr = varROI / (meanROI * meanROI);
}

int main(int argc, char **argv) {
  int rows, cols, size_I, size_R, niter = 10, iter, k;
  float *I, *J, *J0, lambda, q0sqr;
  unsigned int r1, r2, c1, c2;

  if (argc == 9) {
    rows = atoi(argv[1]);
    cols = atoi(argv[2]);
    if ((rows % 16 != 0) || (cols % 16 != 0)) { fprintf(stderr, "rows and cols must be multiples of 16\n"); exit(1); }
    r1 = atoi(argv[3]); r2 = atoi(argv[4]); c1 = atoi(argv[5]); c2 = atoi(argv[6]);
    lambda = atof(argv[7]);
    niter = atoi(argv[8]);
  } else usage(argc, argv);
  size_I = cols * rows;
  size_R = (r2 - r1 + 1) * (c2 - c1 + 1);
  size_t n = (size_t)size_I;
  I = (float *)malloc(sizeof(float) * size_I);
  J = (float *)malloc(sizeof(float) * size_I);
  J0 = (float *)malloc(sizeof(float) * size_I);
  float *E_C = (float *)malloc(sizeof(float) * n), *W_C = (float *)malloc(sizeof(float) * n),
        *N_C = (float *)malloc(sizeof(float) * n), *S_C = (float *)malloc(sizeof(float) * n),
        *C = (float *)malloc(sizeof(float) * n);

  printf("Randomizing the input matrix\n");
  random_matrix(I, rows, cols);
  for (k = 0; k < size_I; k++) J0[k] = (float)exp(I[k]);

  printf("Start the SRAD main loop\n");
  dwarf_init("srad", "structured-grids", "omp");
  dwarf_set_size(rows, cols);
  for (int rep = dwarf_first_rep(); rep < dwarf_reps(); ++rep) {
    memcpy(J, J0, sizeof(float) * size_I); /* untimed */
    dwarf_region_begin();
#pragma omp target enter data map(alloc : J[0 : n], C[0 : n], E_C[0 : n], W_C[0 : n], N_C[0 : n], S_C[0 : n])
    for (iter = 0; iter < niter; iter++) {
      roi_q0(J, cols, r1, r2, c1, c2, size_R, &q0sqr);
#pragma omp target update to(J[0 : n])
      run_kernels(E_C, W_C, N_C, S_C, J, C, cols, rows, lambda, q0sqr, n);
#pragma omp target update from(J[0 : n])
    }
#pragma omp target exit data map(delete : J[0 : n], C[0 : n], E_C[0 : n], W_C[0 : n], N_C[0 : n], S_C[0 : n])
    dwarf_region_end(rep);
  }
  printf("Computation Done\n");

  /* validation: serial evaluation of the same kernels on the host */
  float *Jr = (float *)malloc(sizeof(float) * n);
  memcpy(Jr, J0, sizeof(float) * n);
  {
    int nbx = cols / BS, nby = rows / BS;
    for (iter = 0; iter < niter; iter++) {
      roi_q0(Jr, cols, r1, r2, c1, c2, size_R, &q0sqr);
      for (int by = 0; by < nby; by++) for (int bx = 0; bx < nbx; bx++) for (int ty = 0; ty < BS; ty++) for (int tx = 0; tx < BS; tx++)
        srad1_item(E_C, W_C, N_C, S_C, Jr, C, cols, rows, q0sqr, bx, by, tx, ty, nbx, nby);
      for (int by = 0; by < nby; by++) for (int bx = 0; bx < nbx; bx++) for (int ty = 0; ty < BS; ty++) for (int tx = 0; tx < BS; tx++)
        srad2_item(E_C, W_C, N_C, S_C, Jr, C, cols, rows, lambda, q0sqr, bx, by, tx, ty, nbx, nby);
    }
  }
  long bad = 0;
  for (k = 0; k < size_I; k++) {
    float d = fabsf(Jr[k] - J[k]), m = fmaxf(1.0f, fmaxf(fabsf(Jr[k]), fabsf(J[k])));
    if (!(d / m <= 1e-3f)) bad++;
  }
  printf("%s: %ld pixels differ from the serial host evaluation (rel. tol 1e-3)\n", bad ? "FAILURE" : "SUCCESS", bad);
  free(I); free(J); free(J0); free(Jr); free(E_C); free(W_C); free(N_C); free(S_C); free(C);
  return dwarf_report(bad == 0);
}

void random_matrix(float *I, int rows, int cols) {
  int i, j;
  srand(7);
  for (i = 0; i < rows; i++)
    for (j = 0; j < cols; j++) I[i * cols + j] = rand() / (float)RAND_MAX;
}
