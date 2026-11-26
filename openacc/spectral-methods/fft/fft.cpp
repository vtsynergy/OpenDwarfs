/* Spectral-methods dwarf: FFT -- OpenACC translation of
 * opencl/spectral-methods/fft (SHOC FFT, Cooley-Tukey radix-8/16 kernels with local-memory transposes).
 *
 * Covered: the default one-dimensional path (--pts N, N = 128 ... 8192 points, one transform):
 *   N <= 2048 : one in-place kernel fft1D_N, one work-group;
 *   N = 4096, 8192 : two kernels fft0 (global transpose stage) -> fft1D_N.
 * OpenACC has no barrier inside a gang: each kernel is split at its barriers into vector loops (see the
 * comment above the kernels). Not translated: the two-dimensional path (--2D).
 * [FIX-OCL-24] pow1 = log2(fftn); [FIX-OCL-25] only row jj==0 stores (128/256 points);
 * [FIX-OCL-26] srand(1).  Timed region (as OpenCL): device allocation + H2D + forward FFT + D2H.
 * Validation as OpenCL ([FIX-OCL-27]): double-precision DFT of the same input.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <openacc.h>
#include "dwarf_acc.h"
#include "dwarf_timer.h"

typedef float T;
struct T2 { float x, y; };
#define FFT_128
#define FFT_256
#define FFT_512
#define FFT_1024
#define FFT_2048
#define FFT_4096
#define FFT_8192

#pragma acc routine seq
static inline T2 mk2(float x, float y) { T2 r; r.x = x; r.y = y; return r; }
#pragma acc routine seq
static inline int mad24(int a, int b, int c) { return a * b + c; }
#pragma acc routine seq
static inline int mul24(int a, int b) { return a * b; }
#pragma acc routine seq
static inline int ocl_shr(int x, int s) { return (s >= 32 || s < 0) ? 0 : (x >> s); }
#pragma acc routine seq
static inline int ocl_shl(int x, int s) { return (s >= 32 || s < 0) ? 0 : (x << s); }

#define numIter Radix1/Radix2
#define STORE_REV 1
#define LOAD_REV 2
#define NO_REV 0
#define STORE_S 3
#define LOAD_S 4


#define COS_PI_8  0.923879533f
#define SIN_PI_8  0.382683432f
#define exp_1_8   mk2(  M_SQRT1_2, M_SQRT1_2 )//requires post-multiply by 1/sqrt(2)
#define exp_1_4   mk2(  0, 1 )
#define exp_3_8   mk2( -M_SQRT1_2, M_SQRT1_2 )//requires post-multiply by 1/sqrt(2)
#define exp_1_16  mk2(  COS_PI_8, SIN_PI_8 )
#define exp_3_16  mk2(  M_SQRT1_2, M_SQRT1_2 )
#define exp_5_16  mk2( SIN_PI_8, COS_PI_8 )
#define exp_7_16  mk2( -M_SQRT1_2, M_SQRT1_2 )
#define exp_9_16  mk2( -COS_PI_8,  -SIN_PI_8 )
#pragma acc routine seq
static inline T2 exp_i( T phi ) {
	return mk2( cosf(phi), sinf(phi) );
}
#pragma acc routine seq
static inline T2 cmplx_mul( T2 a, T2 b ) { return mk2( a.x*b.x-a.y*b.y, a.x*b.y+a.y*b.x ); }
#pragma acc routine seq
static inline T2 cm_fl_mul( T2 a, T  b ) { return mk2( b*a.x, b*a.y ); }
#pragma acc routine seq
static inline T2 cmplx_add( T2 a, T2 b ) { return mk2( a.x + b.x, a.y + b.y ); }
#pragma acc routine seq
static inline T2 cmplx_sub( T2 a, T2 b ) { return mk2( a.x - b.x, a.y - b.y ); }
#pragma acc routine seq
static inline T2 cmplx_swap( T2 a) {return mk2( a.y, a.x);}

#define FFT2(a0, a1)                            \
{                                               \
	T2 c0 = *a0;                           \
	*a0 = cmplx_add(c0,*a1);                    \
	*a1 = cmplx_sub(c0,*a1);                    \
}
#define FFT2S(a)                            \
{                                               \
	T2 c0 = *(a);                           \
	*(a) = cmplx_add(c0,*(a+1));                    \
	*(a+1) = cmplx_sub(c0,*(a+1));                    \
}
#define FFT4(a0, a1, a2, a3)                    \
{                                               \
	FFT2( (a0), (a2) );                             \
	FFT2( (a1), (a3) );                             \
	*(a3)  = cmplx_mul(*(a3),  exp_1_4 );		\
	FFT2( (a0), (a1) );                             \
	FFT2( (a2), (a3) );                             \
	T2 c = *(a1); \
	*(a1) = *(a2); \
	*(a2) = c; \
}

#define FFT4S(a)                    \
{                                               \
	FFT2( (a+0), (a+2) );                             \
	FFT2( (a+1), (a+3) );                             \
	*(a+3)  = cmplx_mul(*(a+3),  exp_1_4 );			\
	FFT2( (a+0), (a+1) );                             \
	FFT2( (a+2), (a+3) );                             \
	T2 c = *(a+1); \
	*(a+1) = *(a+2); \
	*(a+2) = c; \
}

#define FFT8(a)                                                 \
{                                                               \
	FFT2( &(a)[0], &(a)[4] );                                       \
	FFT2( &(a)[1], &(a)[5] );                                       \
	FFT2( &(a)[2], &(a)[6] );                                       \
	FFT2( &(a)[3], &(a)[7] );                                       \
	\
	(a)[5] = cmplx_mul( (a)[5],exp_1_8);    \
	(a)[6]  = cmplx_mul((a)[6],  exp_1_4 );			\
	(a)[7] = cmplx_mul( (a)[7],exp_3_8);    \
	FFT2(&(a)[0], &(a)[2]); \
	FFT2(&(a)[1], &(a)[3]); \
	FFT2(&(a)[4], &(a)[6]); \
	FFT2(&(a)[5], &(a)[7]); \
	(a)[3]  = cmplx_mul((a)[3],  exp_1_4 );			\
	(a)[7]  = cmplx_mul((a)[7],  exp_1_4 );			\
	FFT2(&(a)[0], &(a)[1]); \
	FFT2(&(a)[2], &(a)[3]); \
	FFT2(&(a)[4], &(a)[5]); \
	FFT2(&(a)[6], &(a)[7]); \
	T2 c; \
	c = (a)[1]; \
	(a)[1] = (a)[4]; \
	(a)[4] = c; \
	c = (a)[3]; \
	(a)[3] = (a)[6]; \
	(a)[6] = c; \
}



#define FFT16(a)	\
{\
	FFT4( &(a)[0], &(a)[4], &(a)[8], &(a)[12] );			\
	FFT4( &(a)[1], &(a)[5], &(a)[9], &(a)[13] );			\
	FFT4( &(a)[2], &(a)[6], &(a)[10], &(a)[14] );			\
	FFT4( &(a)[3], &(a)[7], &(a)[11], &(a)[15] );			\
	\
	(a)[5]  = cmplx_mul((a)[5], exp_1_16); \
	(a)[6]  = cmplx_mul((a)[6], exp_3_16); \
	(a)[7]  = cmplx_mul((a)[7], exp_5_16); \
	(a)[9]  = cmplx_mul((a)[9], exp_3_16); \
	(a)[10]  = cmplx_mul((a)[10],  exp_1_4 );                 \
	(a)[11] = cmplx_mul((a)[11], exp_7_16); \
	(a)[13] = cmplx_mul((a)[13], exp_5_16); \
	(a)[14] = cmplx_mul((a)[14], exp_7_16); \
	(a)[15] = cmplx_mul((a)[15], exp_9_16); \
	\
	FFT4( &(a)[0],  &(a)[1],  &(a)[2],  &(a)[3] );			\
	FFT4( &(a)[4],  &(a)[5],  &(a)[6],  &(a)[7] );			\
	FFT4( &(a)[8],  &(a)[9],  &(a)[10], &(a)[11] );			\
	FFT4( &(a)[12], &(a)[13], &(a)[14], &(a)[15] );			\
	T2 c; \
	c = (a)[1];  (a)[1]  = (a)[4];  (a)[4]  = c; \
	c = (a)[2];  (a)[2]  = (a)[8];  (a)[8]  = c; \
	c = (a)[3];  (a)[3]  = (a)[12]; (a)[12] = c; \
	c = (a)[6];  (a)[6]  = (a)[9];  (a)[9]  = c; \
	c = (a)[7];  (a)[7]  = (a)[13]; (a)[13] = c; \
	c = (a)[11]; (a)[11] = (a)[14]; (a)[14] = c; \
}




#pragma acc routine seq
static inline void fftKernel(T2 *a, int n)
{
	switch(n)
	{
		case 2:
			FFT2S(a);
			break;
		case 4:
			FFT4S(a);
			break;
		case 8:
			FFT8(a);
			break;
		case 16:
			FFT16(a);
			break;
		default:
			return;
	}
}



#pragma acc routine seq
static inline void globalLoads8(T2 *data, T2 *in, int S){
	for( int i = 0; i < 8; i++ )
		data[i] = in[i*S];
}

#pragma acc routine seq
static inline void globalLoads16(T2 *data, T2 *in, int S){
	for( int i = 0; i < 16; i++ )
		data[i] = in[i*S];
}


#pragma acc routine seq
static inline void globalStores8_r(T2 *data, T2 *out, int S){
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
#ifdef FFT_2D
		out[i*S] = data[reversed[i]];
#else
	out[i*S] = cmplx_swap(data[reversed[i]]);
#endif
}


#pragma acc routine seq
static inline void globalStores8(T2 *data, T2 *out, int S){
	for( int i = 0; i < 8; i++ )
#ifdef FFT_2D
		out[i*S] = data[i];
#else
	out[i*S] = cmplx_swap(data[i]);
#endif
}
#pragma acc routine seq
static inline void globalStores16(T2 *data, T2 *out, int S){
	for( int i = 0; i < 16; i++ )
		out[i*S] = data[i];

}

#pragma acc routine seq
static inline void storex4( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 4; i++ )
		x[i*sx] = a[i].x;
}

#pragma acc routine seq
static inline void storey4( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 4; i++ )
		x[i*sx] = a[i].y;
}

#pragma acc routine seq
static inline void storex16( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 16; i++ )
		x[i*sx] = a[i].x;
}

#pragma acc routine seq
static inline void storey16( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 16; i++ )
		x[i*sx] = a[i].y;
}


#pragma acc routine seq
static inline void storex8( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 8; i++ )
		x[i*sx] = a[i].x;
}

#pragma acc routine seq
static inline void storey8( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 8; i++ )
		x[i*sx] = a[i].y;
}

#pragma acc routine seq
static inline void storex8_r( T2 *a, T *x, int sx ) {
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
		x[i*sx] = a[reversed[i]].x;
}

#pragma acc routine seq
static inline void storey8_r( T2 *a, T *x, int sx ) {
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
		x[i*sx] = a[reversed[i]].y;
}

#pragma acc routine seq
static inline void loadx16( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 16; i++ )
		a[i].x = x[i*sx];
}

#pragma acc routine seq
static inline void loady16( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 16; i++ )
		a[i].y = x[i*sx];
}

#pragma acc routine seq
static inline void loadx8( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 8; i++ )
		a[i].x = x[i*sx];
}

#pragma acc routine seq
static inline void loady8( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 8; i++ )
		a[i].y = x[i*sx];
}

#pragma acc routine seq
static inline void loadx4( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 4; i++ )
		a[i].x = x[i*sx];
}

#pragma acc routine seq
static inline void loady4( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 4; i++ )
		a[i].y = x[i*sx];
}

#pragma acc routine seq
static inline void loadx8_r( T2 *a, T *x, int sx ) {
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
		a[i].x = x[reversed[i]*sx];
}

#pragma acc routine seq
static inline void loady8_r( T2 *a, T *x, int sx ) {
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
		a[i].y = x[reversed[i]*sx];
}

#pragma acc routine seq
static inline void globalLoads(T2 *data, T2 *in, int S, int n){
	for( int i = 0; i < n; i++ )
		data[i] = in[i*S];
}

#pragma acc routine seq
static inline void globalStores(T2 *data, T2 *out, int S, int n, bool rev){
	if(rev){
		for( int i = 0; i < n/2; i++)
			for( int j=0; (j*n/2+i)<n; j++)
				out[(i*2+j)*S] = cmplx_swap(data[j*n/2+i]);
	}
	else{
		for( int i = 0; i < n; i++)
			out[i*S] = cmplx_swap(data[i]);
	}
}

#pragma acc routine seq
static inline void globalStores_s(T2 *data, T2 *out, int S, int n, int offset){
	for( int i = 0; i < n/2; i++)
#ifdef TWIDDLE
		out[i*S] = data[i];
#else
	out[i*S] = cmplx_swap(data[i]);
#endif
	for( int i = 0; i < n/2; i++)
#ifdef TWIDDLE
		out[i*S+offset]= data[i+n/2];
#else
	out[i*S+offset]= cmplx_swap(data[i+n/2]);
#endif
}

#pragma acc routine seq
static inline void storex( T2 *a, T *x, int sx, int n, bool rev, int offset){
	if(rev){
		for( int i = 0; i < n/2; i++)
			for(int j=0; (j*n/2+i)<n; j++)
				x[(i*2+j)*sx] = a[j*n/2+i].x;
	}
	else{
		for( int i = offset; i < n+offset; i++)
			x[i*sx] = a[i].x;
	}
}

#pragma acc routine seq
static inline void storey( T2 *a, T *x, int sx, int n, bool rev, int offset){
	if(rev){
		for( int i = 0; i < n/2; i++)
			for(int j=0; (j*n/2+i)<n; j++)
				x[(i*2+j)*sx] = a[j*n/2+i].y;
	}
	else{
		for( int i = offset; i < n+offset; i++)
			x[i*sx] = a[i].y;
	}
}

#pragma acc routine seq
static inline void loadx( T2 *a, T *x, int sx, int n, bool rev, int offset ) {
	if(rev){
		for( int i = 0; i < n/2; i++)
			for(int j=0; (j*n/2+i)<n; j++)
				a[(i*2+j)].x = x[(j*n/2+i)*sx];
	}
	else{
		for( int i = offset; i < n+offset; i++)
			a[i].x = x[(i)*sx];
	}
}
#pragma acc routine seq
static inline void loady( T2 *a, T *x, int sx, int n, bool rev,int offset ) {
	if(rev){
		for( int i = 0; i < n/2; i++)
			for(int j=0; (j*n/2+i)<n; j++)
				a[(i*2+j)].y = x[(j*n/2+i)*sx];
	}
	else{
		for( int i = offset; i < n+offset; i++)
			a[i].y = x[(i)*sx];
	}
}





#pragma acc routine seq
static inline void twiddle_s(T2 *a, int j, int N, int n, int k )                                              \
{                                                                       \
	for( int i = 0; i < n/2; i++ ){			\
		a[i] = cmplx_mul( a[i],exp_i(((2.0f*M_PI/N)*j)*(k+16*i))); \
	} 							\
	for( int i = n/2; i < n; i++ ){			\
		a[i] = cmplx_mul( a[i],exp_i(((2.0f*M_PI/N)*j)*(k+16*(i-n/2)+1))); \
	}                                                                   \
	\
}




#pragma acc routine seq
static inline void twiddle16(T2 *a, int tt, int n)
{
	T c = (T) tt;
	for( int t = 1; t < 16; t++ ){
		(a)[t] = cmplx_mul( (a)[t],exp_i((2.0f*M_PI*t/n)*c));
	}
}
#pragma acc routine seq
static inline void twiddle8(T2 *a, int tt, int n)
{
	T c = (T) tt;
	for( int t = 1; t < 8; t++ ){
		(a)[t] = cmplx_mul( (a)[t],exp_i((2.0f*M_PI*t/n)*c));
	}
}
#pragma acc routine seq
static inline void  twiddle4(T2 *a, int tt, int n )                                 {
	T c = (T) tt;
	for( int t = 1; t < 4; t++ ){
		(a)[t] = cmplx_mul( (a)[t],exp_i((2.0f*M_PI*t/n)*c));
	}
}


/* ======================= kernels of fft.cl =======================================================
 * One OpenCL work-group = one gang; its L work-items = the vector lanes. __local smem and the
 * per-work-item private registers a[] (here A[lane][]) are gang-private arrays. Every OpenCL
 * barrier(CLK_LOCAL_MEM_FENCE) is the end of one vector loop (PHASE) and the start of the next, so a
 * kernel is its fft.cl body split at the barriers; index expressions that are pure functions of
 * (work-item, work-group) are recomputed in each phase. local_id is the OpenCL get_local_id(0).
 * Each transpose() of fft.cl is four phases: store x | load x | store y | load y.
 * ============================================================================================== */
#define PHASE(...) _Pragma("acc loop vector independent") for (int local_id = 0; local_id < L; local_id++) { T2 *a = A[local_id]; __VA_ARGS__ }

static void launch_fft1D_512(T2 *in_all, size_t n, int G) {
  enum { L = 64 };
  T2 A[L][16]; T smem[576];
#pragma acc parallel loop gang ACC_SHAPE(G, 64) private(A, smem) present(in_all[0 : n])
  for (int group_id = 0; group_id < G; group_id++) {
#define IDX int ii = local_id; T2 *in = in_all + mad24(group_id, 512, ii);
#define S1 (smem + ii)
#define L1 (smem + mad24(ii & 7, 66, ii >> 3))
#define L2 (smem + mad24(ii >> 3, 72, ii & 7))
    PHASE(IDX globalLoads8(a, in, 64); FFT8(a); twiddle8(a, ii, 512); storex8(a, S1, 66);)
    PHASE(IDX loadx8(a, L1, 8);)
    PHASE(IDX storey8(a, S1, 66);)
    PHASE(IDX loady8(a, L1, 8);)
    PHASE(IDX FFT8(a); twiddle8(a, ii >> 3, 64); storex8(a, S1, 72);)
    PHASE(IDX loadx8(a, L2, 8);)
    PHASE(IDX storey8(a, S1, 72);)
    PHASE(IDX loady8(a, L2, 8);)
    PHASE(IDX FFT8(a); globalStores8(a, in, 64);)
#undef IDX
#undef S1
#undef L1
#undef L2
  }
}

static void launch_fft1D_1024(T2 *in_all, size_t n, int G) {
  enum { L = 128 };
  T2 A[L][16]; T smem[1088];
#pragma acc parallel loop gang ACC_SHAPE(G, 128) private(A, smem) present(in_all[0 : n])
  for (int group_id = 0; group_id < G; group_id++) {
#define IDX int ii = local_id; T2 *in = in_all + mad24(group_id, 1024, ii);
#define S1 (smem + ii)
#define L1 (smem + mad24(ii & 7, 130, ii >> 3))
#define L2 (smem + mad24((ii & 63) >> 3, 136, mad24(ii >> 6, 8, ii & 7)))
#define L3 (smem + mad24(ii >> 6, 256, ii & 63))
    PHASE(IDX globalLoads8(a, in, 128); FFT8(a); twiddle8(a, ii, 1024); storex8(a, S1, 130);)
    PHASE(IDX loadx8(a, L1, 16);)
    PHASE(IDX storey8(a, S1, 130);)
    PHASE(IDX loady8(a, L1, 16);)
    /* trans_s(a, S1, 136, L2, 32, 8, -112, 0) */
    PHASE(IDX FFT8(a); twiddle8(a, ii >> 3, 128); storex(a, S1, 136, 4, NO_REV, 0); storex(a, S1 + 0, 136, 4, NO_REV, 4);)
    PHASE(IDX loadx(a, L2, 32, 4, NO_REV, 0); loadx(a, L2 + (-112), 32, 4, NO_REV, 4);)
    PHASE(IDX storey(a, S1, 136, 4, NO_REV, 0); storey(a, S1 + 0, 136, 4, NO_REV, 4);)
    PHASE(IDX loady(a, L2, 32, 4, NO_REV, 0); loady(a, L2 + (-112), 32, 4, NO_REV, 4);)
    /* trans_s(a, S1, 256, L3, 64, 8, 256, -896) */
    PHASE(IDX FFT4S(a + 0); FFT4S(a + 4); twiddle4(a, ii >> 6, 16); twiddle4(a + 4, (128 + ii) >> 6, 16);
          storex(a, S1, 256, 4, NO_REV, 0); storex(a, S1 + (-896), 256, 4, NO_REV, 4);)
    PHASE(IDX loadx(a, L3, 64, 4, NO_REV, 0); loadx(a, L3 + 256, 64, 4, NO_REV, 4);)
    PHASE(IDX storey(a, S1, 256, 4, NO_REV, 0); storey(a, S1 + (-896), 256, 4, NO_REV, 4);)
    PHASE(IDX loady(a, L3, 64, 4, NO_REV, 0); loady(a, L3 + 256, 64, 4, NO_REV, 4);)
    PHASE(IDX FFT4S(a); FFT4S(a + 4); globalStores8(a, in, 128);)
#undef IDX
#undef S1
#undef L1
#undef L2
#undef L3
  }
}

static void launch_fft1D_2048(T2 *in_all, size_t n, int G) {
  enum { L = 256 };
  T2 A[L][16]; T smem[2112];
#pragma acc parallel loop gang ACC_SHAPE(G, 256) private(A, smem) present(in_all[0 : n])
  for (int group_id = 0; group_id < G; group_id++) {
#define IDX int ii = local_id; T2 *in = in_all + mad24(group_id, 2048, ii);
#define S1 (smem + ii)
#define L1 (smem + mad24(ii & 7, 258, ii >> 3))
#define L2 (smem + mad24((ii & 63) >> 3, 264, mad24(ii >> 6, 8, ii & 7)))
#define L3 (smem + mad24(ii >> 6, 256, ii & 63))
    PHASE(IDX globalLoads8(a, in, 256); FFT8(a); twiddle8(a, ii, 2048); storex8(a, S1, 258);)
    PHASE(IDX loadx8(a, L1, 32);)
    PHASE(IDX storey8(a, S1, 258);)
    PHASE(IDX loady8(a, L1, 32);)
    PHASE(IDX FFT8(a); twiddle8(a, ii >> 3, 256); storex8(a, S1, 264);)
    PHASE(IDX loadx8(a, L2, 32);)
    PHASE(IDX storey8(a, S1, 264);)
    PHASE(IDX loady8(a, L2, 32);)
    /* trans_s(a, S1, 256, L3, 64, 8, 768, 0) */
    PHASE(IDX FFT8(a); twiddle8(a, ii >> 6, 32); storex(a, S1, 256, 4, NO_REV, 0); storex(a, S1 + 0, 256, 4, NO_REV, 4);)
    PHASE(IDX loadx(a, L3, 64, 4, NO_REV, 0); loadx(a, L3 + 768, 64, 4, NO_REV, 4);)
    PHASE(IDX storey(a, S1, 256, 4, NO_REV, 0); storey(a, S1 + 0, 256, 4, NO_REV, 4);)
    PHASE(IDX loady(a, L3, 64, 4, NO_REV, 0); loady(a, L3 + 768, 64, 4, NO_REV, 4);)
    PHASE(IDX FFT4S(a); FFT4S(a + 4); globalStores8_r(a, in, 256);)
#undef IDX
#undef S1
#undef L1
#undef L2
#undef L3
  }
}

static void launch_fft1D_256(T2 *in_all, size_t n, int G) {
  enum { L = 64 };
  T2 A[L][16]; T smem[640];
#pragma acc parallel loop gang ACC_SHAPE(G, 64) private(A, smem) present(in_all[0 : n])
  for (int group_id = 0; group_id < G; group_id++) {
#define IDX int ii = local_id & 31; int jj = local_id >> 5;
#define S1 (smem + mad24(jj, 272, ii))
#define L1 (smem + mad24(ii & 7, 34, mad24(jj, 272, ii >> 3)))
#define S2 (smem + mad24(jj, 320, ii))
#define L2 (smem + mad24(ii >> 3, 40, mad24(jj, 320, ii & 7)))
    /* rows jj != 0 load nothing (uninitialised registers in OpenCL; zeroed here) and store nothing */
    PHASE(IDX
          _Pragma("acc loop seq") for (int r = 0; r < 8; r++) a[r] = mk2(0, 0);
          if (jj == 0) globalLoads8(a, in_all + mad24(mad24(group_id, 2, jj), 256, ii), 32);
          FFT8(a); twiddle8(a, ii, 256); storex8(a, S1, 34);)
    PHASE(IDX loadx8(a, L1, 4);)
    PHASE(IDX storey8(a, S1, 34);)
    PHASE(IDX loady8(a, L1, 4);)
    /* trans_s(a, S2, 40, L2, 8, 8, 128, 0) */
    PHASE(IDX FFT8(a); twiddle8(a, ii >> 3, 32); storex(a, S2, 40, 4, NO_REV, 0); storex(a, S2 + 0, 40, 4, NO_REV, 4);)
    PHASE(IDX loadx(a, L2, 8, 4, NO_REV, 0); loadx(a, L2 + 128, 8, 4, NO_REV, 4);)
    PHASE(IDX storey(a, S2, 40, 4, NO_REV, 0); storey(a, S2 + 0, 40, 4, NO_REV, 4);)
    PHASE(IDX loady(a, L2, 8, 4, NO_REV, 0); loady(a, L2 + 128, 8, 4, NO_REV, 4);)
    /* [FIX-OCL-25] only the jj == 0 row stores */
    PHASE(IDX FFT4S(a); FFT4S(a + 4); if (jj == 0) globalStores8_r(a, in_all + mad24(mad24(group_id, 2, jj), 256, ii), 32);)
#undef IDX
#undef S1
#undef L1
#undef S2
#undef L2
  }
}

static void launch_fft1D_128(T2 *in_all, size_t n, int G) {
  enum { L = 64 };
  T2 A[L][16]; T smem[640];
#pragma acc parallel loop gang ACC_SHAPE(G, 64) private(A, smem) present(in_all[0 : n])
  for (int group_id = 0; group_id < G; group_id++) {
#define IDX int ii = local_id & 15; int jj = local_id >> 4;
#define S1 (smem + mad24(jj, 144, ii))
#define L1 (smem + mad24(ii & 7, 18, mad24(jj, 144, ii >> 3)))
#define S2 (smem + mad24(jj, 160, ii))
#define L2 (smem + mad24(ii >> 3, 40, mad24(jj, 160, ii & 7)))
    /* trans_s(a, S1, 18, L1, 4, 8, -14, 0) */
    PHASE(IDX
          _Pragma("acc loop seq") for (int r = 0; r < 8; r++) a[r] = mk2(0, 0);
          if (jj == 0) globalLoads8(a, in_all + mad24(mad24(group_id, 4, jj), 128, ii), 16);
          FFT8(a); twiddle8(a, ii, 128); storex(a, S1, 18, 4, NO_REV, 0); storex(a, S1 + 0, 18, 4, NO_REV, 4);)
    PHASE(IDX loadx(a, L1, 4, 4, NO_REV, 0); loadx(a, L1 + (-14), 4, 4, NO_REV, 4);)
    PHASE(IDX storey(a, S1, 18, 4, NO_REV, 0); storey(a, S1 + 0, 18, 4, NO_REV, 4);)
    PHASE(IDX loady(a, L1, 4, 4, NO_REV, 0); loady(a, L1 + (-14), 4, 4, NO_REV, 4);)
    /* trans_s(a, S2, 40, L2, 8, 8, 48, -144) */
    PHASE(IDX FFT4S(a); FFT4S(a + 4); twiddle4(a, ii >> 3, 16); twiddle4(a + 4, (16 + ii) >> 3, 16);
          storex(a, S2, 40, 4, NO_REV, 0); storex(a, S2 + (-144), 40, 4, NO_REV, 4);)
    PHASE(IDX loadx(a, L2, 8, 4, NO_REV, 0); loadx(a, L2 + 48, 8, 4, NO_REV, 4);)
    PHASE(IDX storey(a, S2, 40, 4, NO_REV, 0); storey(a, S2 + (-144), 40, 4, NO_REV, 4);)
    PHASE(IDX loady(a, L2, 8, 4, NO_REV, 0); loady(a, L2 + 48, 8, 4, NO_REV, 4);)
    /* [FIX-OCL-25] only the jj == 0 row stores */
    PHASE(IDX FFT4S(a); FFT4S(a + 4); if (jj == 0) globalStores8_r(a, in_all + mad24(mad24(group_id, 4, jj), 128, ii), 16);)
#undef IDX
#undef S1
#undef L1
#undef S2
#undef L2
  }
}

/* first (global transpose) stage for N = 4096, 8192 */
static void launch_fft0(T2 *in_all, T2 *out_all, size_t n, int G, int fftn1, int pow1) {
  enum { L = 128 };
  T2 A[L][16]; T smem[2064];
#pragma acc parallel loop gang ACC_SHAPE(G, 128) private(A, smem) present(in_all[0 : n], out_all[0 : n])
  for (int group_id = 0; group_id < G; group_id++) {
#define IDX int b = group_id & (fftn1 / 2048 - 1); int d = ocl_shr(group_id, pow1 - 11); int tid = local_id; int i = tid & 15; int j = tid >> 4;
#define S1 (smem + tid)
#define L1 (smem + mad24(j, 256, i))
#define S2 (smem + mad24(i, 129, j << 1))
#define L2 (smem + mad24(tid >> 7, 129, tid & 127))
    PHASE(IDX
          int Offset_in = mul24(b, 16); Offset_in += ocl_shl(d, pow1); Offset_in += mad24(j, fftn1 / 128, i);
          globalLoads16(a, in_all + Offset_in, fftn1 / 16); FFT16(a); twiddle16(a, j, 128); storex16(a, S1, 128);)
    PHASE(IDX loadx16(a, L1, 16);)
    PHASE(IDX storey16(a, S1, 128);)
    PHASE(IDX loady16(a, L1, 16);)
    /* trans_s(a, S2, 16, L2, 129, 16, 0, -127) */
    PHASE(IDX FFT8(a); FFT8(a + 8); twiddle_s(a, (b << 4) + i, fftn1, 16, j << 1);
          storex(a, S2, 16, 8, NO_REV, 0); storex(a, S2 + (-127), 16, 8, NO_REV, 8);)
    PHASE(IDX loadx(a, L2, 129, 8, NO_REV, 0); loadx(a, L2 + 0, 129, 8, NO_REV, 8);)
    PHASE(IDX storey(a, S2, 16, 8, NO_REV, 0); storey(a, S2 + (-127), 16, 8, NO_REV, 8);)
    PHASE(IDX loady(a, L2, 129, 8, NO_REV, 0); loady(a, L2 + 0, 129, 8, NO_REV, 8);)
    PHASE(IDX
          int Offset_out = mad24(mul24(b, 16), 128, 0); Offset_out += ocl_shl(d, pow1); Offset_out += tid;
          globalStores16(a, out_all + Offset_out, 128);)
#undef IDX
#undef S1
#undef L1
#undef S2
#undef L2
  }
}

static void launch_fft1D_8192(T2 *in_all, T2 *out_all, size_t n, int G) {
  enum { L = 128 };
  T2 A[L][16]; T smem[1024];
#pragma acc parallel loop gang ACC_SHAPE(G, 128) private(A, smem) present(in_all[0 : n], out_all[0 : n])
  for (int group_id = 0; group_id < G; group_id++) {
#define IDX int b = group_id & 7; int d = group_id >> 3; int tid = local_id; int i = tid & 15; int j = tid >> 4;
#define S1 (smem + tid)
#define L1 (smem + mad24(j, 128, i))
    PHASE(IDX
          int Offset_in = mul24(b, 16); Offset_in += (d << 13); Offset_in += mad24(j, 128, i);
          globalLoads8(a, in_all + Offset_in, 1024); FFT8(a); twiddle8(a, j, 64); storex8(a, S1, 128);)
    PHASE(IDX loadx8(a, L1, 16);)
    PHASE(IDX storey8(a, S1, 128);)
    PHASE(IDX loady8(a, L1, 16);)
    PHASE(IDX
          int tid0 = mul24(b, 16); int Offset_out = mad24(tid0 >> 7, 8192, tid0 & 127); Offset_out += (d << 13);
          Offset_out += mad24(j, 128, i);
          FFT8(a); globalStores8(a, out_all + Offset_out, 1024);)
#undef IDX
#undef S1
#undef L1
  }
}

static void launch_fft1D_4096(T2 *in_all, T2 *out_all, size_t n, int G) {
  enum { L = 64 };
  T2 A[L][16]; T smem[512];
#pragma acc parallel loop gang ACC_SHAPE(G, 64) private(A, smem) present(in_all[0 : n], out_all[0 : n])
  for (int group_id = 0; group_id < G; group_id++) {
#define IDX int b = group_id & 7; int d = group_id >> 3; int tid = local_id; int i = tid & 15; int j = tid >> 4;
#define S1 (smem + tid)
#define L1 (smem + mad24(j, 128, i))
    PHASE(IDX
          int Offset_in = mul24(b, 16); Offset_in += (d << 12); Offset_in += mad24(j, 128, i);
          globalLoads8(a, in_all + Offset_in, 512); FFT8(a); twiddle8(a, j, 32); storex8(a, S1, 64);)
    PHASE(IDX loadx8(a, L1, 16);)
    PHASE(IDX storey8(a, S1, 64);)
    PHASE(IDX loady8(a, L1, 16);)
    PHASE(IDX
          int tid0 = mul24(b, 16); int Offset_out = mad24(tid0 >> 7, 4096, tid0 & 127); Offset_out += (d << 12);
          Offset_out += mad24(j, 256, i);
          FFT4S(a); FFT4S(a + 4);
          T2 *out = out_all + Offset_out;
          out[0] = a[0]; out[1024] = a[1]; out[2048] = a[2]; out[3072] = a[3];
          out[128] = a[4]; out[1152] = a[5]; out[2176] = a[6]; out[3200] = a[7];)
#undef IDX
#undef S1
#undef L1
  }
}

/* ---------------- host: forward() of fftlib.cpp for the 1D path ---------------- */
static void forward(T2 *work, T2 *temp, int fftn) {
  size_t n = (size_t)fftn;
  int pow1 = 0; while ((1 << pow1) < fftn) pow1++; /* log2(fftn): [FIX-OCL-24] */
  switch (fftn) {
    case 128: launch_fft1D_128(work, n, 1); break;
    case 256: launch_fft1D_256(work, n, 1); break;
    case 512: launch_fft1D_512(work, n, 1); break;
    case 1024: launch_fft1D_1024(work, n, 1); break;
    case 2048: launch_fft1D_2048(work, n, 1); break;
    case 4096:
      launch_fft0(work, temp, n, 2, fftn, pow1); /* stage 0: 2 groups x 128 */
      launch_fft1D_4096(temp, work, n, 8);       /* stage 1: 8 groups x 64 */
      break;
    case 8192:
      launch_fft0(work, temp, n, 4, fftn, pow1);
      launch_fft1D_8192(temp, work, n, 8);       /* stage 1: 8 groups x 128 */
      break;
  }
}

static int cmp_double(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return (x > y) - (x < y);
}

int main(int argc, char **argv) {
  int sizeIndex = 0, passes = 10;
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--pts") && i + 1 < argc) sizeIndex = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--passes") && i + 1 < argc) passes = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--2D")) { fprintf(stderr, "--2D is not translated in the OpenACC tree\n"); return 3; }
  }
  (void)passes;
  int probSizes[7] = {128, 256, 512, 1024, 2048, 4096, 8192};
  sizeIndex -= 1;
  if (sizeIndex < 0 || sizeIndex >= 7) { fprintf(stderr, "Invalid size index specified\n"); return 1; }
  int fftn = probSizes[sizeIndex];
  unsigned long used_bytes = fftn * sizeof(T2);
  double N = fftn;
  fprintf(stdout, "used_bytes=%lu, N=%g\n", used_bytes, N);
  T2 *source = (T2 *)malloc(used_bytes), *result = (T2 *)malloc(used_bytes);
  T2 *work = (T2 *)malloc(used_bytes), *temp = (T2 *)malloc(used_bytes);
  srand(1); /* [FIX-OCL-26] explicit seed */
  for (int i = 0; i < N; i++) {
    source[i].x = (rand() / (float)RAND_MAX) * 2 - 1;
    source[i].y = (rand() / (float)RAND_MAX) * 2 - 1;
  }
  size_t n = (size_t)fftn;
  acc_init(acc_get_device_type()); /* device initialisation, outside the timed region */
  dwarf_init("fft", "spectral-methods", "openacc");
  dwarf_set_size(fftn, 1);
  for (int rep = dwarf_first_rep(); rep < dwarf_reps(); ++rep) {
    dwarf_region_begin();
    /* as OpenCL: device allocation inside the region, release outside it */
#pragma acc enter data create(work[0 : n], temp[0 : n])
    memcpy(work, source, used_bytes);
#pragma acc update device(work[0 : n])
    forward(work, temp, fftn);
#pragma acc update self(work[0 : n])
    dwarf_region_end(rep);
#pragma acc exit data delete(work[0 : n], temp[0 : n])
  }
  memcpy(result, work, used_bytes);
#ifdef PRINT_RESULT
  for (int i = 0; i < N; i++) fprintf(stdout, "data[%d] (%g, %g)\n", i, result[i].x, result[i].y);
#endif

  /* validation (FIX-OCL-27): sorted magnitudes of the output vs the double-precision DFT of the same input (insensitive to the
   * radix-permuted output order and to swapped real/imaginary parts), relative RMS error < 1e-3. */
  double *refm = (double *)malloc(sizeof(double) * fftn), *gotm = (double *)malloc(sizeof(double) * fftn);
  for (int k = 0; k < fftn; k++) {
    double sr = 0, si = 0;
    for (int t = 0; t < fftn; t++) {
      double ang = -2.0 * M_PI * (double)((long)k * t % fftn) / fftn;
      sr += source[t].x * cos(ang) - source[t].y * sin(ang);
      si += source[t].x * sin(ang) + source[t].y * cos(ang);
    }
    refm[k] = sqrt(sr * sr + si * si);
    gotm[k] = sqrt((double)result[k].x * result[k].x + (double)result[k].y * result[k].y);
  }
  qsort(refm, fftn, sizeof(double), cmp_double);
  qsort(gotm, fftn, sizeof(double), cmp_double);
  double e = 0, nrm = 0;
  for (int k = 0; k < fftn; k++) { e += (refm[k] - gotm[k]) * (refm[k] - gotm[k]); nrm += refm[k] * refm[k]; }
  double rel = sqrt(e / nrm);
  int ok = (rel < 1e-3) && (rel == rel);
  printf("%s: sorted-magnitude relative RMS error vs the double-precision DFT of the input: %.3g\n", ok ? "SUCCESS" : "FAILURE", rel);
  return dwarf_report(ok);
}
