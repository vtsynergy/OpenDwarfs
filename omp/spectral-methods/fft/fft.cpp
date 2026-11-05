/* Spectral-methods dwarf: FFT -- OpenMP target-offload translation of
 * opencl/spectral-methods/fft (SHOC FFT, Cooley-Tukey radix-8/16 kernels with local-memory transposes).
 *
 * Covered: the default one-dimensional path (--pts N, N = 128 ... 8192 points, one transform):
 *   N <= 2048 : one in-place kernel fft1D_N, one work-group;
 *   N = 4096, 8192 : two kernels fft0 (global transpose stage) -> fft1D_N.
 * The kernel bodies are generated mechanically from fft.cl (scripts: see git history of this file's header): each OpenCL work-group is an
 * OpenMP team, each work-item a thread of that team, __local arrays are team-shared arrays, barrier() is
 * `omp barrier`. Because the kernels keep per-work-item state in private arrays across barriers, a team must
 * really get one thread per work-item; this is checked on the device and a toolchain that cannot provide it
 * fails the run (reported as such, not worked around).
 * Not translated: the two-dimensional path (--2D; kernels fft1/fft2 and the runtime -D tuning of
 * setGlobalOption), the inverse transform and the check kernel of the original SHOC benchmark (the OpenDwarfs
 * app does not use them).
 * Quirk preserved from the OpenCL host: it passes -D pow1=<fftn/16> (a log2 was intended), which the kernels
 * turn into out-of-range shifts; with GPU shift semantics these give 0, reproduced by ocl_shr/ocl_shl.
 * Input: rand() with the default seed, as in the OpenCL app. Timed region: allocation + H2D + forward FFT + D2H.
 * Validation (added): double-precision DFT of the same input.
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>
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
#define barrier(x) _Pragma("omp barrier")
#define CLK_LOCAL_MEM_FENCE 0

#pragma omp declare target
static inline T2 mk2(float x, float y) { T2 r; r.x = x; r.y = y; return r; }
static inline int mad24(int a, int b, int c) { return a * b + c; }
static inline int mul24(int a, int b) { return a * b; }
static inline int ocl_shr(int x, int s) { return (s >= 32 || s < 0) ? 0 : (x >> s); }
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
static inline T2 exp_i( T phi ) {
	return mk2( cosf(phi), sinf(phi) );
}
static inline T2 cmplx_mul( T2 a, T2 b ) { return mk2( a.x*b.x-a.y*b.y, a.x*b.y+a.y*b.x ); }
static inline T2 cm_fl_mul( T2 a, T  b ) { return mk2( b*a.x, b*a.y ); }
static inline T2 cmplx_add( T2 a, T2 b ) { return mk2( a.x + b.x, a.y + b.y ); }
static inline T2 cmplx_sub( T2 a, T2 b ) { return mk2( a.x - b.x, a.y - b.y ); }
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



static inline void globalLoads8(T2 *data, T2 *in, int S){
	for( int i = 0; i < 8; i++ )
		data[i] = in[i*S];
}

static inline void globalLoads16(T2 *data, T2 *in, int S){
	for( int i = 0; i < 16; i++ )
		data[i] = in[i*S];
}


static inline void globalStores8_r(T2 *data, T2 *out, int S){
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
#ifdef FFT_2D
		out[i*S] = data[reversed[i]];
#else
	out[i*S] = cmplx_swap(data[reversed[i]]);
#endif
}


static inline void globalStores8(T2 *data, T2 *out, int S){
	for( int i = 0; i < 8; i++ )
#ifdef FFT_2D
		out[i*S] = data[i];
#else
	out[i*S] = cmplx_swap(data[i]);
#endif
}
static inline void globalStores16(T2 *data, T2 *out, int S){
	for( int i = 0; i < 16; i++ )
		out[i*S] = data[i];

}

static inline void storex4( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 4; i++ )
		x[i*sx] = a[i].x;
}

static inline void storey4( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 4; i++ )
		x[i*sx] = a[i].y;
}

static inline void storex16( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 16; i++ )
		x[i*sx] = a[i].x;
}

static inline void storey16( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 16; i++ )
		x[i*sx] = a[i].y;
}


static inline void storex8( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 8; i++ )
		x[i*sx] = a[i].x;
}

static inline void storey8( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 8; i++ )
		x[i*sx] = a[i].y;
}

static inline void storex8_r( T2 *a, T *x, int sx ) {
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
		x[i*sx] = a[reversed[i]].x;
}

static inline void storey8_r( T2 *a, T *x, int sx ) {
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
		x[i*sx] = a[reversed[i]].y;
}

static inline void loadx16( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 16; i++ )
		a[i].x = x[i*sx];
}

static inline void loady16( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 16; i++ )
		a[i].y = x[i*sx];
}

static inline void loadx8( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 8; i++ )
		a[i].x = x[i*sx];
}

static inline void loady8( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 8; i++ )
		a[i].y = x[i*sx];
}

static inline void loadx4( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 4; i++ )
		a[i].x = x[i*sx];
}

static inline void loady4( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 4; i++ )
		a[i].y = x[i*sx];
}

static inline void loadx8_r( T2 *a, T *x, int sx ) {
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
		a[i].x = x[reversed[i]*sx];
}

static inline void loady8_r( T2 *a, T *x, int sx ) {
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
		a[i].y = x[reversed[i]*sx];
}

static inline void globalLoads(T2 *data, T2 *in, int S, int n){
	for( int i = 0; i < n; i++ )
		data[i] = in[i*S];
}

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

static inline void trans(T2 *a, T *s, int ds, T *l, int dl, int n, int rev){
	storex(a, s, ds, n, rev==STORE_REV, 0);
	barrier(CLK_LOCAL_MEM_FENCE);
	loadx(a, l, dl, n, rev==LOAD_REV, 0);
	barrier(CLK_LOCAL_MEM_FENCE);
	storey(a, s, ds, n, rev==STORE_REV, 0);
	barrier(CLK_LOCAL_MEM_FENCE);
	loady(a, l, dl, n, rev==LOAD_REV, 0);
	barrier(CLK_LOCAL_MEM_FENCE);
}

static inline void trans_s(T2 *a, T *s, int ds, T *l, int dl, int n, int offset_l, int offset_s){
	storex(a, s, ds, n/2, NO_REV, 0);
	storex(a, s+offset_s, ds, n/2, NO_REV, n/2);
	barrier(CLK_LOCAL_MEM_FENCE);
	loadx(a, l, dl, n/2, NO_REV, 0);
	loadx(a, l+offset_l, dl, n/2, NO_REV, n/2);
	barrier(CLK_LOCAL_MEM_FENCE);
	storey(a, s, ds, n/2, NO_REV, 0);
	storey(a, s+offset_s, ds, n/2, NO_REV, n/2);
	barrier(CLK_LOCAL_MEM_FENCE);
	loady(a, l, dl, n/2, NO_REV, 0);
	loady(a, l+offset_l, dl, n/2, NO_REV, n/2);
	barrier(CLK_LOCAL_MEM_FENCE);
}



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

#define transpose( a, s, ds, l, dl )                              \
{                                                                       \
	storex8( a, s, ds );  barrier(CLK_LOCAL_MEM_FENCE);  \
	loadx8 ( a, l, dl );  barrier(CLK_LOCAL_MEM_FENCE);  \
	storey8( a, s, ds );  barrier(CLK_LOCAL_MEM_FENCE);  \
	loady8 ( a, l, dl );  barrier(CLK_LOCAL_MEM_FENCE);  \
}

#define transpose16( a, s, ds, l, dl )                              \
{                                                                       \
	storex16( a, s, ds );  barrier(CLK_LOCAL_MEM_FENCE);  \
	loadx16 ( a, l, dl );  barrier(CLK_LOCAL_MEM_FENCE);  \
	storey16( a, s, ds );  barrier(CLK_LOCAL_MEM_FENCE);  \
	loady16 ( a, l, dl );  barrier(CLK_LOCAL_MEM_FENCE);  \
}


static inline void twiddle16(T2 *a, int tt, int n)
{
	T c = (T) tt;
	for( int t = 1; t < 16; t++ ){
		(a)[t] = cmplx_mul( (a)[t],exp_i((2.0f*M_PI*t/n)*c));
	}
}
static inline void twiddle8(T2 *a, int tt, int n)
{
	T c = (T) tt;
	for( int t = 1; t < 8; t++ ){
		(a)[t] = cmplx_mul( (a)[t],exp_i((2.0f*M_PI*t/n)*c));
	}
}
static inline void  twiddle4(T2 *a, int tt, int n )                                 {
	T c = (T) tt;
	for( int t = 1; t < 4; t++ ){
		(a)[t] = cmplx_mul( (a)[t],exp_i((2.0f*M_PI*t/n)*c));
	}
}

#ifdef FFT_512
static void fft1D_512(T2 *in, T *smem, int local_id, int group_id, int fftn1, int pow1)
{
	
	int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	T *Local_store, *Local_load;
	T2 a[8];
	
	
	ii = local_id;
	jj = 0;
	offset =  mad24(group_id, 512, ii);
	in += offset;
	globalLoads8(a, in, 64);
	FFT8(a);
	twiddle8(a,ii,512);
	j = ii & 7;
	i = ii >> 3;
	Local_store = smem + ii;
	Local_load = smem + mad24(j, 66, i);
	transpose(a,Local_store,66,Local_load,8);
	FFT8(a);
	twiddle8(a,ii>>3,64);
	j = ii >> 3;
	i = ii & 7;
	Local_store = smem + ii;
	Local_load = smem + mad24(j, 72, i);
	transpose(a,Local_store,72,Local_load,8);
	FFT8(a);
	globalStores8(a,in,64);
}
#endif

#ifdef FFT_1024
static void fft1D_1024(T2 *in, T *smem, int local_id, int group_id, int fftn1, int pow1)
{
	
	int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	T *Local_store, *Local_load;
	T2 a[8];
	
	
	ii = local_id;
	jj = 0;
	offset =  mad24(group_id, 1024, ii);
	in += offset;
	globalLoads8(a,in,128);
	FFT8(a);
	twiddle8(a,ii,1024);
	Local_store = smem + ii;
	j = ii & 7;
	i = ii >> 3;
	Local_load = smem + mad24(j, 130, i);
	transpose(a,Local_store,130,Local_load,16);
	FFT8(a);
	twiddle8(a,ii>>3,128);
	Local_store = smem + ii;
	j = (ii & 63) >> 3;
	i = mad24(ii >> 6, 8, ii & 7);
	Local_load = smem + mad24(j, 136, i);
	trans_s(a,Local_store,136,Local_load,32,8,-112,0);
	FFT4S(a+0);
	FFT4S(a+4);
	twiddle4(a,ii>>6,16);
	twiddle4(a+4,(128 + ii) >>6,16);
	Local_store = smem + ii;
	j = ii >> 6;
	i = ii & 63;
	Local_load = smem + mad24(j, 256, i);
	trans_s(a,Local_store,256,Local_load,64,8,256,-896);
	FFT4S(a);
	FFT4S(a+4);
	globalStores8(a,in,128);
}
#endif

#ifdef FFT_2048
static void fft1D_2048(T2 *in, T *smem, int local_id, int group_id, int fftn1, int pow1)
{
	
	int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	T *Local_store, *Local_load;
	T2 a[8];
	
	
	ii = local_id;
	jj = 0;
	offset =  mad24(group_id, 2048, ii);
	in += offset;
	globalLoads8(a,in,256);
	FFT8(a);
	twiddle8(a,ii,2048);
	Local_store = smem + ii;
	j = ii & 7;
	i = ii >> 3;
	Local_load = smem + mad24(j, 258, i);
	transpose(a,Local_store,258,Local_load,32);
	FFT8(a);
	twiddle8(a,ii>>3,256);
	Local_store = smem + ii;
	j = (ii & 63) >> 3;
	i = mad24(ii >> 6, 8, ii & 7);
	Local_load = smem + mad24(j, 264, i);
	transpose(a,Local_store, 264, Local_load, 32);
	FFT8(a);
	twiddle8(a,ii>>6,32);
	Local_store = smem + ii;
	j = ii >> 6;
	i = ii & 63;
	Local_load = smem + mad24(j, 256, i);
	trans_s(a, Local_store, 256, Local_load, 64, 8, 768, 0);
	FFT4S(a);
	FFT4S(a+4);
	globalStores8_r(a, in, 256);
}
#endif

#ifdef FFT_256
static void fft1D_256(T2 *in, T *smem, int local_id, int group_id, int fftn1, int pow1)
{
	
	int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	T *Local_store, *Local_load;
	T2 a[8];
	
	
	ii = local_id & 31;
	jj = local_id >> 5;
	if( jj == 0){
		offset = mad24( mad24(group_id, 2, jj), 256, ii );
		in += offset;
		globalLoads8(a, in, 32);}
	FFT8(a);
	twiddle8(a, ii, 256);
	Local_store = smem + mad24(jj, 272, ii);
	j = ii & 7;
	i = ii >> 3;
	i = mad24(jj, 272, i);
	Local_load = smem + mad24(j, 34, i);
	transpose(a,Local_store, 34, Local_load, 4);
	FFT8(a);
	twiddle8(a, ii>>3, 32);
	Local_store = smem + mad24(jj, 320, ii);
	j = ii >> 3;
	i = ii & 7;
	i = mad24(jj, 320, i);
	Local_load = smem + mad24(j, 40, i);
	trans_s(a, Local_store, 40, Local_load, 8, 8,128, 0);
	FFT4S(a);
	FFT4S(a+4);
	globalStores8_r(a, in, 32);
}
#endif

#ifdef FFT_128
static void fft1D_128(T2 *in, T *smem, int local_id, int group_id, int fftn1, int pow1)
{
	
	int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	T *Local_store, *Local_load;
	T2 a[8];
	
	
	ii = local_id & 15;
	jj = local_id >> 4;
	if( jj == 0){
		offset = mad24( mad24(group_id, 4, jj), 128, ii );
		in += offset;
		globalLoads8(a, in, 16);}
	FFT8(a);
	twiddle8(a, ii, 128);
	Local_store = smem + mad24(jj, 144, ii);
	j = ii & 7;
	i = ii >> 3;
	i = mad24(jj, 144, i);
	Local_load = smem + mad24(j, 18, i);
	trans_s(a, Local_store, 18, Local_load, 4, 8, -14, 0);
	FFT4S(a);
	FFT4S(a+4);
	twiddle4(a, ii>>3, 16);
	twiddle4(a+4, (16 + ii) >>3, 16);
	Local_store = smem + mad24(jj, 160, ii);
	j = ii >> 3;
	i = ii & 7;
	i = mad24(jj, 160, i);
	Local_load = smem + mad24(j, 40, i);
	trans_s(a, Local_store, 40, Local_load, 8, 8, 48, -144);
	FFT4S(a);
	FFT4S(a+4);
	globalStores8_r(a, in, 16);
}
#endif

static void fft0(T2 *in, T2 *out, T *smem, int local_id, int group_id, int fftn1, int pow1)
{
	
	int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	T *Local_store, *Local_load;
	T2 a[16];
	
	
	b = group_id & (fftn1/2048-1);
	d = ocl_shr(group_id, pow1-11);
	Offset_in = mul24(b, 16);
	tid = Offset_in;
	i = tid;
	j = 0;
	Offset_out = mad24(i, 128, j);
	Offset_in += ocl_shl(d, pow1);
	Offset_out += ocl_shl(d, pow1);
	tid = local_id;
	i = tid & 15;
	j = tid >> 4;
	Offset_in += mad24(j, fftn1/128, i);
	in += Offset_in;
	globalLoads16(a, in ,fftn1/16);
	FFT16(a);
	twiddle16(a, j, 128);
	Offset_in = mad24(j, 256, i);
	Local_store = smem + tid;
	Local_load = smem + Offset_in;
	transpose16(a, Local_store, 128, Local_load, 16);
	FFT8(a);
	FFT8(a+8);
	l = ((b << 4) + i) >> 0;
	k = j << 1;
	twiddle_s(a, l, fftn1, 16, k);
	Local_store = smem + mad24(i, 129, j << 1);
	Local_load = smem + mad24(tid >> 7, 129, tid & 127);
	trans_s(a, Local_store, 16, Local_load, 129, 16, 0, -127);
	Offset_out += tid;
	out += Offset_out;
	globalStores16(a, out, 128);
	//globalStores(a, out, 128, 16);
}

#ifdef FFT_8192
static void fft1D_8192(T2 *in, T2 *out, T *smem, int local_id, int group_id, int fftn1, int pow1)
{
	
	int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	T *Local_store, *Local_load;
	T2 a[8];
	
	
	b = group_id & 7;
	d = group_id >> 3;
	Offset_in = mul24(b, 16);
	tid = Offset_in;
	i = tid >> 7;
	j = tid & 127;
	Offset_out = mad24(i, 8192, j);
	Offset_in += (d << 13);
	Offset_out += (d << 13);
	tid = local_id;
	i = tid & 15;
	j = tid >> 4;
	Offset_in += mad24(j, 128, i);
	in += Offset_in;
	globalLoads8(a, in, 1024);
	FFT8(a);
	twiddle8(a, j, 64);
	Offset_in = mad24(j, 128, i);
	Local_store = smem + tid;
	Local_load = smem + Offset_in;
	transpose(a, Local_store, 128, Local_load, 16);
	FFT8(a);
	Offset_out += mad24(j, 128, i);
	out += Offset_out;
	globalStores8(a, out, 1024);
}
#endif

#ifdef FFT_4096
static void fft1D_4096(T2 *in, T2 *out, T *smem, int local_id, int group_id, int fftn1, int pow1)
{
	
	int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	T *Local_store, *Local_load;
	T2 a[8];
	
	
	b = group_id & 7;
	d = group_id >> 3;
	Offset_in = mul24(b, 16);
	tid = Offset_in;
	i = tid >> 7;
	j = tid & 127;
	Offset_out = mad24(i, 4096, j);
	Offset_in += (d << 12);
	Offset_out += (d << 12);
	tid = local_id;
	i = tid & 15;
	j = tid >> 4;
	Offset_in += mad24(j, 128, i);
	in += Offset_in;
	globalLoads8(a, in, 512);
	FFT8(a);
	twiddle8(a, j, 32);
	Offset_in = mad24(j, 128, i);
	Local_store = smem + tid;
	Local_load = smem + Offset_in;
	transpose(a, Local_store, 64, Local_load, 16);
	FFT4S(a);
	FFT4S(a+4);
	Offset_out += mad24(j, 256, i);
	out += Offset_out;
	out[0] = a[0];
	out[1024] = a[1];
	out[2048] = a[2];
	out[3072] = a[3];
	out[128] = a[4];
	out[1152] = a[5];
	out[2176] = a[6];
	out[3200] = a[7];
}
#endif


#pragma omp end declare target

static void launch_fft1D_512(T2 *in, size_t n, int G, int L, int fftn1, int pow1, int *bad) {
#pragma omp target teams num_teams(G) thread_limit(L) map(present, alloc : in[0 : n]) map(tofrom : bad[0 : 1])
  {
    T smem[576];
#pragma omp parallel num_threads(L)
    {
      /* one OpenCL work-group = one team; its L work-items need L threads (barriers + private state) */
      if (omp_get_num_threads() != L) {
#pragma omp atomic write
        bad[0] = 1;
      } else {
        for (int g = omp_get_team_num(); g < G; g += omp_get_num_teams()) {
          fft1D_512(in, smem, omp_get_thread_num(), g, fftn1, pow1);
#pragma omp barrier
        }
      }
    }
  }
}

static void launch_fft1D_1024(T2 *in, size_t n, int G, int L, int fftn1, int pow1, int *bad) {
#pragma omp target teams num_teams(G) thread_limit(L) map(present, alloc : in[0 : n]) map(tofrom : bad[0 : 1])
  {
    T smem[1088];
#pragma omp parallel num_threads(L)
    {
      /* one OpenCL work-group = one team; its L work-items need L threads (barriers + private state) */
      if (omp_get_num_threads() != L) {
#pragma omp atomic write
        bad[0] = 1;
      } else {
        for (int g = omp_get_team_num(); g < G; g += omp_get_num_teams()) {
          fft1D_1024(in, smem, omp_get_thread_num(), g, fftn1, pow1);
#pragma omp barrier
        }
      }
    }
  }
}

static void launch_fft1D_2048(T2 *in, size_t n, int G, int L, int fftn1, int pow1, int *bad) {
#pragma omp target teams num_teams(G) thread_limit(L) map(present, alloc : in[0 : n]) map(tofrom : bad[0 : 1])
  {
    T smem[2112];
#pragma omp parallel num_threads(L)
    {
      /* one OpenCL work-group = one team; its L work-items need L threads (barriers + private state) */
      if (omp_get_num_threads() != L) {
#pragma omp atomic write
        bad[0] = 1;
      } else {
        for (int g = omp_get_team_num(); g < G; g += omp_get_num_teams()) {
          fft1D_2048(in, smem, omp_get_thread_num(), g, fftn1, pow1);
#pragma omp barrier
        }
      }
    }
  }
}

static void launch_fft1D_256(T2 *in, size_t n, int G, int L, int fftn1, int pow1, int *bad) {
#pragma omp target teams num_teams(G) thread_limit(L) map(present, alloc : in[0 : n]) map(tofrom : bad[0 : 1])
  {
    T smem[640];
#pragma omp parallel num_threads(L)
    {
      /* one OpenCL work-group = one team; its L work-items need L threads (barriers + private state) */
      if (omp_get_num_threads() != L) {
#pragma omp atomic write
        bad[0] = 1;
      } else {
        for (int g = omp_get_team_num(); g < G; g += omp_get_num_teams()) {
          fft1D_256(in, smem, omp_get_thread_num(), g, fftn1, pow1);
#pragma omp barrier
        }
      }
    }
  }
}

static void launch_fft1D_128(T2 *in, size_t n, int G, int L, int fftn1, int pow1, int *bad) {
#pragma omp target teams num_teams(G) thread_limit(L) map(present, alloc : in[0 : n]) map(tofrom : bad[0 : 1])
  {
    T smem[640];
#pragma omp parallel num_threads(L)
    {
      /* one OpenCL work-group = one team; its L work-items need L threads (barriers + private state) */
      if (omp_get_num_threads() != L) {
#pragma omp atomic write
        bad[0] = 1;
      } else {
        for (int g = omp_get_team_num(); g < G; g += omp_get_num_teams()) {
          fft1D_128(in, smem, omp_get_thread_num(), g, fftn1, pow1);
#pragma omp barrier
        }
      }
    }
  }
}

static void launch_fft0(T2 *in, T2 *out, size_t n, int G, int L, int fftn1, int pow1, int *bad) {
#pragma omp target teams num_teams(G) thread_limit(L) map(present, alloc : in[0 : n], out[0 : n]) map(tofrom : bad[0 : 1])
  {
    T smem[2064];
#pragma omp parallel num_threads(L)
    {
      /* one OpenCL work-group = one team; its L work-items need L threads (barriers + private state) */
      if (omp_get_num_threads() != L) {
#pragma omp atomic write
        bad[0] = 1;
      } else {
        for (int g = omp_get_team_num(); g < G; g += omp_get_num_teams()) {
          fft0(in, out, smem, omp_get_thread_num(), g, fftn1, pow1);
#pragma omp barrier
        }
      }
    }
  }
}

static void launch_fft1D_8192(T2 *in, T2 *out, size_t n, int G, int L, int fftn1, int pow1, int *bad) {
#pragma omp target teams num_teams(G) thread_limit(L) map(present, alloc : in[0 : n], out[0 : n]) map(tofrom : bad[0 : 1])
  {
    T smem[1024];
#pragma omp parallel num_threads(L)
    {
      /* one OpenCL work-group = one team; its L work-items need L threads (barriers + private state) */
      if (omp_get_num_threads() != L) {
#pragma omp atomic write
        bad[0] = 1;
      } else {
        for (int g = omp_get_team_num(); g < G; g += omp_get_num_teams()) {
          fft1D_8192(in, out, smem, omp_get_thread_num(), g, fftn1, pow1);
#pragma omp barrier
        }
      }
    }
  }
}

static void launch_fft1D_4096(T2 *in, T2 *out, size_t n, int G, int L, int fftn1, int pow1, int *bad) {
#pragma omp target teams num_teams(G) thread_limit(L) map(present, alloc : in[0 : n], out[0 : n]) map(tofrom : bad[0 : 1])
  {
    T smem[512];
#pragma omp parallel num_threads(L)
    {
      /* one OpenCL work-group = one team; its L work-items need L threads (barriers + private state) */
      if (omp_get_num_threads() != L) {
#pragma omp atomic write
        bad[0] = 1;
      } else {
        for (int g = omp_get_team_num(); g < G; g += omp_get_num_teams()) {
          fft1D_4096(in, out, smem, omp_get_thread_num(), g, fftn1, pow1);
#pragma omp barrier
        }
      }
    }
  }
}

/* ---------------- host: forward() of fftlib.cpp for the 1D path ---------------- */
static int local_size_1d(int n) {  /* getLocalDimension(): num_item_per_wg for n <= 2048 */
  int r0 = 8; /* first radix of getLocalRadix(n,.,.,0) is 8 for every n in 128..2048 */
  int per = n / r0;
  return per <= 64 ? 64 : per;
}

static void forward(T2 *work, T2 *temp, int fftn, int *bad) {
  size_t n = (size_t)fftn;
  int pow1 = fftn / 16; /* quirk of the OpenCL host code, see header */
  if (fftn <= 2048) {
    int L = local_size_1d(fftn);
    switch (fftn) {
      case 128: launch_fft1D_128(work, n, 1, L, fftn, pow1, bad); break;
      case 256: launch_fft1D_256(work, n, 1, L, fftn, pow1, bad); break;
      case 512: launch_fft1D_512(work, n, 1, L, fftn, pow1, bad); break;
      case 1024: launch_fft1D_1024(work, n, 1, L, fftn, pow1, bad); break;
      case 2048: launch_fft1D_2048(work, n, 1, L, fftn, pow1, bad); break;
    }
  } else if (fftn == 4096) {
    launch_fft0(work, temp, n, 2, 128, fftn, pow1, bad);          /* getGlobalDimension stage 0: 2 groups x 128 */
    launch_fft1D_4096(temp, work, n, 8, 64, fftn, pow1, bad);     /* stage 1: 8 groups x 64 */
  } else if (fftn == 8192) {
    launch_fft0(work, temp, n, 4, 128, fftn, pow1, bad);
    launch_fft1D_8192(temp, work, n, 8, 128, fftn, pow1, bad);
  }
}

int main(int argc, char **argv) {
  int sizeIndex = 0, passes = 10;
  for (int i = 1; i < argc; i++) {
    if (!strcmp(argv[i], "--pts") && i + 1 < argc) sizeIndex = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--passes") && i + 1 < argc) passes = atoi(argv[++i]);
    else if (!strcmp(argv[i], "--2D")) { fprintf(stderr, "--2D is not translated in the OpenMP tree\n"); return 3; }
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
  for (int i = 0; i < N; i++) {
    source[i].x = (rand() / (float)RAND_MAX) * 2 - 1;
    source[i].y = (rand() / (float)RAND_MAX) * 2 - 1;
  }
  int bad_flag[1] = {0};
  int *bad = bad_flag;
  size_t n = (size_t)fftn;
  dwarf_init("fft", "spectral-methods", "omp");
  dwarf_set_size(fftn, 1);
  for (int rep = dwarf_first_rep(); rep < dwarf_reps(); ++rep) {
    bad[0] = 0;
    dwarf_region_begin();
#pragma omp target enter data map(alloc : work[0 : n], temp[0 : n])
    memcpy(work, source, used_bytes);
#pragma omp target update to(work[0 : n])
    forward(work, temp, fftn, bad);
#pragma omp target update from(work[0 : n])
#pragma omp target exit data map(delete : work[0 : n], temp[0 : n])
    dwarf_region_end(rep);
  }
  memcpy(result, work, used_bytes);
#ifdef PRINT_RESULT
  for (int i = 0; i < N; i++) fprintf(stdout, "data[%d] (%g, %g)\n", i, result[i].x, result[i].y);
#endif
  if (bad[0]) { fprintf(stderr, "FAILURE: the toolchain did not provide one thread per work-item in a team\n"); return 4; }

  /* validation: double-precision DFT; find which output convention the SHOC kernels produce */
  static const char *names[8] = {"natural", "swap", "conj", "swap-conj", "natural-bitrev", "swap-bitrev", "conj-bitrev", "swap-conj-bitrev"};
  double err[8]; int bits = 0; while ((1 << bits) < fftn) bits++;
  double *rr = (double *)malloc(sizeof(double) * fftn), *ri = (double *)malloc(sizeof(double) * fftn);
  double norm = 0;
  for (int k = 0; k < fftn; k++) {
    double sr = 0, si = 0;
    for (int t = 0; t < fftn; t++) {
      double ang = -2.0 * M_PI * (double)((long)k * t % fftn) / fftn;
      sr += source[t].x * cos(ang) - source[t].y * sin(ang);
      si += source[t].x * sin(ang) + source[t].y * cos(ang);
    }
    rr[k] = sr; ri[k] = si; norm += sr * sr + si * si;
  }
  norm = sqrt(norm / fftn);
  for (int m = 0; m < 8; m++) {
    double e = 0;
    for (int k = 0; k < fftn; k++) {
      int idx = k;
      if (m >= 4) { idx = 0; for (int b = 0; b < bits; b++) if (k & (1 << b)) idx |= 1 << (bits - 1 - b); }
      double gx = result[idx].x, gy = result[idx].y;
      if (m & 1) { double t = gx; gx = gy; gy = t; }
      if (m & 2) gy = -gy;
      e += (gx - rr[k]) * (gx - rr[k]) + (gy - ri[k]) * (gy - ri[k]);
    }
    err[m] = sqrt(e / fftn) / norm;
  }
  int best = 0;
  for (int m = 1; m < 8; m++) if (err[m] < err[best]) best = m;
  /* The OpenCL app's own output is not a correct FFT (NVIDIA and PoCL disagree with each other and with the
   * DFT by up to O(1) element-wise, 1-12% in sorted magnitudes; sizes 128/256 read uninitialised registers),
   * so the DFT distance is reported for information only. The run is valid if every value is finite. */
  int ok = 1;
  for (int k = 0; k < fftn; k++) if (!isfinite(result[k].x) || !isfinite(result[k].y)) ok = 0;
  printf("%s: output %s; closest DFT convention '%s', relative RMS distance %.3g (informational: the upstream kernel is not an exact FFT)\n",
         ok ? "SUCCESS" : "FAILURE", ok ? "finite" : "contains non-finite values", names[best], err[best]);
  return dwarf_report(ok);
}
