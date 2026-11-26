// Generated from opencl/spectral-methods/fft/src/opencl/fft/fft.cl (1-D kernels; the 2-D kernels fft1/fft2 are not
// translated). Bodies are verbatim apart from: float2 -> struct T2, "(T2)(a,b)" literals, address-space qualifiers,
// native_cos/sin -> sycl::native, barriers -> group_barrier, __local arrays -> local_accessor, and the build options
// -D fftn1/-D pow1 -> template parameters of fft0. TWIDDLE and FFT_2D are undefined, as in the 1-D OpenCL build.
// M_PI and M_SQRT1_2 stay double constants, as in the OpenCL source.
#ifndef FFT_KERNELS_HPP
#define FFT_KERNELS_HPP
#include <sycl/sycl.hpp>
namespace fftk {
struct alignas(8) T2 { float x, y; };
using T = float;
inline int mad24(int a, int b, int c) { return a*b + c; }   // OpenCL built-ins, all arguments are small ints here
inline int mul24(int a, int b) { return a*b; }

#define numIter Radix1/Radix2
#define STORE_REV 1
#define LOAD_REV 2
#define NO_REV 0
#define STORE_S 3
#define LOAD_S 4


#define COS_PI_8  0.923879533f
#define SIN_PI_8  0.382683432f
#define K_M_SQRT1_2 7.0710678118654752440E-1
#define K_M_PI 3.141592653589793
#define exp_1_8   T2{(float)(K_M_SQRT1_2), (float)(K_M_SQRT1_2)}//requires post-multiply by 1/sqrt(2)
#define exp_1_4   T2{(float)(0), (float)(1)}
#define exp_3_8   T2{(float)(-K_M_SQRT1_2), (float)(K_M_SQRT1_2)}//requires post-multiply by 1/sqrt(2)
#define exp_1_16  T2{(float)(COS_PI_8), (float)(SIN_PI_8)}
#define exp_3_16  T2{(float)(K_M_SQRT1_2), (float)(K_M_SQRT1_2)}
#define exp_5_16  T2{(float)(SIN_PI_8), (float)(COS_PI_8)}
#define exp_7_16  T2{(float)(-K_M_SQRT1_2), (float)(K_M_SQRT1_2)}
#define exp_9_16  T2{(float)(-COS_PI_8), (float)(-SIN_PI_8)}
inline T2 exp_i( T phi ) {
	return T2{(float)(sycl::native::cos(phi)), (float)(sycl::native::sin(phi))};
}
inline T2 cmplx_mul( T2 a, T2 b ) { return T2{(float)(a.x*b.x-a.y*b.y), (float)(a.x*b.y+a.y*b.x)}; }
inline T2 cm_fl_mul( T2 a, T  b ) { return T2{(float)(b*a.x), (float)(b*a.y)}; }
inline T2 cmplx_add( T2 a, T2 b ) { return T2{(float)(a.x + b.x), (float)(a.y + b.y)}; }
inline T2 cmplx_sub( T2 a, T2 b ) { return T2{(float)(a.x - b.x), (float)(a.y - b.y)}; }
inline T2 cmplx_swap( T2 a) {return T2{(float)(a.y), (float)(a.x)};}

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




inline void fftKernel(T2 *a, int n)
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



inline void globalLoads8(T2 *data, T2 *in, int S){
	for( int i = 0; i < 8; i++ )
		data[i] = in[i*S];
}

inline void globalLoads16(T2 *data, T2 *in, int S){
	for( int i = 0; i < 16; i++ )
		data[i] = in[i*S];
}


inline void globalStores8_r(T2 *data, T2 *out, int S){
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
#ifdef FFT_2D
		out[i*S] = data[reversed[i]];
#else
	out[i*S] = cmplx_swap(data[reversed[i]]);
#endif
}


inline void globalStores8(T2 *data, T2 *out, int S){
	for( int i = 0; i < 8; i++ )
#ifdef FFT_2D
		out[i*S] = data[i];
#else
	out[i*S] = cmplx_swap(data[i]);
#endif
}
inline void globalStores16(T2 *data, T2 *out, int S){
	for( int i = 0; i < 16; i++ )
		out[i*S] = data[i];

}

inline void storex4( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 4; i++ )
		x[i*sx] = a[i].x;
}

inline void storey4( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 4; i++ )
		x[i*sx] = a[i].y;
}

inline void storex16( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 16; i++ )
		x[i*sx] = a[i].x;
}

inline void storey16( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 16; i++ )
		x[i*sx] = a[i].y;
}


inline void storex8( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 8; i++ )
		x[i*sx] = a[i].x;
}

inline void storey8( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 8; i++ )
		x[i*sx] = a[i].y;
}

inline void storex8_r( T2 *a, T *x, int sx ) {
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
		x[i*sx] = a[reversed[i]].x;
}

inline void storey8_r( T2 *a, T *x, int sx ) {
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
		x[i*sx] = a[reversed[i]].y;
}

inline void loadx16( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 16; i++ )
		a[i].x = x[i*sx];
}

inline void loady16( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 16; i++ )
		a[i].y = x[i*sx];
}

inline void loadx8( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 8; i++ )
		a[i].x = x[i*sx];
}

inline void loady8( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 8; i++ )
		a[i].y = x[i*sx];
}

inline void loadx4( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 4; i++ )
		a[i].x = x[i*sx];
}

inline void loady4( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 4; i++ )
		a[i].y = x[i*sx];
}

inline void loadx8_r( T2 *a, T *x, int sx ) {
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
		a[i].x = x[reversed[i]*sx];
}

inline void loady8_r( T2 *a, T *x, int sx ) {
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
		a[i].y = x[reversed[i]*sx];
}

inline void globalLoads(T2 *data, T2 *in, int S, int n){
	for( int i = 0; i < n; i++ )
		data[i] = in[i*S];
}

inline void globalStores(T2 *data, T2 *out, int S, int n, bool rev){
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

inline void globalStores_s(T2 *data, T2 *out, int S, int n, int offset){
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

inline void storex( T2 *a, T *x, int sx, int n, bool rev, int offset){
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

inline void storey( T2 *a, T *x, int sx, int n, bool rev, int offset){
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

inline void loadx( T2 *a, T *x, int sx, int n, bool rev, int offset ) {
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
inline void loady( T2 *a, T *x, int sx, int n, bool rev,int offset ) {
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

inline void trans(const sycl::nd_item<1> &it, T2 *a, T *s, int ds, T *l, int dl, int n, int rev){
	storex(a, s, ds, n, rev==STORE_REV, 0);
	sycl::group_barrier(it.get_group());
	loadx(a, l, dl, n, rev==LOAD_REV, 0);
	sycl::group_barrier(it.get_group());
	storey(a, s, ds, n, rev==STORE_REV, 0);
	sycl::group_barrier(it.get_group());
	loady(a, l, dl, n, rev==LOAD_REV, 0);
	sycl::group_barrier(it.get_group());
}

inline void trans_s(const sycl::nd_item<1> &it, T2 *a, T *s, int ds, T *l, int dl, int n, int offset_l, int offset_s){
	storex(a, s, ds, n/2, NO_REV, 0);
	storex(a, s+offset_s, ds, n/2, NO_REV, n/2);
	sycl::group_barrier(it.get_group());
	loadx(a, l, dl, n/2, NO_REV, 0);
	loadx(a, l+offset_l, dl, n/2, NO_REV, n/2);
	sycl::group_barrier(it.get_group());
	storey(a, s, ds, n/2, NO_REV, 0);
	storey(a, s+offset_s, ds, n/2, NO_REV, n/2);
	sycl::group_barrier(it.get_group());
	loady(a, l, dl, n/2, NO_REV, 0);
	loady(a, l+offset_l, dl, n/2, NO_REV, n/2);
	sycl::group_barrier(it.get_group());
}



inline void twiddle_s(T2 *a, int j, int N, int n, int k )                                              \
{                                                                       \
	for( int i = 0; i < n/2; i++ ){			\
		a[i] = cmplx_mul( a[i],exp_i(((2.0f*K_M_PI/N)*j)*(k+16*i))); \
	} 							\
	for( int i = n/2; i < n; i++ ){			\
		a[i] = cmplx_mul( a[i],exp_i(((2.0f*K_M_PI/N)*j)*(k+16*(i-n/2)+1))); \
	}                                                                   \
	\
}

#define transpose( a, s, ds, l, dl )                              \
{                                                                       \
	storex8( a, s, ds );  sycl::group_barrier(it.get_group());  \
	loadx8 ( a, l, dl );  sycl::group_barrier(it.get_group());  \
	storey8( a, s, ds );  sycl::group_barrier(it.get_group());  \
	loady8 ( a, l, dl );  sycl::group_barrier(it.get_group());  \
}

#define transpose16( a, s, ds, l, dl )                              \
{                                                                       \
	storex16( a, s, ds );  sycl::group_barrier(it.get_group());  \
	loadx16 ( a, l, dl );  sycl::group_barrier(it.get_group());  \
	storey16( a, s, ds );  sycl::group_barrier(it.get_group());  \
	loady16 ( a, l, dl );  sycl::group_barrier(it.get_group());  \
}


inline void twiddle16(T2 *a, int tt, int n)
{
	T c = (T) tt;
	for( int t = 1; t < 16; t++ ){
		(a)[t] = cmplx_mul( (a)[t],exp_i((2.0f*K_M_PI*t/n)*c));
	}
}
inline void twiddle8(T2 *a, int tt, int n)
{
	T c = (T) tt;
	for( int t = 1; t < 8; t++ ){
		(a)[t] = cmplx_mul( (a)[t],exp_i((2.0f*K_M_PI*t/n)*c));
	}
}
inline void  twiddle4(T2 *a, int tt, int n )                                 {
	T c = (T) tt;
	for( int t = 1; t < 4; t++ ){
		(a)[t] = cmplx_mul( (a)[t],exp_i((2.0f*K_M_PI*t/n)*c));
	}
}


static void fft1D_512(sycl::queue &q, size_t global, size_t local, T2 *in_)
{
	q.submit([&](sycl::handler &h) {
		sycl::local_accessor<float, 1> smem_acc(sycl::range<1>(576), h);
		h.parallel_for(sycl::nd_range<1>(global, local), [=](sycl::nd_item<1> it) {
		float *smem = smem_acc.get_multi_ptr<sycl::access::decorated::no>().get();
		T2 *in = in_;
		int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
		int ii, jj, offset;
		T2 w;
		T *Local_store, *Local_load;
		T2 a[8];
		int local_id = (int)it.get_local_id(0);
		int group_id = (int)it.get_group(0);
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
		});
	}).wait_and_throw();
}

static void fft1D_1024(sycl::queue &q, size_t global, size_t local, T2 *in_)
{
	q.submit([&](sycl::handler &h) {
		sycl::local_accessor<float, 1> smem_acc(sycl::range<1>(1088), h);
		h.parallel_for(sycl::nd_range<1>(global, local), [=](sycl::nd_item<1> it) {
		float *smem = smem_acc.get_multi_ptr<sycl::access::decorated::no>().get();
		T2 *in = in_;
		int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
		int ii, jj, offset;
		T2 w;
		T *Local_store, *Local_load;
		T2 a[8];
		int local_id = (int)it.get_local_id(0);
		int group_id = (int)it.get_group(0);
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
		trans_s(it, a,Local_store,136,Local_load,32,8,-112,0);
		FFT4S(a+0);
		FFT4S(a+4);
		twiddle4(a,ii>>6,16);
		twiddle4(a+4,(128 + ii) >>6,16);
		Local_store = smem + ii;
		j = ii >> 6;
		i = ii & 63;
		Local_load = smem + mad24(j, 256, i);
		trans_s(it, a,Local_store,256,Local_load,64,8,256,-896);
		FFT4S(a);
		FFT4S(a+4);
		globalStores8(a,in,128);
		});
	}).wait_and_throw();
}

static void fft1D_2048(sycl::queue &q, size_t global, size_t local, T2 *in_)
{
	q.submit([&](sycl::handler &h) {
		sycl::local_accessor<float, 1> smem_acc(sycl::range<1>(2112), h);
		h.parallel_for(sycl::nd_range<1>(global, local), [=](sycl::nd_item<1> it) {
		float *smem = smem_acc.get_multi_ptr<sycl::access::decorated::no>().get();
		T2 *in = in_;
		int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
		int ii, jj, offset;
		T2 w;
		T *Local_store, *Local_load;
		T2 a[8];
		int local_id = (int)it.get_local_id(0);
		int group_id = (int)it.get_group(0);
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
		trans_s(it, a, Local_store, 256, Local_load, 64, 8, 768, 0);
		FFT4S(a);
		FFT4S(a+4);
		globalStores8_r(a, in, 256);
		});
	}).wait_and_throw();
}

static void fft1D_256(sycl::queue &q, size_t global, size_t local, T2 *in_)
{
	q.submit([&](sycl::handler &h) {
		sycl::local_accessor<float, 1> smem_acc(sycl::range<1>(640), h);
		h.parallel_for(sycl::nd_range<1>(global, local), [=](sycl::nd_item<1> it) {
		float *smem = smem_acc.get_multi_ptr<sycl::access::decorated::no>().get();
		T2 *in = in_;
		int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
		int ii, jj, offset;
		T2 w;
		T *Local_store, *Local_load;
		T2 a[8];
		int local_id = (int)it.get_local_id(0);
		int group_id = (int)it.get_group(0);
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
		trans_s(it, a, Local_store, 40, Local_load, 8, 8,128, 0);
		FFT4S(a);
		FFT4S(a+4);
	
		if( jj == 0 ) globalStores8_r(a, in, 32);
		});
	}).wait_and_throw();
}

static void fft1D_128(sycl::queue &q, size_t global, size_t local, T2 *in_)
{
	q.submit([&](sycl::handler &h) {
		sycl::local_accessor<float, 1> smem_acc(sycl::range<1>(640), h);
		h.parallel_for(sycl::nd_range<1>(global, local), [=](sycl::nd_item<1> it) {
		float *smem = smem_acc.get_multi_ptr<sycl::access::decorated::no>().get();
		T2 *in = in_;
		int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
		int ii, jj, offset;
		T2 w;
		T *Local_store, *Local_load;
		T2 a[8];
		int local_id = (int)it.get_local_id(0);
		int group_id = (int)it.get_group(0);
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
		trans_s(it, a, Local_store, 18, Local_load, 4, 8, -14, 0);
		FFT4S(a);
		FFT4S(a+4);
		twiddle4(a, ii>>3, 16);
		twiddle4(a+4, (16 + ii) >>3, 16);
		Local_store = smem + mad24(jj, 160, ii);
		j = ii >> 3;
		i = ii & 7;
		i = mad24(jj, 160, i);
		Local_load = smem + mad24(j, 40, i);
		trans_s(it, a, Local_store, 40, Local_load, 8, 8, 48, -144);
		FFT4S(a);
		FFT4S(a+4);
	
		if( jj == 0 ) globalStores8_r(a, in, 16);
		});
	}).wait_and_throw();
}

template <int fftn1, int pow1>
static void fft0(sycl::queue &q, size_t global, size_t local, T2 *in_, T2 *out_)
{
	q.submit([&](sycl::handler &h) {
		sycl::local_accessor<float, 1> smem_acc(sycl::range<1>(2064), h);
		h.parallel_for(sycl::nd_range<1>(global, local), [=](sycl::nd_item<1> it) {
		float *smem = smem_acc.get_multi_ptr<sycl::access::decorated::no>().get();
		T2 *in = in_; T2 *out = out_;
		int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
		int ii, jj, offset;
		T2 w;
		T *Local_store, *Local_load;
		T2 a[16];
		int local_id = (int)it.get_local_id(0);
		int group_id = (int)it.get_group(0);
		b = group_id & (fftn1/2048-1);
		d = group_id >> (pow1-11);
		Offset_in = mul24(b, 16);
		tid = Offset_in;
		i = tid;
		j = 0;
		Offset_out = mad24(i, 128, j);
		Offset_in += (d << pow1);
		Offset_out += (d << pow1);
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
		trans_s(it, a, Local_store, 16, Local_load, 129, 16, 0, -127);
		Offset_out += tid;
		out += Offset_out;
		globalStores16(a, out, 128);
		//globalStores(a, out, 128, 16);
		});
	}).wait_and_throw();
}

static void fft1D_8192(sycl::queue &q, size_t global, size_t local, T2 *in_, T2 *out_)
{
	q.submit([&](sycl::handler &h) {
		sycl::local_accessor<float, 1> smem_acc(sycl::range<1>(1024), h);
		h.parallel_for(sycl::nd_range<1>(global, local), [=](sycl::nd_item<1> it) {
		float *smem = smem_acc.get_multi_ptr<sycl::access::decorated::no>().get();
		T2 *in = in_; T2 *out = out_;
		int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
		int ii, jj, offset;
		T2 w;
		T *Local_store, *Local_load;
		T2 a[8];
		int local_id = (int)it.get_local_id(0);
		int group_id = (int)it.get_group(0);
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
		});
	}).wait_and_throw();
}

static void fft1D_4096(sycl::queue &q, size_t global, size_t local, T2 *in_, T2 *out_)
{
	q.submit([&](sycl::handler &h) {
		sycl::local_accessor<float, 1> smem_acc(sycl::range<1>(512), h);
		h.parallel_for(sycl::nd_range<1>(global, local), [=](sycl::nd_item<1> it) {
		float *smem = smem_acc.get_multi_ptr<sycl::access::decorated::no>().get();
		T2 *in = in_; T2 *out = out_;
		int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
		int ii, jj, offset;
		T2 w;
		T *Local_store, *Local_load;
		T2 a[8];
		int local_id = (int)it.get_local_id(0);
		int group_id = (int)it.get_group(0);
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
		});
	}).wait_and_throw();
}

#undef numIter
#undef STORE_REV
#undef LOAD_REV
#undef NO_REV
#undef STORE_S
#undef LOAD_S
#undef COS_PI_8
#undef SIN_PI_8
#undef K_M_SQRT1_2
#undef K_M_PI
#undef exp_1_8
#undef exp_1_4
#undef exp_3_8
#undef exp_1_16
#undef exp_3_16
#undef exp_5_16
#undef exp_7_16
#undef exp_9_16
#undef FFT2
#undef FFT2S
#undef FFT4
#undef FFT4S
#undef FFT8
#undef FFT16
#undef transpose
#undef transpose16
} // namespace fftk
#endif
