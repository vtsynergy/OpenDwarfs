/* CUDA translation of the 1-D part of fft.cl (SHOC). Mechanical: T2 = float2, __local -> __shared__, barrier -> __syncthreads,
 * mad24/mul24 -> plain integer arithmetic. The -D fftn1/pow1 build options of fft0 are template parameters. FFT_2D and TWIDDLE are
 * undefined, as in the OpenCL 1-D build; the 2-D kernels are not translated. */
#ifndef FFT_NATIVE_COS
#define FFT_NATIVE_COS __cosf
#define FFT_NATIVE_SIN __sinf
#endif
#define mad24(a,b,c) ((a)*(b)+(c))
#define mul24(a,b) ((a)*(b))

#define T float
#define T2 float2
#define MK2(...) make_float2(__VA_ARGS__)

#define numIter Radix1/Radix2
#define STORE_REV 1
#define LOAD_REV 2
#define NO_REV 0
#define STORE_S 3
#define LOAD_S 4


#define COS_PI_8  0.923879533f
#define SIN_PI_8  0.382683432f
#undef M_SQRT1_2
#define M_SQRT1_2 7.0710678118654752440E-1
#undef M_PI
#define M_PI      3.141592653589793
#define exp_1_8   make_float2((T)(M_SQRT1_2), (T)(M_SQRT1_2))//requires post-multiply by 1/sqrt(2)
#define exp_1_4   make_float2((T)(0), (T)(1))
#define exp_3_8   make_float2((T)(-M_SQRT1_2), (T)(M_SQRT1_2))//requires post-multiply by 1/sqrt(2)
#define exp_1_16  make_float2((T)(COS_PI_8), (T)(SIN_PI_8))
#define exp_3_16  make_float2((T)(M_SQRT1_2), (T)(M_SQRT1_2))
#define exp_5_16  make_float2((T)(SIN_PI_8), (T)(COS_PI_8))
#define exp_7_16  make_float2((T)(-M_SQRT1_2), (T)(M_SQRT1_2))
#define exp_9_16  make_float2((T)(-COS_PI_8), (T)(-SIN_PI_8))
__device__ static __forceinline__ T2 exp_i( T phi ) {
	return make_float2( FFT_NATIVE_COS(phi), FFT_NATIVE_SIN(phi) );
}
__device__ static __forceinline__ T2 cmplx_mul( T2 a, T2 b ) { return make_float2( a.x*b.x-a.y*b.y, a.x*b.y+a.y*b.x ); }
__device__ static __forceinline__ T2 cm_fl_mul( T2 a, T  b ) { return make_float2( b*a.x, b*a.y ); }
__device__ static __forceinline__ T2 cmplx_add( T2 a, T2 b ) { return make_float2( a.x + b.x, a.y + b.y ); }
__device__ static __forceinline__ T2 cmplx_sub( T2 a, T2 b ) { return make_float2( a.x - b.x, a.y - b.y ); }
__device__ static __forceinline__ T2 cmplx_swap( T2 a) {return make_float2( a.y, a.x);}

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




__device__ static __forceinline__ void fftKernel(T2 *a, int n)
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



__device__ static __forceinline__ void globalLoads8(T2 *data, T2 *in, int S){
	for( int i = 0; i < 8; i++ )
		data[i] = in[i*S];
}

__device__ static __forceinline__ void globalLoads16(T2 *data, T2 *in, int S){
	for( int i = 0; i < 16; i++ )
		data[i] = in[i*S];
}


__device__ static __forceinline__ void globalStores8_r(T2 *data, T2 *out, int S){
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
#ifdef FFT_2D
		out[i*S] = data[reversed[i]];
#else
	out[i*S] = cmplx_swap(data[reversed[i]]);
#endif
}


__device__ static __forceinline__ void globalStores8(T2 *data, T2 *out, int S){
	for( int i = 0; i < 8; i++ )
#ifdef FFT_2D
		out[i*S] = data[i];
#else
	out[i*S] = cmplx_swap(data[i]);
#endif
}
__device__ static __forceinline__ void globalStores16(T2 *data, T2 *out, int S){
	for( int i = 0; i < 16; i++ )
		out[i*S] = data[i];

}

__device__ static __forceinline__ void storex4( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 4; i++ )
		x[i*sx] = a[i].x;
}

__device__ static __forceinline__ void storey4( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 4; i++ )
		x[i*sx] = a[i].y;
}

__device__ static __forceinline__ void storex16( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 16; i++ )
		x[i*sx] = a[i].x;
}

__device__ static __forceinline__ void storey16( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 16; i++ )
		x[i*sx] = a[i].y;
}


__device__ static __forceinline__ void storex8( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 8; i++ )
		x[i*sx] = a[i].x;
}

__device__ static __forceinline__ void storey8( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 8; i++ )
		x[i*sx] = a[i].y;
}

__device__ static __forceinline__ void storex8_r( T2 *a, T *x, int sx ) {
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
		x[i*sx] = a[reversed[i]].x;
}

__device__ static __forceinline__ void storey8_r( T2 *a, T *x, int sx ) {
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
		x[i*sx] = a[reversed[i]].y;
}

__device__ static __forceinline__ void loadx16( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 16; i++ )
		a[i].x = x[i*sx];
}

__device__ static __forceinline__ void loady16( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 16; i++ )
		a[i].y = x[i*sx];
}

__device__ static __forceinline__ void loadx8( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 8; i++ )
		a[i].x = x[i*sx];
}

__device__ static __forceinline__ void loady8( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 8; i++ )
		a[i].y = x[i*sx];
}

__device__ static __forceinline__ void loadx4( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 4; i++ )
		a[i].x = x[i*sx];
}

__device__ static __forceinline__ void loady4( T2 *a, T *x, int sx ) {
	for( int i = 0; i < 4; i++ )
		a[i].y = x[i*sx];
}

__device__ static __forceinline__ void loadx8_r( T2 *a, T *x, int sx ) {
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
		a[i].x = x[reversed[i]*sx];
}

__device__ static __forceinline__ void loady8_r( T2 *a, T *x, int sx ) {
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
		a[i].y = x[reversed[i]*sx];
}

__device__ static __forceinline__ void globalLoads(T2 *data, T2 *in, int S, int n){
	for( int i = 0; i < n; i++ )
		data[i] = in[i*S];
}

__device__ static __forceinline__ void globalStores(T2 *data, T2 *out, int S, int n, bool rev){
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

__device__ static __forceinline__ void globalStores_s(T2 *data, T2 *out, int S, int n, int offset){
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

__device__ static __forceinline__ void storex( T2 *a, T *x, int sx, int n, bool rev, int offset){
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

__device__ static __forceinline__ void storey( T2 *a, T *x, int sx, int n, bool rev, int offset){
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

__device__ static __forceinline__ void loadx( T2 *a, T *x, int sx, int n, bool rev, int offset ) {
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
__device__ static __forceinline__ void loady( T2 *a, T *x, int sx, int n, bool rev,int offset ) {
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

__device__ static __forceinline__ void trans(T2 *a, T *s, int ds, T *l, int dl, int n, int rev){
	storex(a, s, ds, n, rev==STORE_REV, 0);
	__syncthreads();
	loadx(a, l, dl, n, rev==LOAD_REV, 0);
	__syncthreads();
	storey(a, s, ds, n, rev==STORE_REV, 0);
	__syncthreads();
	loady(a, l, dl, n, rev==LOAD_REV, 0);
	__syncthreads();
}

__device__ static __forceinline__ void trans_s(T2 *a, T *s, int ds, T *l, int dl, int n, int offset_l, int offset_s){
	storex(a, s, ds, n/2, NO_REV, 0);
	storex(a, s+offset_s, ds, n/2, NO_REV, n/2);
	__syncthreads();
	loadx(a, l, dl, n/2, NO_REV, 0);
	loadx(a, l+offset_l, dl, n/2, NO_REV, n/2);
	__syncthreads();
	storey(a, s, ds, n/2, NO_REV, 0);
	storey(a, s+offset_s, ds, n/2, NO_REV, n/2);
	__syncthreads();
	loady(a, l, dl, n/2, NO_REV, 0);
	loady(a, l+offset_l, dl, n/2, NO_REV, n/2);
	__syncthreads();
}



__device__ static __forceinline__ void twiddle_s(T2 *a, int j, int N, int n, int k )                                              \
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
	storex8( a, s, ds );  __syncthreads();  \
	loadx8 ( a, l, dl );  __syncthreads();  \
	storey8( a, s, ds );  __syncthreads();  \
	loady8 ( a, l, dl );  __syncthreads();  \
}

#define transpose16( a, s, ds, l, dl )                              \
{                                                                       \
	storex16( a, s, ds );  __syncthreads();  \
	loadx16 ( a, l, dl );  __syncthreads();  \
	storey16( a, s, ds );  __syncthreads();  \
	loady16 ( a, l, dl );  __syncthreads();  \
}


__device__ static __forceinline__ void twiddle16(T2 *a, int tt, int n)
{
	T c = (T) tt;
	for( int t = 1; t < 16; t++ ){
		(a)[t] = cmplx_mul( (a)[t],exp_i((2.0f*M_PI*t/n)*c));
	}
}
__device__ static __forceinline__ void twiddle8(T2 *a, int tt, int n)
{
	T c = (T) tt;
	for( int t = 1; t < 8; t++ ){
		(a)[t] = cmplx_mul( (a)[t],exp_i((2.0f*M_PI*t/n)*c));
	}
}
__device__ static __forceinline__ void  twiddle4(T2 *a, int tt, int n )                                 {
	T c = (T) tt;
	for( int t = 1; t < 4; t++ ){
		(a)[t] = cmplx_mul( (a)[t],exp_i((2.0f*M_PI*t/n)*c));
	}
}

__global__ void fft1D_512(T2 *in)
{
	__shared__ T smem[576];
	int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	T *Local_store, *Local_load;
	T2 a[8];
	int local_id = threadIdx.x;
	int group_id = blockIdx.x;
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

__global__ void fft1D_1024(T2 *in)
{
	__shared__ T smem[1088];
	int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	T *Local_store, *Local_load;
	T2 a[8];
	int local_id = threadIdx.x;
	int group_id = blockIdx.x;
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

__global__ void fft1D_2048(T2 *in)
{
	__shared__ T smem[2112];
	int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	T *Local_store, *Local_load;
	T2 a[8];
	int local_id = threadIdx.x;
	int group_id = blockIdx.x;
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

__global__ void fft1D_256(T2 *in)
{
	__shared__ T smem[640];
	int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	T *Local_store, *Local_load;
	T2 a[8];
	int local_id = threadIdx.x;
	int group_id = blockIdx.x;
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
	/* [FIX-OCL-25] only the jj == 0 row of work-items loads data (OpenDwarfs runs one transform per group); the other rows
	 * still executed this store with `in` un-offset and uninitialised registers, overwriting the result of row 0 (fft1D_256) */
	if( jj == 0 ) globalStores8_r(a, in, 32);
}

__global__ void fft1D_128(T2 *in)
{
	__shared__ T smem[640];
	int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	T *Local_store, *Local_load;
	T2 a[8];
	int local_id = threadIdx.x;
	int group_id = blockIdx.x;
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
	/* [FIX-OCL-25] as in fft1D_256: rows jj != 0 must not store (their `in` is un-offset and their data uninitialised) */
	if( jj == 0 ) globalStores8_r(a, in, 16);
}

template <int fftn1, int pow1>
__global__ void fft0(T2 *in, T2 *out)
{
	__shared__ T smem[2064];
	int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	T *Local_store, *Local_load;
	T2 a[16];
	int local_id = threadIdx.x;
	int group_id = blockIdx.x;
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
	trans_s(a, Local_store, 16, Local_load, 129, 16, 0, -127);
	Offset_out += tid;
	out += Offset_out;
	globalStores16(a, out, 128);
	//globalStores(a, out, 128, 16);
}

__global__ void fft1D_8192(T2 *in, T2 *out)
{
	__shared__ T smem[1024];
	int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	T *Local_store, *Local_load;
	T2 a[8];
	int local_id = threadIdx.x;
	int group_id = blockIdx.x;
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

__global__ void fft1D_4096(T2 *in, T2 *out)
{
	__shared__ T smem[512];
	int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	T *Local_store, *Local_load;
	T2 a[8];
	int local_id = threadIdx.x;
	int group_id = blockIdx.x;
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

