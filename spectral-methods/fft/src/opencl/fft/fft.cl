#define T float
#define T2 float2

#define numIter Radix1/Radix2
#define STORE_REV 1
#define LOAD_REV 2
#define NO_REV 0
#define STORE_S 3
#define LOAD_S 4


#define COS_PI_8  0.923879533f
#define SIN_PI_8  0.382683432f
#define M_SQRT1_2 7.0710678118654752440E-1
#define M_PI      3.141592653589793
#define exp_1_8   (T2)(  M_SQRT1_2, M_SQRT1_2 )//requires post-multiply by 1/sqrt(2)
#define exp_1_4   (T2)(  0, 1 )
#define exp_3_8   (T2)( -M_SQRT1_2, M_SQRT1_2 )//requires post-multiply by 1/sqrt(2)
#define NULL      0
#define exp_1_16  (T2)(  COS_PI_8, SIN_PI_8 )
#define exp_3_16  (T2)(  M_SQRT1_2, M_SQRT1_2 )
#define exp_5_16  (T2)( SIN_PI_8, COS_PI_8 )
#define exp_7_16  (T2)( -M_SQRT1_2, M_SQRT1_2 )
#define exp_9_16  (T2)( -COS_PI_8,  -SIN_PI_8 )
inline T2 exp_i( T phi ) {
	return (T2)( native_cos(phi), native_sin(phi) );
}
inline T2 cmplx_mul( T2 a, T2 b ) { return (T2)( a.x*b.x-a.y*b.y, a.x*b.y+a.y*b.x ); }
inline T2 cm_fl_mul( T2 a, T  b ) { return (T2)( b*a.x, b*a.y ); }
inline T2 cmplx_add( T2 a, T2 b ) { return (T2)( a.x + b.x, a.y + b.y ); }
inline T2 cmplx_sub( T2 a, T2 b ) { return (T2)( a.x - b.x, a.y - b.y ); }
inline T2 cmplx_swap( T2 a) {return (T2)( a.y, a.x);}

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



inline void globalLoads8(T2 *data, __global T2 *in, int S){
	for( int i = 0; i < 8; i++ )
		data[i] = in[i*S];
}

inline void globalLoads16(T2 *data, __global T2 *in, int S){
	for( int i = 0; i < 16; i++ )
		data[i] = in[i*S];
}


inline void globalStores8_r(T2 *data, __global T2 *out, int S){
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
#ifdef FFT_2D
		out[i*S] = data[reversed[i]];
#else
	out[i*S] = cmplx_swap(data[reversed[i]]);
#endif
}


inline void globalStores8(T2 *data, __global T2 *out, int S){
	for( int i = 0; i < 8; i++ )
#ifdef FFT_2D
		out[i*S] = data[i];
#else
	out[i*S] = cmplx_swap(data[i]);
#endif
}
inline void globalStores16(T2 *data, __global T2 *out, int S){
	for( int i = 0; i < 16; i++ )
		out[i*S] = data[i];

}

inline void storex4( T2 *a, __local T *x, int sx ) {
	for( int i = 0; i < 4; i++ )
		x[i*sx] = a[i].x;
}

inline void storey4( T2 *a, __local T *x, int sx ) {
	for( int i = 0; i < 4; i++ )
		x[i*sx] = a[i].y;
}

inline void storex16( T2 *a, __local T *x, int sx ) {
	for( int i = 0; i < 16; i++ )
		x[i*sx] = a[i].x;
}

inline void storey16( T2 *a, __local T *x, int sx ) {
	for( int i = 0; i < 16; i++ )
		x[i*sx] = a[i].y;
}


inline void storex8( T2 *a, __local T *x, int sx ) {
	for( int i = 0; i < 8; i++ )
		x[i*sx] = a[i].x;
}

inline void storey8( T2 *a, __local T *x, int sx ) {
	for( int i = 0; i < 8; i++ )
		x[i*sx] = a[i].y;
}

inline void storex8_r( T2 *a, __local T *x, int sx ) {
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
		x[i*sx] = a[reversed[i]].x;
}

inline void storey8_r( T2 *a, __local T *x, int sx ) {
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
		x[i*sx] = a[reversed[i]].y;
}

inline void loadx16( T2 *a, __local T *x, int sx ) {
	for( int i = 0; i < 16; i++ )
		a[i].x = x[i*sx];
}

inline void loady16( T2 *a, __local T *x, int sx ) {
	for( int i = 0; i < 16; i++ )
		a[i].y = x[i*sx];
}

inline void loadx8( T2 *a, __local T *x, int sx ) {
	for( int i = 0; i < 8; i++ )
		a[i].x = x[i*sx];
}

inline void loady8( T2 *a, __local T *x, int sx ) {
	for( int i = 0; i < 8; i++ )
		a[i].y = x[i*sx];
}

inline void loadx4( T2 *a, __local T *x, int sx ) {
	for( int i = 0; i < 4; i++ )
		a[i].x = x[i*sx];
}

inline void loady4( T2 *a, __local T *x, int sx ) {
	for( int i = 0; i < 4; i++ )
		a[i].y = x[i*sx];
}

inline void loadx8_r( T2 *a, __local T *x, int sx ) {
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
		a[i].x = x[reversed[i]*sx];
}

inline void loady8_r( T2 *a, __local T *x, int sx ) {
	int reversed[] = {0,4,1,5,2,6,3,7};
	for( int i = 0; i < 8; i++ )
		a[i].y = x[reversed[i]*sx];
}

inline void globalLoads(T2 *data, __global T2 *in, int S, int n){
	for( int i = 0; i < n; i++ )
		data[i] = in[i*S];
}

inline void globalStores(T2 *data, __global T2 *out, int S, int n, bool rev){
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

inline void globalStores_s(T2 *data, __global T2 *out, int S, int n, int offset){
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

inline void storex( T2 *a, __local T *x, int sx, int n, bool rev, int offset){
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

inline void storey( T2 *a, __local T *x, int sx, int n, bool rev, int offset){
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

inline void loadx( T2 *a, __local T *x, int sx, int n, bool rev, int offset ) {
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
inline void loady( T2 *a, __local T *x, int sx, int n, bool rev,int offset ) {
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

inline void trans(T2 *a, __local T *s, int ds, __local T *l, int dl, int n, int rev){
	storex(a, s, ds, n, rev==STORE_REV, 0);
	barrier(CLK_LOCAL_MEM_FENCE);
	loadx(a, l, dl, n, rev==LOAD_REV, 0);
	barrier(CLK_LOCAL_MEM_FENCE);
	storey(a, s, ds, n, rev==STORE_REV, 0);
	barrier(CLK_LOCAL_MEM_FENCE);
	loady(a, l, dl, n, rev==LOAD_REV, 0);
	barrier(CLK_LOCAL_MEM_FENCE);
}

inline void trans_s(T2 *a, __local T *s, int ds, __local T *l, int dl, int n, int offset_l, int offset_s){
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



inline void twiddle_s(T2 *a, int j, int N, int n, int k )                                              \
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


inline void twiddle16(T2 *a, int tt, int n)
{
	T c = (T) tt;
	for( int t = 1; t < 16; t++ ){
		(a)[t] = cmplx_mul( (a)[t],exp_i((2.0f*M_PI*t/n)*c));
	}
}
inline void twiddle8(T2 *a, int tt, int n)
{
	T c = (T) tt;
	for( int t = 1; t < 8; t++ ){
		(a)[t] = cmplx_mul( (a)[t],exp_i((2.0f*M_PI*t/n)*c));
	}
}
inline void  twiddle4(T2 *a, int tt, int n )                                 {
	T c = (T) tt;
	for( int t = 1; t < 4; t++ ){
		(a)[t] = cmplx_mul( (a)[t],exp_i((2.0f*M_PI*t/n)*c));
	}
}

#ifdef FFT_512
__kernel void fft1D_512(__global T2 *in)
{
	__local T smem[576];
	int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	__local T *Local_store, *Local_load;
	T2 a[8];
	int local_id = get_local_id( 0 );
	int group_id = get_group_id( 0 );
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
__kernel void fft1D_1024(__global T2 *in)
{
	__local T smem[1088];
	int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	__local T *Local_store, *Local_load;
	T2 a[8];
	int local_id = get_local_id( 0 );
	int group_id = get_group_id( 0 );
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
__kernel void fft1D_2048(__global T2 *in)
{
	__local T smem[2112];
	int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	__local T *Local_store, *Local_load;
	T2 a[8];
	int local_id = get_local_id( 0 );
	int group_id = get_group_id( 0 );
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
__kernel void fft1D_256(__global T2 *in)
{
	__local T smem[640];
	int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	__local T *Local_store, *Local_load;
	T2 a[8];
	int local_id = get_local_id( 0 );
	int group_id = get_group_id( 0 );
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
__kernel void fft1D_128(__global T2 *in)
{
	__local T smem[640];
	int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	__local T *Local_store, *Local_load;
	T2 a[8];
	int local_id = get_local_id( 0 );
	int group_id = get_group_id( 0 );
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

__kernel void fft0(__global T2 *in, __global T2 *out)
{
	__local T smem[2064];
	int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	__local T *Local_store, *Local_load;
	T2 a[16];
	int local_id = get_local_id( 0 );
	int group_id = get_group_id( 0 );
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

#ifdef FFT_8192
__kernel void fft1D_8192(__global T2 *in, __global T2 *out)
{
	__local T smem[1024];
	int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	__local T *Local_store, *Local_load;
	T2 a[8];
	int local_id = get_local_id( 0 );
	int group_id = get_group_id( 0 );
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
__kernel void fft1D_4096(__global T2 *in, __global T2 *out)
{
	__local T smem[512];
	int i, j, r, Offset_in, Offset_out,  tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	__local T *Local_store, *Local_load;
	T2 a[8];
	int local_id = get_local_id( 0 );
	int group_id = get_group_id( 0 );
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

#ifdef FFT_2D
__kernel void fft1(__global T2 *in, __global T2 *out)
{
	__local T smem[2048];
	int i, j, r, Offset_in, Offset_out, tid, b, d, k, l;
	int ii, jj, offset;
	T2 a[16];
	int local_id = get_local_id( 0 );
	int group_id = get_group_id( 0 );
	d = group_id >> (pow1+pow2-11);
	group_id = group_id & (fftn1*fftn2/2048 -1);
	Offset_in = mad24(group_id, 16, d << (pow1+pow2));
	tid = mul24(group_id, 16);
	i = tid >> pow1;
	j = tid & (fftn1-1);
	Offset_out = mad24(i, fftn1*128, j + (d << (pow1+pow2)));
	b = group_id;
	tid = local_id;
	i = tid & 15;
	j = tid >> 4;
	Offset_in += mad24(j, S1, i);
	in += Offset_in;
	globalLoads16(a, in, fftn1*fftn2/16);
	FFT16(a);
	twiddle16(a, j, 128);
	Offset_in = mad24(j, 256, i);
	transpose16(a, &smem[tid], 128, &smem[Offset_in], 16);
	barrier(CLK_LOCAL_MEM_FENCE);
	FFT8(a);
	FFT8(a+8);
#ifdef TWIDDLE
	l = ((b << 4) + i) >> (pow1);
	k = j << 1;
	twiddle_s(a, l, fftn2, 16, k);
#endif
	Offset_out += mad24(j, fftn1*2, i);
	out += Offset_out;
	globalStores_s(a, out, fftn1*16, 16, fftn1);
}
__kernel void fft2(__global T2 *in, __global T2 *out)
{
	__local T smem[256];
	int i, j, r, Offset_in, Offset_out, tid, b, d, k, l;
	int ii, jj, offset;
	T2 w;
	T ang;
	T2 a[8];
	int local_id = get_local_id( 0 );
	int group_id = get_group_id( 0 );
	d = group_id >> log_numBlocks;
	group_id = group_id & (numBlocks-1);
	Offset_in = mad24(group_id, batchSize, d << (pow1+pow2));
	tid = mul24(group_id, batchSize);
	i = tid >> lgStrideO;
	j = tid & (SO - 1);
	Offset_out = mad24(i, SI2, j + (d << (pow1+pow2)));
	b = group_id;
	tid = local_id;
	i = tid & (batchSize -1);
	j = tid >> lgBatchSize;
	Offset_in += mad24(j, SI1, i);
	in += Offset_in;
	for(int x =0; x<Radix1; x++)
		a[x] = in[x*SI1*Radix2];
	fftKernel(a, Radix1);
	if(Radix2 > 1){
		for(int x =1; x<Radix1; x++)
		{
			ang = 2.0f*M_PI*x/batchSize*j;
			w = (T2)(native_cos(ang), native_sin(ang));
			a[x] = cmplx_mul(a[x],w);
		}
		Offset_in = mad24(j, numIter*Radix2*batchSize, i);
		trans(a, &smem[tid], Radix2*batchSize, &smem[Offset_in],16,Radix2,NO_REV);
		fftKernel(a,Radix2);
	}
	Offset_out += mad24(j, numIter*SO, i);
	out += Offset_out;
	globalStores(a, out,SO,Radix1,NO_REV);
}
#endif
