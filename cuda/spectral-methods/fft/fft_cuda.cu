/* CUDA version of the 1-D FFT benchmark: opencl/spectral-methods/fft (fft.cpp dump1D + fftlib.cpp). The SHOC host scaffolding
 * (OptionParser, device enumeration) is replaced by a minimal parser for the same options (--pts N, --passes N; --passes does not
 * change the 1-D path in OpenCL either). The 2-D path is not translated. */
#include <algorithm>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "dwarf_timer.h"
#include "dwarf_cuda.h"
#include "fft_kernel.cuh"

/* NVIDIA OpenCL reports 256 for these kernels */
#define FFT_KERNEL_WG_SIZE 256
struct cplxflt { float x; float y; };
static int Radix1 = 0, Radix2 = 0; /* globals of fftlib.cpp; only the 2-D set-up assigns them */
void getLocalRadix(unsigned int n, unsigned int *radix, unsigned int *num_radix, unsigned int maxRadix);
void getGlobalRadix(int n, int *radix, int *Radix1, int *Radix2, int *num_radix);
static inline int imin_(int a, int b) { return a < b ? a : b; }
#define min(a,b) imin_((int)(a),(int)(b))

void getLocalRadix(unsigned int n, unsigned int *radix, unsigned int *num_radix, unsigned int maxRadix)
{
	if(maxRadix > 1)
	{
		maxRadix = min(n, maxRadix);
		unsigned int cnt = 0;
		while(n > maxRadix)
		{
			radix[cnt++] = maxRadix;
			n /= maxRadix;
		}
		radix[cnt++] = n;
		*num_radix = cnt;
		return;
	}
	int i = 0;
	while(n > 1)
	{
		if(n > 8)
		{
			n = n / 8;
			radix[i] = 8;
			i++;			
		}
		else if(n >= 4)
		{
			n = n / 4;
			radix[i] = 4;
			i++;
		}
		else
		{
			radix[i] = n;
			i++;
			*num_radix = i;
			return;
		}
	} 
}

void getGlobalRadix(int n, int *radix, int *Radix1, int *Radix2, int *num_radix)
{
	int baseRadix = min(n, 128);

	int numR = 0;
	int N = n;
	while(N > baseRadix) 
	{
		N /= baseRadix;
		numR++;
	}

	for(int i = 0; i < numR; i++)
		radix[i] = baseRadix;

	radix[numR] = N;
	numR++;
	*num_radix = numR;

	for(int i = 0; i < numR; i++)
	{
		int B = radix[i];
		if(B <= 8)
		{
			Radix1[i] = B;
			Radix2[i] = 1;
			continue;
		}

		int r1 = 2; 
		int r2 = B / r1;
		while(r2 > r1)
		{
			r1 *=2;
			r2 = B / r1;
		}
		Radix1[i] = r1;
		Radix2[i] = r2;
	}	
}

void getLocalDimension(size_t &localsz, size_t &globalsz, int fftn1, int fftn2)
{
	size_t wg_size = FFT_KERNEL_WG_SIZE; /* CL_KERNEL_WORK_GROUP_SIZE */
	unsigned int radix[5];
	unsigned int numRadix;
	getLocalRadix(fftn1, radix, &numRadix, 0);
	if(fftn1/radix[0] > wg_size)
		getLocalRadix(fftn1, radix, &numRadix, 32);
	//	for(int i = 0; i < numRadix; i++)
	//	{      
	//		assert( radix[i] && !( (radix[i] - 1) & radix[i] ) );
	//	}

	int batchsize = Radix2 == 1 ? wg_size : min(1, 16);
	int num_item_per_wg = ((fftn1/radix[0]) <=64) ? 64 : (fftn1/radix[0]);
	int wg_num = num_item_per_wg/(fftn1/radix[0]);
	localsz = num_item_per_wg;
	if(fftn1<4096)
	{
		batchsize = batchsize * fftn2;
		wg_num = batchsize % wg_num ? batchsize / wg_num + 1 : batchsize/ wg_num;
		globalsz = wg_num * localsz;
	}
	else
	{
		globalsz = wg_num * num_item_per_wg * fftn2;
	}

}

void getGlobalDimension(size_t  &localsz, size_t &globalsz, int BS, int n, int n_g)
{	
	size_t wg_size = FFT_KERNEL_WG_SIZE; /* CL_KERNEL_WORK_GROUP_SIZE */
	int radixArr[5], Radix1Arr[5], Radix2Arr[5];
	int radix, Radix1, Radix2;
	int num_radix;
	bool v = (BS==1) ? 0 : 1;	
	int R = v ? BS : 1;
	int SI = R;
	getGlobalRadix(n, radixArr, Radix1Arr, Radix2Arr, &num_radix);
	for(int i = 0; i < num_radix; i++)
		if(n_g != i)
			SI *= radixArr[i];

	Radix1 = Radix1Arr[n_g];
	Radix2 = Radix2Arr[n_g];



	int batchSize = 16;

	batchSize = v ? min(BS, batchSize) : batchSize;
	batchSize = Radix2 == 1 ? wg_size : batchSize;
	batchSize = min(batchSize, SI);

	int threads = batchSize * Radix2;
	threads = threads > wg_size ? wg_size : threads;
	batchSize = threads / Radix2;

	int blocks = SI / batchSize;
	if(!v)
		blocks *= BS;
	else
		blocks *= 1;
	globalsz=blocks;
	localsz=threads;
	//	printf("globalsz: %d localsz: %d\n",globalsz,localsz);

}
#undef min

#define FFT_LAUNCH(kern, global_, local_, ...) do { kern<<<(unsigned int)((global_)/(local_)), (unsigned int)(local_)>>>(__VA_ARGS__); CUDA_SYNC_LAUNCH(); } while (0)

static void forward(float2 *workp, float2 *temp, int fftn)
{
	size_t localsz0 = 0, globalsz0 = 0, localsz, globalsz;
	if(fftn <= 2048){
		getLocalDimension(localsz,globalsz,fftn, 1);
	}
	else
	{
		getGlobalDimension(localsz0, globalsz0, 1, fftn, 0);
		globalsz0 = globalsz0*localsz0;
		getGlobalDimension(localsz, globalsz, 1, fftn,1);
		globalsz = globalsz*localsz;
	}
	static int shown = 0;
	if(!shown) { shown = 1; if(fftn > 2048) printf("NDRange fft0: global %zu local %zu\n", globalsz0, localsz0); printf("NDRange fft1D_%d: global %zu local %zu\n", fftn, globalsz, localsz); }
	switch(fftn) {
		case 128:  FFT_LAUNCH(fft1D_128, globalsz, localsz, workp); break;
		case 256:  FFT_LAUNCH(fft1D_256, globalsz, localsz, workp); break;
		case 512:  FFT_LAUNCH(fft1D_512, globalsz, localsz, workp); break;
		case 1024: FFT_LAUNCH(fft1D_1024, globalsz, localsz, workp); break;
		case 2048: FFT_LAUNCH(fft1D_2048, globalsz, localsz, workp); break;
		case 4096: FFT_LAUNCH((fft0<4096,12>), globalsz0, localsz0, workp, temp); FFT_LAUNCH(fft1D_4096, globalsz, localsz, temp, workp); break;
		case 8192: FFT_LAUNCH((fft0<8192,13>), globalsz0, localsz0, workp, temp); FFT_LAUNCH(fft1D_8192, globalsz, localsz, temp, workp); break;
	}
}

int main(int argc, char *argv[])
{
	dwarf_cuda_init(&argc, &argv);
	int pts = 0, i, fftn, fft_valid = 1;
	for (i = 1; i < argc; i++) {
		if (!strcmp(argv[i], "--pts") && i + 1 < argc) pts = atoi(argv[++i]);
		else if (!strcmp(argv[i], "--passes") && i + 1 < argc) ++i;
		else { fprintf(stderr, "unsupported option %s (supported: --pts 1..7 --passes N; the 2-D path is not translated)\n", argv[i]); return 1; }
	}
	int probSizes[7] = { 128, 256, 512, 1024, 2048, 4096, 8192 };
	if (pts < 1 || pts > 7) { fprintf(stderr, "Invalid size index specified\n"); return 1; }
	fftn = probSizes[pts-1];
	unsigned long used_bytes = fftn * sizeof(cplxflt);
	double N = fftn;
	fprintf(stdout, "used_bytes=%lu, N=%g\n", used_bytes, N);

	cplxflt *source, *result;
	CUDA_CHECK(cudaMallocHost((void**)&source, used_bytes)); /* CL_MEM_ALLOC_HOST_PTR + map */
	CUDA_CHECK(cudaMallocHost((void**)&result, used_bytes));

	/* [FIX-OCL-26] seed right before generating the input */
	srand(1);
	for (i = 0; i < N; i++) {
		source[i].x = (rand()/(float)RAND_MAX)*2-1;
		source[i].y = (rand()/(float)RAND_MAX)*2-1;
	}

	/* common timer: alloc + H2D + forward transform + D2H per repetition */
	dwarf_init("fft", "spectral-methods", "cuda");
	dwarf_set_size(fftn, 1);
	float2 *work = 0, *temp = 0;
	for (int rep = dwarf_first_rep(); rep < dwarf_reps(); rep++) {
		dwarf_region_begin();
		CUDA_CHECK(cudaMalloc((void**)&work, used_bytes));
		CUDA_CHECK(cudaMalloc((void**)&temp, used_bytes));
		CUDA_CHECK(cudaMemcpy(work, source, used_bytes, cudaMemcpyHostToDevice));
		CUDA_CHECK(cudaDeviceSynchronize());
		forward(work, temp, fftn);
		CUDA_CHECK(cudaMemcpy(result, work, used_bytes, cudaMemcpyDeviceToHost));
		CUDA_CHECK(cudaDeviceSynchronize());
		dwarf_region_end(rep);
		if (rep != dwarf_reps() - 1) { cudaFree(work); cudaFree(temp); }
	}
#ifdef PRINT_RESULT
	for (i = 0; i < N; i++) fprintf(stdout, "data[%d] (%g, %g)\n", i, result[i].x, result[i].y);
#endif
	/* [FIX-OCL-27] no validation existed. Added: the transform of this very input in double precision (naive DFT); the SHOC
	 * kernels emit the spectrum in a radix-permuted order and with swapped real/imaginary parts, so the check compares the
	 * sorted magnitudes (permutation- and swap-insensitive), relative RMS error < 1e-3. */
	{
		double *ref = (double*) malloc(sizeof(double) * fftn), *got = (double*) malloc(sizeof(double) * fftn);
		for (int k = 0; k < fftn; k++) {
			double sr = 0, si = 0;
			for (int t = 0; t < fftn; t++) {
				double ang = -2.0 * M_PI * (double)((long)k * t % fftn) / fftn;
				sr += source[t].x * cos(ang) - source[t].y * sin(ang);
				si += source[t].x * sin(ang) + source[t].y * cos(ang);
			}
			ref[k] = sqrt(sr * sr + si * si);
			got[k] = sqrt((double)result[k].x * result[k].x + (double)result[k].y * result[k].y);
		}
		std::sort(ref, ref + fftn); std::sort(got, got + fftn);
		double e = 0, nrm = 0;
		for (int k = 0; k < fftn; k++) { e += (ref[k] - got[k]) * (ref[k] - got[k]); nrm += ref[k] * ref[k]; }
		double rel = sqrt(e / nrm);
		fft_valid = (rel < 1e-3) && (rel == rel);
		printf("%s: sorted-magnitude relative RMS error vs the double-precision DFT of the input: %.3g\n", fft_valid ? "SUCCESS" : "FAILURE", rel);
		free(ref); free(got);
	}
	CUDA_CHECK(cudaFree(work));
	CUDA_CHECK(cudaFree(temp));
	CUDA_CHECK(cudaFreeHost(source));
	CUDA_CHECK(cudaFreeHost(result));
	return dwarf_report(fft_valid);
}
