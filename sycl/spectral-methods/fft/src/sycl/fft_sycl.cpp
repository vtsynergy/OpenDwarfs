// fft_sycl.cpp -- 1-D FFT driver, translated from fft.cpp and fftlib.cpp (SHOC / OpenDwarfs)
#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <string>
#include "dwarf_sycl.hpp"
#include "OptionParser.h"
#include "fft_kernels.hpp"
using namespace std;
using std::min;

struct cplxflt { float x; float y; };
int fft_valid = 1;   /* result of the 1-D validation, read by main() */
static int Radix1, Radix2;      /* globals, zero-initialised, as in fftlib.cpp */
static size_t g_wg_size;        /* CL_KERNEL_WORK_GROUP_SIZE of the OpenCL tree */

static void getLocalRadix(unsigned int n, unsigned int *radix, unsigned int *num_radix, unsigned int maxRadix);
static void getGlobalRadix(int n, int *radix, int *Radix1, int *Radix2, int *num_radix);

/* ---- work-size selection, verbatim from fftlib.cpp ---- */
static void getLocalDimension(size_t &localsz, size_t &globalsz, int fftn1, int fftn2)
{
	size_t wg_size;
	wg_size = g_wg_size;
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
static void getGlobalDimension(size_t  &localsz, size_t &globalsz, int BS, int n, int n_g)
{	
	size_t wg_size;
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
	wg_size = g_wg_size;

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
static void getLocalRadix(unsigned int n, unsigned int *radix, unsigned int *num_radix, unsigned int maxRadix)
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
static void getGlobalRadix(int n, int *radix, int *Radix1, int *Radix2, int *num_radix)
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



static sycl::queue &init(int fftn)
{
	static sycl::queue q = dsycl::make_queue();
	g_wg_size = q.get_device().get_info<sycl::info::device::max_work_group_size>();
	(void)fftn;
	dsycl::prebuild(q);   /* kernel JIT, as clBuildProgram: outside the timed region */
	return q;
}

static void forward(sycl::queue &q, void *workp_, void *temp_, int fftn)
{
	using fftk::T2;
	T2 *workp = (T2*)workp_, *temp = (T2*)temp_;
	size_t localsz0 = 0, globalsz0 = 0, localsz = 0, globalsz = 0;
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
	static bool shown = false;
	if(!shown) { printf("FFT launch: stage0 global=%zu local=%zu; stage1 global=%zu local=%zu\n", globalsz0, localsz0, globalsz, localsz); shown = true; }
	switch(fftn) {
		case 128:  fftk::fft1D_128(q, globalsz, localsz, workp); break;
		case 256:  fftk::fft1D_256(q, globalsz, localsz, workp); break;
		case 512:  fftk::fft1D_512(q, globalsz, localsz, workp); break;
		case 1024: fftk::fft1D_1024(q, globalsz, localsz, workp); break;
		case 2048: fftk::fft1D_2048(q, globalsz, localsz, workp); break;
		case 4096: fftk::fft0<4096,12>(q, globalsz0, localsz0, workp, temp); fftk::fft1D_4096(q, globalsz, localsz, temp, workp); break;
		case 8192: fftk::fft0<8192,13>(q, globalsz0, localsz0, workp, temp); fftk::fft1D_8192(q, globalsz, localsz, temp, workp); break;
	}
}

template <class T2> void dump1D(OptionParser& op);

void addBenchmarkSpecOptions(OptionParser &op)
{
	op.addOption("pts", OPT_INT, "0", "data size (in megabytes)");
	op.addOption("pts1", OPT_INT, "0", "data size (in megabytes)");
	op.addOption("pts2", OPT_INT, "0", "data size (in megabytes)");
	op.addOption("2D", OPT_BOOL, "false", "2D FFT");
}

void RunBenchmark(OptionParser &op)
{
	if(op.getOptionBool("2D")) { cerr << "--2D is not available in the SYCL tree (the OpenCL 2-D path is unverified; FIXES.md)\n"; exit(1); }
	dump1D<cplxflt>(op);
}

	template <class T2> 
void dump1D(OptionParser& op)
{	
	int i;
	int fftn;
	char *work, *temp;
	T2* source, * result;
	unsigned long bytes = 0;

	int probSizes[7] = { 128, 256, 512, 1024, 2048, 4096, 8192 };
	int sizeIndex = op.getOptionInt("pts")-1;
	if (sizeIndex < 0 || sizeIndex >= 7) {
		cerr << "Invalid size index specified\n";
		exit(-1);
	}
	fftn = probSizes[sizeIndex];

	// Convert to MB
	unsigned long used_bytes = fftn * sizeof(T2);

	sycl::queue &q = init(fftn);

	// now determine how much available memory will be used
	//int half_n_ffts = bytes / (fftn*sizeof(T2)*2);
	int n_ffts = 1;
	double N = fftn;

	fprintf(stdout, "used_bytes=%lu, N=%g\n", used_bytes, N);

	// allocate host and device memory
	source = (T2*)sycl::malloc_host(used_bytes, q);
	result = (T2*)sycl::malloc_host(used_bytes, q);

	/* [FIX-OCL-26] the input came from rand() without a seed. The OpenCL runtimes (NVIDIA driver, PoCL, ...) may consume rand()
	 * themselves while loading, so the "random" input differed between implementations and runs and the result could never be
	 * compared or checked. Seed explicitly right before generating the data. */
	srand(1);
	// init host memory...
	for (i = 0; i < N; i++) {
		source[i].x = (rand()/(float)RAND_MAX)*2-1;
		source[i].y = (rand()/(float)RAND_MAX)*2-1;
	}

	/* common timer: alloc + H2D + forward transform + D2H per repetition (JIT is in init() above) */
	dwarf_init("fft", "spectral-methods", "sycl");
	dwarf_set_size(fftn, 1);
	for (int rep = dwarf_first_rep(); rep < dwarf_reps(); rep++) {
	dwarf_region_begin();
	// alloc device memory
	work = dsycl::alloc<char>(q, used_bytes);
	temp = dsycl::alloc<char>(q, used_bytes);
	dsycl::write(q, work, source, used_bytes);


	forward(q, work, temp, fftn);

	dsycl::read(q, result, work, used_bytes);
	dwarf_region_end(rep);
	if (rep != dwarf_reps() - 1) { dsycl::release(q, work); dsycl::release(q, temp); }
	}
#ifdef PRINT_RESULT
	for (i = 0; i < N; i++) {
		fprintf(stdout, "data[%d] (%g, %g)\n", i, result[i].x, result[i].y);
	}
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
	dsycl::release(q, work);
	dsycl::release(q, temp);
	sycl::free(source, q);
	sycl::free(result, q);
}
