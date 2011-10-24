#include <cfloat>
#include <iostream>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include <unistd.h>
#include <time.h>
#include "support.h"
#include "ResultDatabase.h"
#include "Event.h"
#include "OptionParser.h"
#include "Timer.h"

#include "fftlib.h"
using namespace std;

// ****************************************************************************
// Function: addBenchmarkSpecOptions
//
// Purpose:
//   Add benchmark specific options parsing.  The user is allowed to specify
//   the size of the input data in megabytes.
//
// Arguments:
//   op: the options parser / parameter database
//
// Programmer: Collin McCurdy
// Creation: September 08, 2009
// Returns:  nothing
//
// ****************************************************************************
	void 
addBenchmarkSpecOptions(OptionParser &op) 
{
	op.addOption("pts", OPT_INT, "0", "data size (in megabytes)");
	op.addOption("pts1", OPT_INT, "0", "data size (in megabytes)");
	op.addOption("pts2", OPT_INT, "0", "data size (in megabytes)");
	op.addOption("2D", OPT_BOOL, "false", "2D FFT");
}


// ****************************************************************************
// Function: RunBenchmark
//
// Purpose:
//   Calls single precision and, if viable, double precision FFT
//   benchmark.  Optionally dumps data arrays for correctness check.
//
// Arguments:
//  resultDB: the benchmark stores its results in this ResultDatabase
//  op: the options parser / parameter database
//
// Returns:  nothing
//
// Programmer: Collin McCurdy
// Creation: September 08, 2009
//
// Modifications:
//
// ****************************************************************************

template <class T2> void runTest(const string& name, cl_device_id id, 
		cl_context ctx, cl_command_queue queue,
		ResultDatabase &resultDB, OptionParser& op);
template <class T2> void dump1D(OptionParser& op);
template <class T2> void dump2D(OptionParser& op);


	void
RunBenchmark(OptionParser &op)
{
	// convert from C++ bindings to C bindings
	// TODO propagate use of C++ bindings



	if(op.getOptionBool("2D"))
		dump2D<cplxflt>(op);
	else
		dump1D<cplxflt>(op);

}


// ****************************************************************************
// Function: runTest
//
// Purpose:
//   This benchmark measures the performance of a single or double
//   precision fast fourier transform (FFT).  Data transfer time over
//   the PCIe bus is not included in this measurement.
//
// Arguments:
//  resultDB: the benchmark stores its results in this ResultDatabase
//  op: the options parser / parameter database
//
// Returns:  nothing
//
// Programmer: Collin McCurdy
// Creation: September 08, 2009
//
// Modifications:
//   Jeremy Meredith, Thu Aug 19 15:43:55 EDT 2010
//   Added PCIe timings.  Added size index bounds check.
//
// ****************************************************************************

template <class T2> inline bool dp(void);
template <> inline bool dp<cplxflt>(void) { return false; }
template <> inline bool dp<cplxdbl>(void) { return true; }


// ****************************************************************************
// Function: dump
//
// Purpose:
//   Dump result array to stdout after FFT and IFFT.  For correctness 
//   checking.
//
// Arguments:
//  op: the options parser / parameter database
//
// Returns:  nothing
//
// Programmer: Collin McCurdy
// Creation: September 30, 2010
//
// Modifications:
//
// ****************************************************************************
	template <class T2> 
void dump2D(OptionParser& op)
{	
	int i;
	void* work, *temp;
	T2* source, * result;
	unsigned long bytes = 0;

	int probSizes[7] = { 128, 256, 512, 1024, 2048, 4096, 8192};
	int sizeIndex = op.getOptionInt("pts1")-1;
	int sizeIndey = op.getOptionInt("pts2")-1;
	if (sizeIndex < 0 || sizeIndex >= 7) {
		cerr << "Invalid size index specified\n";
		exit(-1);
	}
	if (sizeIndey < 0 || sizeIndey >= 7) {
		cerr << "Invalid size index specified\n";
		exit(-1);
	}

	int FFTN1=probSizes[sizeIndex],FFTN2=probSizes[sizeIndey];
	//int FFTN1=8192,FFTN2=512;
	unsigned long used_bytes = FFTN1*FFTN2*sizeof(T2);

	bool do_dp = dp<T2>();
	init2(op, do_dp, FFTN1, FFTN2);

	int n_ffts = 1;
	double N = FFTN1*FFTN2;


	// allocate host and device memory
	allocHostBuffer((void**)&source, used_bytes);
	allocHostBuffer((void**)&result, used_bytes);

	// init host memory...
	for (i = 0; i < N; i++) {
		source[i].x = (rand()/(float)RAND_MAX)*2-1;
		source[i].y = (rand()/(float)RAND_MAX)*2-1;
	}

	// alloc device memory
	allocDeviceBuffer(&work, used_bytes);
	allocDeviceBuffer(&temp, used_bytes);

	copyToDevice(work, source, used_bytes);

	forward2(work, temp, n_ffts, FFTN1, FFTN2);
	copyFromDevice(result, work, used_bytes);

#ifdef PRINT_RESULT
	for (i = 0; i < N; i++) {
		fprintf(stdout, "data[%d] (%g, %g) \n",i, result[i].x, result[i].y);
	}
#endif
	freeDeviceBuffer(work);
	freeDeviceBuffer(temp);
	freeHostBuffer(source);
	freeHostBuffer(result);
}

	template <class T2> 
void dump1D(OptionParser& op)
{	
	int i;
	int fftn;
	void* work, *temp;
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

	bool do_dp = dp<T2>();
	init(op, do_dp, fftn);

	// now determine how much available memory will be used
	//int half_n_ffts = bytes / (fftn*sizeof(T2)*2);
	int n_ffts = 1;
	double N = fftn;

	fprintf(stderr, "used_bytes=%lu, N=%g\n", used_bytes, N);

	// allocate host and device memory
	allocHostBuffer((void**)&source, used_bytes);
	allocHostBuffer((void**)&result, used_bytes);

	// init host memory...
	for (i = 0; i < N; i++) {
		source[i].x = (rand()/(float)RAND_MAX)*2-1;
		source[i].y = (rand()/(float)RAND_MAX)*2-1;
	}

	// alloc device memory
	allocDeviceBuffer(&work, used_bytes);
	allocDeviceBuffer(&temp, used_bytes);
	copyToDevice(work, source, used_bytes);


	forward(work, temp, n_ffts, fftn);

	copyFromDevice(result, work, used_bytes);
#ifdef PRINT_RESULT
	for (i = 0; i < N; i++) {
		fprintf(stdout, "data[%d] (%g, %g)\n", i, result[i].x, result[i].y);
	}
#endif
	freeDeviceBuffer(work);
	freeDeviceBuffer(temp);
	freeHostBuffer(source);
	freeHostBuffer(result);
}
