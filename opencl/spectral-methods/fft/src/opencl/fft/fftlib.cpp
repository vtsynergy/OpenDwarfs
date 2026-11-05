#include <stdio.h>
#include <assert.h>
#include "Event.h"
#include "OptionParser.h"
#include "ResultDatabase.h"
#include "support.h"
#include "fftlib.h"
#include "OpenCLDeviceInfo.h"
#include <math.h>
#include "../../../../include/rdtsc.h"
#include "../../../../include/common_args.h"


Event fftEvent("FFT");
const char *cl_source_fft;
int Radix1, Radix2, SI;
static cl_kernel fftKrnl, fftKrnl1, fftKrnl2, fftKrnl0;
static cl_program fftProg;
static bool do_dp;

void setGlobalOption(string &arg, int fftn1, int fftn2)
{
	cl_int err;
	char *opt = (char *) malloc (sizeof(char)*500);
	int radixArr[10], Radix1Arr[10], Radix2Arr[10];
	int num_radix;
	getGlobalRadix(fftn2, radixArr, Radix1Arr, Radix2Arr, &num_radix);
	SI=fftn1*radixArr[0];
	int S0=SI;
	int batchSize = (fftn2==2048) ? 16 : 512;
	size_t wg_size[3];
	err = clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(size_t)*3, &wg_size, NULL);
	CHKERR(err, "Failed to get device info!");
	//	printf("wg_size: %d radixArr[0]: %d radixArr[1]: %d\n",wg_size[0],radixArr[0], radixArr[1]);
	batchSize = (batchSize > wg_size[0]) ? wg_size[0] : batchSize;	
	Radix1 = Radix1Arr[1];
	Radix2 = Radix2Arr[1];
	int S = radixArr[1]*fftn1;
	for(int i = 0; i < num_radix-1; i++)
		S *= radixArr[i];

	int numBlocks = SI / batchSize;
	sprintf(opt,"-D FFT_%d -D FFT_2D -D fftn1=%d -D pow1=%d -D fftn2=%d -D pow2=%d -D numBlocks=%d -D log_numBlocks=%d -D batchSize=%d -D lgBatchSize=%d -D lgStrideO=%d -D SO=%d -D SI1=%d -D Radix1=%d -D Radix2=%d -D SI2=%d -D S1=%d", fftn1, fftn1, (int)log2(fftn1),fftn2, (int)log2(fftn2),numBlocks,(int)log2(numBlocks),batchSize,(int)log2(batchSize),(int)log2(S0),S0, SI, Radix1, Radix2, S, radixArr[1]*fftn1);
	arg += opt;
	if(fftn2>128)
		arg+= " -D TWIDDLE";
	free(opt);
}

void createKernelWithSource()
{
	cl_int err;
	FILE *kernelFile=NULL;
	size_t kernelLength;
	size_t lengthRead;

    if(_deviceType == 3) 
	    kernelFile = fopen("fft.aocx", "rb");
    else
	    kernelFile = fopen("fft.cl", "r");

	if(kernelFile==NULL){
		printf("kernel file do not exist!");
		exit(0);}

	fseek(kernelFile, 0, SEEK_END);
	kernelLength = (size_t) ftell(kernelFile);
	cl_source_fft = (const char *) malloc(sizeof(char)*kernelLength);
	rewind(kernelFile);
	lengthRead = fread((void *) cl_source_fft, kernelLength, 1, kernelFile);
	fclose(kernelFile);
    if(_deviceType == 3) 
	    fftProg = clCreateProgramWithBinary(context, 1, &device_id, &kernelLength, (const unsigned char **) &cl_source_fft,NULL, &err);
    else
	    fftProg = clCreateProgramWithSource(context, 1, &cl_source_fft, &kernelLength, &err);
	CHKERR(err, "Failed to create program with source!");

}

void getLocalDimension(size_t &localsz, size_t &globalsz, int fftn1, int fftn2)
{
	cl_int err;
	size_t wg_size;
	err = clGetKernelWorkGroupInfo(fftKrnl, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &wg_size, NULL);
	CHKERR(err, "Failed to get kernel workgroup info!");
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
	cl_int err;
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


	err = clGetKernelWorkGroupInfo(fftKrnl, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), &wg_size, NULL);
	CHKERR(err, "Failed to get kernel workgroup info!");

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

	void
getGlobalRadix(int n, int *radix, int *Radix1, int *Radix2, int *num_radix)
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


	void 
init2(OptionParser& op, bool _do_dp, int fftn1, int fftn2)
{
	cl_int err;

	do_dp = _do_dp;

	ocd_initCL();

	createKernelWithSource();
	// ...and build it
	string args = " -cl-mad-enable ";

	setGlobalOption(args, fftn1, fftn2);
	printf(args.c_str());
    if(_deviceType == 3) 
	    err = clBuildProgram(fftProg, 1, &device_id, "-DOPENCL -I.", NULL, NULL);
    else 
	    err = clBuildProgram(fftProg, 0, NULL, args.c_str(), NULL, NULL);
#ifdef ENABLE_DEBUG
	{
		char* log = NULL;
		size_t bytesRequired = 0;
		err = clGetProgramBuildInfo(fftProg, 
				device_id, 
				CL_PROGRAM_BUILD_LOG,
				0,
				NULL,
				&bytesRequired );
		log = (char*)malloc( bytesRequired + 1 );
		err = clGetProgramBuildInfo(fftProg, 
				device_id, 
				CL_PROGRAM_BUILD_LOG,
				bytesRequired,
				log,
				NULL );
		//std::cout << log << std::endl;
		//		FILE *f;
		//		f=fopen("log.txt","w");
		//		fprintf(f,log);
		//		fclose(f);	

		free( log );
	}
#endif
	if (err != CL_SUCCESS) {
		char log[50000];
		size_t retsize = 0;
		err = clGetProgramBuildInfo(fftProg, device_id, CL_PROGRAM_BUILD_LOG,
				50000*sizeof(char),  log, &retsize);
		CHKERR(err, "Failed to get program build info!");
		cout << "Retsize: " << retsize << endl;
		cout << "Log: " << log << endl;
		dumpPTXCode(context, fftProg, "oclFFT");
		exit(-1);
	}
	else {
		// dumpPTXCode(context, fftProg, "oclFFT");
	}
	char kernel_name[20];
	sprintf(kernel_name,"fft1D_%d",fftn1);
	fftKrnl = clCreateKernel(fftProg, kernel_name, &err);
	fftKrnl0 = clCreateKernel(fftProg, "fft0", &err);
	fftKrnl1 = clCreateKernel(fftProg, "fft1", &err);
	fftKrnl2 = clCreateKernel(fftProg, "fft2", &err);
	CHKERR(err, "Failed to create kernels!");
}

	void 
init(OptionParser& op, bool _do_dp, int fftn)
{
	cl_int err;
	do_dp = _do_dp;

	ocd_initCL();

	createKernelWithSource();

	string args = " -cl-mad-enable ";

	char opt[400];
	sprintf(opt,"-D FFT_%d -D fftn1=%d -D pow1=%d", fftn, fftn, fftn/16, log2(fftn));
	args += opt;
	err = clBuildProgram(fftProg, 0, NULL, args.c_str(), NULL, NULL);
#ifdef ENABLE_DEBUG
	{
		char* log = NULL;
		size_t bytesRequired = 0;
		err = clGetProgramBuildInfo(fftProg, 
				device_id, 
				CL_PROGRAM_BUILD_LOG,
				0,
				NULL,
				&bytesRequired );
		log = (char*)malloc( bytesRequired + 1 );
		err = clGetProgramBuildInfo(fftProg, 
				device_id, 
				CL_PROGRAM_BUILD_LOG,
				bytesRequired,
				log,
				NULL );
		//		FILE *f;
		//		f=fopen("log.txt","w");
		//		fprintf(f,log);
		//		fclose(f);	
		free( log );
	}
#endif
	if (err != CL_SUCCESS) {
		char log[50000];
		size_t retsize = 0;
		err = clGetProgramBuildInfo(fftProg, device_id, CL_PROGRAM_BUILD_LOG,
				50000*sizeof(char),  log, &retsize);
		CL_CHECK_ERROR(err);
		cout << "Retsize: " << retsize << endl;
		cout << "Log: " << log << endl;
		dumpPTXCode(context, fftProg, "oclFFT");
		exit(-1);
	}
	else {
		// dumpPTXCode(context, fftProg, "oclFFT");
	}
	char kernel_name[20];
	sprintf(kernel_name,"fft1D_%d",fftn);
	fftKrnl = clCreateKernel(fftProg, kernel_name, &err);
	CHKERR(err, "Failed to create kernel!");
	fftKrnl1 = clCreateKernel(fftProg, "fft0", &err);
	CHKERR(err, "Failed to create kernel!");
}


	void 
forward(void* workp, void *temp, int n_ffts, int fftn)
{
	cl_int err;
	size_t localsz0;
	size_t globalsz0;
	size_t localsz;
	size_t globalsz;
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
	if(fftn> 2048)
	{
		clSetKernelArg(fftKrnl1, 0, sizeof(cl_mem), workp);
		clSetKernelArg(fftKrnl1, 1, sizeof(cl_mem), temp);
		clSetKernelArg(fftKrnl, 0, sizeof(cl_mem), temp);
		clSetKernelArg(fftKrnl, 1, sizeof(cl_mem), workp);
	}
	else{
		clSetKernelArg(fftKrnl, 0, sizeof(cl_mem), workp);
	}
	START_KERNEL
		if(fftn> 2048){
			//		printf("local size0: %d global size0 %d \n",localsz0,globalsz0);
			err = clEnqueueNDRangeKernel(commands, fftKrnl1, 1, NULL, 
					&globalsz0, &localsz0, 0, 
					NULL, &ocdTempEvent);
			err = clWaitForEvents(1, &ocdTempEvent);
			START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "FFT Kernels fftKrnl1", ocdTempTimer)
				END_TIMER(ocdTempTimer)

		}
	//	printf("local size: %d global size %d \n",localsz,globalsz);
	err = clEnqueueNDRangeKernel(commands, fftKrnl, 1, NULL, 
			&globalsz, &localsz, 0, 
			NULL, &ocdTempEvent);
	err = clWaitForEvents(1, &ocdTempEvent);
	START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "FFT Kernels fftKrnl", ocdTempTimer)
		END_TIMER(ocdTempTimer)
		CHKERR(err, "Failed to enqueue kernel!");
	END_KERNEL
}

	void 
forward2(void* workp, void* temp, int n_ffts, int fftn1, int fftn2)
{
	cl_int err;

	size_t localsz;
	size_t globalsz;
	size_t localsz0;
	size_t globalsz0;
	size_t localsz1;
	size_t globalsz1;
	size_t localsz2;
	size_t globalsz2;
	int i = 0;
	if(fftn1<4096){
		getLocalDimension(localsz,globalsz,fftn1, fftn2);
	}
	else{
		getGlobalDimension(localsz0, globalsz0, 1, fftn1, 0);
		globalsz0 = globalsz0 * localsz0*fftn2;
		getGlobalDimension(localsz, globalsz, 1, fftn1, 1);
		globalsz = globalsz * localsz*fftn2;
	}
	getGlobalDimension(localsz1, globalsz1, fftn1, fftn2, 0);
	globalsz1 = globalsz1 * localsz1;
	if(fftn2>128){
		getGlobalDimension(localsz2, globalsz2, fftn1, fftn2, 1);
		globalsz2 = globalsz2 * localsz2;
		clSetKernelArg(fftKrnl2, 0, sizeof(cl_mem), temp);
		clSetKernelArg(fftKrnl2, 1, sizeof(cl_mem), workp);
		clSetKernelArg(fftKrnl1, 0, sizeof(cl_mem), workp);
		clSetKernelArg(fftKrnl1, 1, sizeof(cl_mem), temp);
	}
	else
	{
		clSetKernelArg(fftKrnl1, 0, sizeof(cl_mem), workp);
		clSetKernelArg(fftKrnl1, 1, sizeof(cl_mem), workp);
	}
	if(fftn1>=4096){
		clSetKernelArg(fftKrnl, 0, sizeof(cl_mem), temp);
		clSetKernelArg(fftKrnl, 1, sizeof(cl_mem), workp);
		clSetKernelArg(fftKrnl0, 0, sizeof(cl_mem), workp);
		clSetKernelArg(fftKrnl0, 1, sizeof(cl_mem), temp);
	}
	else
	{
		clSetKernelArg(fftKrnl, 0, sizeof(cl_mem), workp);

	}
	START_KERNEL
		//Use a dual timer composed of single timers
		cl_event *firstEvent;
	if(fftn1>=4096){
		//	printf("local size0: %d global size0: %d \n",localsz0,globalsz0);
		err = clEnqueueNDRangeKernel(commands, fftKrnl0, 1, NULL,
				&globalsz0, &localsz0, 0,
				NULL, &fftEvent.CLEvent());
		firstEvent = &fftEvent.CLEvent();
		err = clWaitForEvents(1, &fftEvent.CLEvent());
		START_TIMER(fftEvent.CLEvent(), OCD_TIMER_KERNEL, "FFT Kernels fftKrnl0", ocdTempTimer)
		END_TIMER(ocdTempTimer)
		CHKERR(err, "Failed to enqueue kernel!");
	}

	//	printf("local size: %d global size: %d \n",localsz,globalsz);
	err = clEnqueueNDRangeKernel(commands, fftKrnl, 1, NULL, 
			&globalsz, &localsz, 0, 
			NULL, &fftEvent.CLEvent());
	if (fftn1<4096) firstEvent = &fftEvent.CLEvent();//inserted for dual timer support
	err = clWaitForEvents(1, &fftEvent.CLEvent());
	START_TIMER(fftEvent.CLEvent(), OCD_TIMER_KERNEL, "FFT Kernels fftKrnl1", ocdTempTimer)
	END_TIMER(ocdTempTimer)
	CHKERR(err, "Failed to wait for events!");

	//	printf("local size1: %d global size1: %d \n",localsz1,globalsz1);
	err = clEnqueueNDRangeKernel(commands, fftKrnl1, 1, NULL, 
			&globalsz1, &localsz1, 0, 
			NULL, &fftEvent.CLEvent());
	err = clWaitForEvents(1, &fftEvent.CLEvent());
	START_TIMER(fftEvent.CLEvent(), OCD_TIMER_KERNEL, "FFT Kernels fftKrnl1", ocdTempTimer)
	END_TIMER(ocdTempTimer)
	CHKERR(err, "Failed to wait for events!");

	if(fftn2>128){
		//	printf("local size2: %d global size2: %d \n",localsz2,globalsz2);
		err = clEnqueueNDRangeKernel(commands, fftKrnl2, 1, NULL, 
				&globalsz2, &localsz2, 0, 
				NULL, &fftEvent.CLEvent());
		err = clWaitForEvents(1, &fftEvent.CLEvent());
		START_TIMER(fftEvent.CLEvent(), OCD_TIMER_KERNEL, "FFT Kernels fftKrnl2", ocdTempTimer)
		END_TIMER(ocdTempTimer)
		CHKERR(err, "Failed to wait for events!");
	}
	//set a dual timer to check the entire range
	START_DUAL_TIMER(*firstEvent, fftEvent.CLEvent(), "FFT Kernels (Span)", ocdTempDualTimer)
	END_DUAL_TIMER(ocdTempDualTimer)
	END_KERNEL
}


#include <map>
static map<void*, cl_mem> memobjmap;

	void
allocHostBuffer(void** bufp, unsigned long bytes)
{
#if 1 // pinned memory?
	cl_int err;
	cl_mem memobj = clCreateBuffer(context,CL_MEM_READ_WRITE|CL_MEM_ALLOC_HOST_PTR,
			bytes, NULL, &err);
	CHKERR(err, "Failed to create buffer!");

	*bufp = clEnqueueMapBuffer(commands, memobj, true,
			CL_MAP_READ|CL_MAP_WRITE,
			0,bytes,0,NULL,NULL,&err);

	memobjmap[*bufp] = memobj;
	CHKERR(err, "Failed to enqueue map buffer!");
#else
	*bufp = malloc(bytes);
#endif
}


	void
freeHostBuffer(void* buf)
{
#if 1 // pinned memory?
	cl_int err;
	cl_mem memobj = memobjmap[buf];
	err = clReleaseMemObject(memobj);
	CHKERR(err, "Failed to release memo object!");
#else
	free(buf);
#endif
}


	void
allocDeviceBuffer(void** objp, unsigned long bytes)
{
	cl_int err;

	*(cl_mem**)objp = new cl_mem;
	**(cl_mem**)objp = clCreateBuffer(context, CL_MEM_READ_WRITE, bytes, 
			NULL, &err);
	CHKERR(err, "Failed to create buffer!");
}


	void
freeDeviceBuffer(void* buffer)
{
	clReleaseMemObject(*(cl_mem*)buffer);
}


	void
copyToDevice(void* to_device, void* from_host, unsigned long bytes)
{
	cl_int err = clEnqueueWriteBuffer(commands, *(cl_mem*)to_device, CL_TRUE, 
			0, bytes, from_host, 0, NULL, &ocdTempEvent);
	clFinish(commands);
	START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "FFT Data Copy", ocdTempTimer)
	END_TIMER(ocdTempTimer)
	CHKERR(err, "Failed to enqueue write buffer!");
}


	void
copyFromDevice(void* to_host, void* from_device, unsigned long bytes)
{
	cl_int err = clEnqueueReadBuffer(commands, *(cl_mem*)from_device, CL_TRUE, 
			0, bytes, to_host, 0, NULL, &ocdTempEvent);
	clFinish(commands);
	START_TIMER(ocdTempEvent, OCD_TIMER_D2H, "FFT Data Copy", ocdTempTimer)
	END_TIMER(ocdTempTimer)
	CHKERR(err, "Failed to enqueue read buffer!");
}

