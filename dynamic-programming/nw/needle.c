#define LIMIT -999
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "needle.h"
#include <sys/time.h>
#include "../../include/rdtsc.h"
#include "../../include/common_args.h"

#define CHECKERR(err) \
    if (err != CL_SUCCESS) \
    { \
        fprintf(stderr, "Error: %d\n", err);\
        exit(1); \
    }
//#define USEGPU 1  
////////////////////////////////////////////////////////////////////////////////
// declaration, forward
void runTest( int argc, char** argv);

int platform_id = PLATFORM_ID, n_device = DEVICE_ID;

int blosum62[24][24] = {
{ 4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4},
{-1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4},
{-2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4},
{-2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4},
{ 0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4},
{-1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4},
{-1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
{ 0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4},
{-2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4},
{-1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4},
{-1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4},
{-1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4},
{-1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4},
{-2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4},
{-1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4},
{ 1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4},
{ 0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4},
{-3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4},
{-2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4},
{ 0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4},
{-2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4},
{-1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4},
{ 0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4},
{-4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1}
};

double gettime() {
  struct timeval t;
  gettimeofday(&t,NULL);
  return t.tv_sec+t.tv_usec*1e-6;
}

int 
maximum( int a,
		 int b,
		 int c){

int k;
if( a <= b )
k = b;
else 
k = a;

if( k <=c )
return(c);
else
return(k);

}

////////////////////////////////////////////////////////////////////////////////
// Program main
////////////////////////////////////////////////////////////////////////////////
int
main( int argc, char** argv) 
{
    OCD_INIT
	runTest( argc, argv);
        OCD_FINISH
    return EXIT_SUCCESS;
}

void usage(int argc, char **argv)
{
	fprintf(stderr, "Usage: %s <max_rows/max_cols> <penalty> [platform & device]\n", argv[0]);
	fprintf(stderr, "\t<dimension>  - x and y dimensions\n");
	fprintf(stderr, "\t<penalty> - penalty(positive integer)\n");
	fprintf(stderr, "\t[platform] - index of platform)\n");
	fprintf(stderr, "\t[device] - index of device)\n");
	exit(1);
}

void runTest( int argc, char** argv) 
{
    int max_rows, max_cols, penalty;
    int *input_itemsets, *output_itemsets, *referrence;
	cl_mem matrix_cuda, matrix_cuda_out, referrence_cuda;
	int size;
	
    int i, j;
   
	ocd_init(&argc, &argv, NULL);
	ocd_options opts = ocd_get_options();
	platform_id = opts.platform_id;
	n_device = opts.device_id;

    // the lengths of the two sequences should be able to divided by 16.
	// And at current stage  max_rows needs to equal max_cols
	if (argc == 3)
	{
		max_rows = atoi(argv[1]);
		max_cols = atoi(argv[1]);
		penalty = atoi(argv[2]);
	}
    else{
	usage(argc, argv);
    }
	
	if(atoi(argv[1])%16!=0){
	fprintf(stderr,"The dimension values must be a multiple of 16\n");
	exit(1);
	}
	

	max_rows = max_rows + 1;
	max_cols = max_cols + 1;
	referrence = (int *)malloc( max_rows * max_cols * sizeof(int) );
    input_itemsets = (int *)malloc( max_rows * max_cols * sizeof(int) );
	output_itemsets = (int *)malloc( max_rows * max_cols * sizeof(int) );
	

	if (!input_itemsets)
		fprintf(stderr, "error: can not allocate memory");

    srand ( 7 );
	
	
    for (i = 0 ; i < max_cols; i++){
		for (j = 0 ; j < max_rows; j++){
			input_itemsets[i*max_cols+j] = 0;
		}
	}
	
	printf("Start Needleman-Wunsch\n");
	
	for(i=1; i< max_rows ; i++){    //please define your own sequence. 
       input_itemsets[i*max_cols] = rand() % 10 + 1;
	}
    for(j=1; j< max_cols ; j++){    //please define your own sequence.
       input_itemsets[j] = rand() % 10 + 1;
	}


	for (i = 1 ; i < max_cols; i++){
		for (j = 1 ; j < max_rows; j++){
		referrence[i*max_cols+j] = blosum62[input_itemsets[i*max_cols]][input_itemsets[j]];
		}
	}

    for(i = 1; i< max_rows ; i++)
       input_itemsets[i*max_cols] = -i * penalty;
	for(j = 1; j< max_cols ; j++)
       input_itemsets[j] = -j * penalty;

    cl_device_id clDevice;
    cl_context clContext;
    cl_command_queue clCommands;
    cl_program clProgram;
    cl_kernel clKernel_nw1;
    cl_kernel clKernel_nw2;

    cl_int errcode;

    FILE *kernelFile;
    char *kernelSource;
    size_t kernelLength;
	
    clDevice = GetDevice(platform_id, n_device);
  
 
    clContext = clCreateContext(NULL, 1, &clDevice, NULL, NULL, &errcode);
    CHECKERR(errcode);

    clCommands = clCreateCommandQueue(clContext, clDevice, CL_QUEUE_PROFILING_ENABLE, &errcode);
    CHECKERR(errcode);

    kernelFile = fopen("needle_kernel.cl", "r");
    fseek(kernelFile, 0, SEEK_END);
    kernelLength = (size_t) ftell(kernelFile);
    kernelSource = (char *) malloc(sizeof(char)*kernelLength);
    rewind(kernelFile);
    fread((void *) kernelSource, kernelLength, 1, kernelFile);
    fclose(kernelFile);

    clProgram = clCreateProgramWithSource(clContext, 1, (const char **) &kernelSource, &kernelLength, &errcode);
    CHECKERR(errcode);

    free(kernelSource);

    errcode = clBuildProgram(clProgram, 1, &clDevice, NULL, NULL, NULL);
    if (errcode == CL_BUILD_PROGRAM_FAILURE)
    {
        char *log;
        size_t logLength;
        errcode = clGetProgramBuildInfo(clProgram, clDevice, CL_PROGRAM_BUILD_LOG, 0, NULL, &logLength);
        log = (char *) malloc(sizeof(char)*logLength);
        errcode = clGetProgramBuildInfo(clProgram, clDevice, CL_PROGRAM_BUILD_LOG, logLength, (void *) log, NULL);
        fprintf(stderr, "Kernel build error! Log:\n%s", log);
        free(log);
        return;
    }
    CHECKERR(errcode);

    clKernel_nw1 = clCreateKernel(clProgram, "needle_opencl_shared_1", &errcode);
    CHECKERR(errcode);
    clKernel_nw2 = clCreateKernel(clProgram, "needle_opencl_shared_2", &errcode);
    CHECKERR(errcode);

    size = max_cols * max_rows;
    referrence_cuda = clCreateBuffer(clContext, CL_MEM_READ_ONLY, sizeof(int)*size, NULL, &errcode);
    CHECKERR(errcode);
    matrix_cuda = clCreateBuffer(clContext, CL_MEM_READ_WRITE, sizeof(int)*size, NULL, &errcode);
    CHECKERR(errcode);
    matrix_cuda_out = clCreateBuffer(clContext, CL_MEM_READ_WRITE, sizeof(int)*size, NULL, &errcode);
    CHECKERR(errcode);
    errcode = clEnqueueWriteBuffer(clCommands, referrence_cuda, CL_TRUE, 0, sizeof(int)*size, (void *) referrence, 0, NULL, &ocdTempEvent);

    START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "NW Reference Copy", ocdTempTimer)
    END_TIMER(ocdTempTimer)
    CHECKERR(errcode);
    errcode = clEnqueueWriteBuffer(clCommands, matrix_cuda, CL_TRUE, 0, sizeof(int)*size, (void *) input_itemsets, 0, NULL, &ocdTempEvent);
    clFinish(clCommands);
    START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "NW Item Set Copy", ocdTempTimer)
	CHECKERR(errcode);
	END_TIMER(ocdTempTimer)
    size_t localWorkSize[2] = {BLOCK_SIZE, 1};
    size_t globalWorkSize[2];
	int block_width = ( max_cols - 1 )/BLOCK_SIZE;

	printf("Processing top-left matrix\n");
	//process top-left matrix
	for(i = 1 ; i <= block_width ; i++){
        globalWorkSize[0] = i*localWorkSize[0];
        globalWorkSize[1] = 1*localWorkSize[1];
        errcode = clSetKernelArg(clKernel_nw1, 0, sizeof(cl_mem), (void *) &referrence_cuda);
        errcode |= clSetKernelArg(clKernel_nw1, 1, sizeof(cl_mem), (void *) &matrix_cuda);
        errcode |= clSetKernelArg(clKernel_nw1, 2, sizeof(cl_mem), (void *) &matrix_cuda_out);
        errcode |= clSetKernelArg(clKernel_nw1, 3, sizeof(int), (void *) &max_cols);
        errcode |= clSetKernelArg(clKernel_nw1, 4, sizeof(int), (void *) &penalty);
        errcode |= clSetKernelArg(clKernel_nw1, 5, sizeof(int), (void *) &i);
        errcode |= clSetKernelArg(clKernel_nw1, 6, sizeof(int), (void *) &block_width);
        CHECKERR(errcode);
        errcode = clEnqueueNDRangeKernel(clCommands, clKernel_nw1, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
        clFinish(clCommands);
	START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "NW Kernels", ocdTempTimer)
        END_TIMER(ocdTempTimer)
	CHECKERR(errcode);
	}
	printf("Processing bottom-right matrix\n");
    //process bottom-right matrix
	for(i = block_width - 1  ; i >= 1 ; i--){
        globalWorkSize[0] = i*localWorkSize[0];
        globalWorkSize[1] = 1*localWorkSize[1];
        errcode = clSetKernelArg(clKernel_nw2, 0, sizeof(cl_mem), (void *) &referrence_cuda);
        errcode |= clSetKernelArg(clKernel_nw2, 1, sizeof(cl_mem), (void *) &matrix_cuda);
        errcode |= clSetKernelArg(clKernel_nw2, 2, sizeof(cl_mem), (void *) &matrix_cuda_out);
        errcode |= clSetKernelArg(clKernel_nw2, 3, sizeof(int), (void *) &max_cols);
        errcode |= clSetKernelArg(clKernel_nw2, 4, sizeof(int), (void *) &penalty);
        errcode |= clSetKernelArg(clKernel_nw2, 5, sizeof(int), (void *) &i);
        errcode |= clSetKernelArg(clKernel_nw2, 6, sizeof(int), (void *) &block_width);
        CHECKERR(errcode);
	errcode = clEnqueueNDRangeKernel(clCommands, clKernel_nw2, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
        clFinish(clCommands);
	START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "NW Kernels", ocdTempTimer)
        END_TIMER(ocdTempTimer)
        CHECKERR(errcode);
	}

    errcode = clEnqueueReadBuffer(clCommands, matrix_cuda, CL_TRUE, 0, sizeof(float)*size, (void *) output_itemsets, 0, NULL, &ocdTempEvent);
    clFinish(clCommands);
    	START_TIMER(ocdTempEvent, OCD_TIMER_D2H, "NW Item Set Copy", ocdTempTimer)
	END_TIMER(ocdTempTimer)
	CHECKERR(errcode);
	
	
#ifdef TRACE

	printf("print traceback value GPU:\n");
    
	for (i = max_rows - 2,  j = max_rows - 2; i>=0, j>=0;){

		int nw, n, w, traceback;

		if ( i == max_rows - 2 && j == max_rows - 2 )
			printf("%d ", output_itemsets[ i * max_cols + j]); //print the first element
            

		if ( i == 0 && j == 0 )
           break;


		if ( i > 0 && j > 0 ){
			nw = output_itemsets[(i - 1) * max_cols + j - 1];
		    w  = output_itemsets[ i * max_cols + j - 1 ];
            n  = output_itemsets[(i - 1) * max_cols + j];
		}
		else if ( i == 0 ){
		    nw = n = LIMIT;
		    w  = output_itemsets[ i * max_cols + j - 1 ];
		}
		else if ( j == 0 ){
		    nw = w = LIMIT;
            n  = output_itemsets[(i - 1) * max_cols + j];
		}
		else{
		}

		traceback = maximum(nw, w, n);
		
		printf("%d ", traceback);

		if(traceback == nw )
		{i--; j--; continue;}

        else if(traceback == w )
		{j--; continue;}

        else if(traceback == n )
		{i--; continue;}

		else
		;
	}
    printf("\n");

#endif

	clReleaseMemObject(referrence_cuda);
	clReleaseMemObject(matrix_cuda);
	clReleaseMemObject(matrix_cuda_out);
    clReleaseKernel(clKernel_nw1);
    clReleaseKernel(clKernel_nw2);
    clReleaseProgram(clProgram);
    clReleaseCommandQueue(clCommands);
    clReleaseContext(clContext);

}

