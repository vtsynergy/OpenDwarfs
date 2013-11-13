#define LIMIT -999
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "needle.h"
#include <sys/time.h>
#include "../../include/rdtsc.h"
#include "../../include/common_args.h"

//#define TRACE

void runTest( int argc, char** argv);

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
	ocd_init(&argc, &argv, NULL);
	ocd_initCL();
	runTest( argc, argv);
    ocd_finalize();
    return EXIT_SUCCESS;
}

void usage(int argc, char **argv)
{
	fprintf(stderr, "Usage: %s <max_rows/max_cols> <penalty> \n", argv[0]);
	fprintf(stderr, "\t<dimension>  - x and y dimensions\n");
	fprintf(stderr, "\t<penalty> - penalty(positive integer)\n");
	exit(1);
}

void runTest( int argc, char** argv) 
{
    int max_rows, max_cols, penalty;
    int *input_itemsets, *output_itemsets, *referrence;
	cl_mem matrix_cuda,referrence_cuda;
	int size;
    int i, j;
   
   	cl_int errcode;
   	
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
	//Initialize to zero
    for (i = 0 ; i < max_cols; i++){
		for (j = 0 ; j < max_rows; j++){
			input_itemsets[i*max_cols+j] = 0;
		}
	}
	
	printf("Start Needleman-Wunsch\n");
	//Used for input itemsets, then after the referrence matrix is constructed it is used/rewritten to store scores.
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

	//Fill first row and first column with initial scores -10, -20, -30, ... (upper left corner is set to zero).
    for(i = 1; i< max_rows ; i++)
       input_itemsets[i*max_cols] = -i * penalty;
	for(j = 1; j< max_cols ; j++)
       input_itemsets[j] = -j * penalty;

    //cl_device_id device_id;
    //cl_context context;
    //cl_command_queue commands;
    cl_program clProgram;
    cl_kernel clKernel_nw1;
    cl_kernel clKernel_nw2;

    
    FILE *kernelFile;
    char *kernelSource;
    size_t kernelLength;
	
    kernelFile = fopen("needle_kernel.cl", "r");
    fseek(kernelFile, 0, SEEK_END);
    kernelLength = (size_t) ftell(kernelFile);
    kernelSource = (char *) malloc(sizeof(char)*kernelLength);
    rewind(kernelFile);
    fread((void *) kernelSource, kernelLength, 1, kernelFile);
    fclose(kernelFile);

    clProgram = clCreateProgramWithSource(context, 1, (const char **) &kernelSource, &kernelLength, &errcode);
  	CHKERR(errcode, "Failed to create program with source!");

    free(kernelSource);

    errcode = clBuildProgram(clProgram, 1, &device_id, NULL, NULL, NULL);
    if (errcode == CL_BUILD_PROGRAM_FAILURE)
    {
        char *log;
        size_t logLength;
        errcode = clGetProgramBuildInfo(clProgram, device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &logLength);
        log = (char *) malloc(sizeof(char)*logLength);
        errcode = clGetProgramBuildInfo(clProgram, device_id, CL_PROGRAM_BUILD_LOG, logLength, (void *) log, NULL);
        fprintf(stderr, "Kernel build error! Log:\n%s", log);
        free(log);
        return;
    }
  	CHKERR(errcode, "Failed to get program build info!");

    clKernel_nw1 = clCreateKernel(clProgram, "needle_opencl_shared_1", &errcode);
  	CHKERR(errcode, "Failed to create kernel!");
    clKernel_nw2 = clCreateKernel(clProgram, "needle_opencl_shared_2", &errcode);
  	CHKERR(errcode, "Failed to create kernel!");

    size = max_cols * max_rows;
    referrence_cuda = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(int)*size, NULL, &errcode);
  	CHKERR(errcode, "Failed to create buffer!");
    matrix_cuda = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(int)*size, NULL, &errcode);
  	CHKERR(errcode, "Failed to create buffer!");
  	
    errcode = clEnqueueWriteBuffer(commands, referrence_cuda, CL_TRUE, 0, sizeof(int)*size, (void *) referrence, 0, NULL, &ocdTempEvent);
    clFinish(commands);
    START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "NW Reference Copy", ocdTempTimer)
    END_TIMER(ocdTempTimer)
  	CHKERR(errcode, "Failed to enqueue write buffer!");
  	
    errcode = clEnqueueWriteBuffer(commands, matrix_cuda, CL_TRUE, 0, sizeof(int)*size, (void *) input_itemsets, 0, NULL, &ocdTempEvent);
    clFinish(commands);
    START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "NW Item Set Copy", ocdTempTimer)
	END_TIMER(ocdTempTimer)
	CHKERR(errcode, "Failed to enqueue write buffer!");
	
    size_t localWorkSize[2] = {BLOCK_SIZE, 1}; //BLOCK_SIZE work items per work-group in 1D only.
    size_t globalWorkSize[2];
	int block_width = ( max_cols - 1 )/BLOCK_SIZE;

	printf("Processing top-left matrix\n");
	//process top-left matrix
	//Does what the 1st kernel loop does in a higher (block) level. i.e., takes care of blocks of BLOCK_SIZExBLOCK_SIZE in a wave-front pattern upwards
	//the main anti-diagonal (on block-level).
	//Each iteration takes care of 1, 2, 3, ... blocks that can be computed in parallel w/o dependencies
	//E.g. first block [0][0], then blocks [0][1] and [1][0], then [0][2], [1][1], [2][0], etc.
	for(i = 1 ; i <= block_width ; i++){
        globalWorkSize[0] = i*localWorkSize[0]; //i.e., for 1st iteration BLOCK_SIZE total (=1 W.G.), for 2nd iteration 2*BLOCK_SIZE total work items
        // (=2 W.G.)
        globalWorkSize[1] = 1*localWorkSize[1];
        errcode = clSetKernelArg(clKernel_nw1, 0, sizeof(cl_mem), (void *) &referrence_cuda);
        errcode |= clSetKernelArg(clKernel_nw1, 1, sizeof(cl_mem), (void *) &matrix_cuda);
        errcode |= clSetKernelArg(clKernel_nw1, 2, sizeof(int), (void *) &max_cols);
        errcode |= clSetKernelArg(clKernel_nw1, 3, sizeof(int), (void *) &penalty);
        errcode |= clSetKernelArg(clKernel_nw1, 4, sizeof(int), (void *) &i);
        errcode |= clSetKernelArg(clKernel_nw1, 5, sizeof(int), (void *) &block_width);
  		CHKERR(errcode, "Failed to set kernel arguments!");
        errcode = clEnqueueNDRangeKernel(commands, clKernel_nw1, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
        clFinish(commands);
		START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "NW Kernel nw1", ocdTempTimer)
        END_TIMER(ocdTempTimer)
  		CHKERR(errcode, "Failed to enqueue kernel!");
	}
	printf("Processing bottom-right matrix\n");
    //process bottom-right matrix
    //Does what the 2nd kernel loop does in a higher (block) level. i.e., takes care of blocks of BLOCK_SIZExBLOCK_SIZE in a wave-front pattern downwards
    //the main anti-diagonal.
	//Each iteration takes care of ..., 3, 2, 1 blocks that can be computed in parallel w/o dependencies
	for(i = block_width - 1  ; i >= 1 ; i--){
        globalWorkSize[0] = i*localWorkSize[0];
        globalWorkSize[1] = 1*localWorkSize[1];
        errcode = clSetKernelArg(clKernel_nw2, 0, sizeof(cl_mem), (void *) &referrence_cuda);
        errcode |= clSetKernelArg(clKernel_nw2, 1, sizeof(cl_mem), (void *) &matrix_cuda);
        errcode |= clSetKernelArg(clKernel_nw2, 2, sizeof(int), (void *) &max_cols);
        errcode |= clSetKernelArg(clKernel_nw2, 3, sizeof(int), (void *) &penalty);
        errcode |= clSetKernelArg(clKernel_nw2, 4, sizeof(int), (void *) &i);
        errcode |= clSetKernelArg(clKernel_nw2, 5, sizeof(int), (void *) &block_width);
  		CHKERR(errcode, "Failed to set kernel arguments!");
		errcode = clEnqueueNDRangeKernel(commands, clKernel_nw2, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
        clFinish(commands);
		START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "NW Kernel nw2", ocdTempTimer)
        END_TIMER(ocdTempTimer)
  		CHKERR(errcode, "Failed to enqueue kernel!");
	}

    errcode = clEnqueueReadBuffer(commands, matrix_cuda, CL_TRUE, 0, sizeof(float)*size, (void *) output_itemsets, 0, NULL, &ocdTempEvent);
    clFinish(commands);
    START_TIMER(ocdTempEvent, OCD_TIMER_D2H, "NW Item Set Copy", ocdTempTimer)
	END_TIMER(ocdTempTimer)
  	CHKERR(errcode, "Failed to enqueue read buffer!");
	
	
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
    clReleaseKernel(clKernel_nw1);
    clReleaseKernel(clKernel_nw2);
    clReleaseProgram(clProgram);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);

}

