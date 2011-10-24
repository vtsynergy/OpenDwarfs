/*
 * File:   main.c
 * Author: ben
 *
 * Created on October 21, 2010, 8:18 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "../../include/rdtsc.h"

int platform_id = PLATFORM_ID, n_device = DEVICE_ID;

int SHARED_SIZE_LIMIT=512;

#define CHKERR(err, str) \
    if (err != CL_SUCCESS) \
    { \
        fprintf(stderr, "CL Error %d: %s\n", err, str); \
        exit(1); \
    }
//#define USEGPU 1

int compare(const void * a, const void * b) {
    return ( *(int*) a - *(int*) b);
}

int main(int argc, char** argv) {
    INI_TIMER
	cl_int err;
    unsigned int *h_InputKey, *h_InputVal, *h_OutputKeyGPU, *h_OutputValGPU;
    cl_mem d_InputKey, d_InputVal, d_OutputKey, d_OutputVal;
    unsigned int hTimer, batchSize;
    unsigned int i, j, factorizationRemainder, blockCount, threadCount, arrayLength = 64U, glo = 128;
    int flag = 1, dir = 1;
    unsigned int log2L;

    size_t local_size, global_size, local;

    cl_device_id device_id;
    cl_context context;
    cl_command_queue commands;
    cl_program program;
    cl_kernel kernel;

    if(argc == 3)
	{
		platform_id = atoi(argv[1]);
		n_device = atoi(argv[2]);
	}
    FILE *kernelFile;
    char *kernelSource;
    size_t kernelLength;
    printf(" Starting...\n\n");

    unsigned int N = 512;
    const unsigned int DIR = 0;
    const unsigned int numValues = 65536;
    const unsigned int numIterations = 1;
    printf("Allocating and initializing host arrays...\n\n");
    h_InputKey = (unsigned int *) malloc(N * sizeof (unsigned int));
    h_InputVal = (unsigned int *) malloc(N * sizeof (unsigned int));
    h_OutputKeyGPU = (unsigned int *) malloc(N * sizeof (unsigned int));
    h_OutputValGPU = (unsigned int *) malloc(N * sizeof (unsigned int));
    srand(2001);
    for (i = 0; i < N; i++) {
        h_InputKey[i] = rand() % numValues;
        h_InputVal[i] = i;
    }
	#ifndef ENABLE_TIMER
    for (i = 0; i < N; i++)
        printf("%u  ", h_InputKey[i]);
	#endif
    printf("Running GPU oddEvenMergeSort...\n\n");

    dir = (dir != 0);
//    batchSize = N / arrayLength;
//    blockCount = /*batchSize * arrayLength*/N / SHARED_SIZE_LIMIT;
//    threadCount = SHARED_SIZE_LIMIT / 2;
//    //assert((batchSize * arrayLength) % SHARED_SIZE_LIMIT == 0);
//    glo = threadCount*blockCount;
//    global_size = glo;
//    local = threadCount;

    //glo/=2;

device_id = GetDevice(platform_id, n_device);

    size_t max_worksize[3];
    err = clGetDeviceInfo(device_id, CL_DEVICE_MAX_WORK_ITEM_SIZES,sizeof(size_t)*3, &max_worksize, NULL);
 while(SHARED_SIZE_LIMIT/2>max_worksize[0])
	 {
		SHARED_SIZE_LIMIT = SHARED_SIZE_LIMIT/2;
    		N = N/2;
	 }

    context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
    CHKERR(err, "Failed to create a compute context!");


    commands = clCreateCommandQueue(context, device_id, TIMER_ENABLE, &err);
    CHKERR(err, "Failed to create a command queue!");


    kernelFile = fopen("OddEvenSort.cl", "r");
    fseek(kernelFile, 0, SEEK_END);
    kernelLength = (size_t) ftell(kernelFile);
    kernelSource = (char *) malloc(sizeof (char) *kernelLength);
    rewind(kernelFile);
    fread((void *) kernelSource, kernelLength, 1, kernelFile);
    fclose(kernelFile);


    program = clCreateProgramWithSource(context, 1, (const char **) &kernelSource, &kernelLength, &err);
    CHKERR(err, "Failed to create a compute program!");

    free(kernelSource);

    char arg[100];
	sprintf(arg, "-D SHARED_SIZE_LIMIT=%d",SHARED_SIZE_LIMIT);
    err = clBuildProgram(program, 0, NULL, arg, NULL, NULL);
    if (err == CL_BUILD_PROGRAM_FAILURE) {
        char *log;
        size_t logLen;
        err = clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &logLen);
        log = (char *) malloc(sizeof (char) *logLen);
        err = clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, logLen, (void *) log, NULL);
        fprintf(stderr, "CL Error %d: Failed to build program! Log:\n%s", err, log);
        free(log);
        exit(1);
    }
    CHKERR(err, "Failed to build program!");


    kernel = clCreateKernel(program, "oddEvenMergeSortShared", &err);
    CHKERR(err, "Failed to create a compute kernel!");

    d_InputKey = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof (unsigned int) *N, NULL, &err);
    CHKERR(err, "Failed to Create device buffer!");

    d_InputVal = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof (unsigned int) *N, NULL, &err);
    CHKERR(err, "Failed to Create device buffer!");

    d_OutputKey = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof (unsigned int) *N, NULL, &err);
    CHKERR(err, "Failed to Create device buffer!");

    d_OutputVal = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof (unsigned int) *N, NULL, &err);
    CHKERR(err, "Failed to Create  device buffer!");
	START_TIMER
    err = clEnqueueWriteBuffer(commands, d_InputKey, CL_TRUE, 0, sizeof (unsigned int) *N, (void *) h_InputKey, 0, NULL, &myEvent);
    END_TIMER
	COUNT_H2D
    CHKERR(err, "Failed to Enqueue  device buffer!");

    err = clEnqueueWriteBuffer(commands, d_InputVal, CL_TRUE, 0, sizeof (unsigned int) *N, (void *) h_InputVal, 0, NULL, &myEvent);
	CL_FINISH(commands)
    CHKERR(err, "Failed to Enqueue  device buffer!");
    END_TIMER
	COUNT_H2D
	for (arrayLength = 64; arrayLength <= N; arrayLength *= 2) {
        batchSize = N / arrayLength;
        blockCount = /*batchSize * arrayLength*/N / SHARED_SIZE_LIMIT;
        threadCount = SHARED_SIZE_LIMIT / 2;
        //assert((batchSize * arrayLength) % SHARED_SIZE_LIMIT == 0);
        glo = threadCount*blockCount;
        global_size = glo;
        local = threadCount;

        err = clSetKernelArg(kernel, 0, sizeof (cl_mem), (void *) &d_OutputKey);
        err |= clSetKernelArg(kernel, 1, sizeof (cl_mem), (void *) &d_OutputVal);
        err |= clSetKernelArg(kernel, 2, sizeof (cl_mem), (void *) &d_InputKey);
        err |= clSetKernelArg(kernel, 3, sizeof (cl_mem), (void *) &d_InputVal);
        err |= clSetKernelArg(kernel, 4, sizeof (unsigned int), (void *) &arrayLength);
        err |= clSetKernelArg(kernel, 5, sizeof (unsigned int), (void *) &dir);
        CHKERR(err, "Failed to set kernel arguments!");

        err = clGetKernelWorkGroupInfo(kernel, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof (size_t), (void *) &local_size, NULL);
        CHKERR(err, "Failed to retrieve kernel work group info!");
	START_TIMER	
        err = clEnqueueNDRangeKernel(commands, kernel, 1, NULL, &global_size, &local, 0, NULL, &myEvent);
        clFinish(commands);
	END_TIMER
	COUNT_K
	CHKERR(err, "Failed to execute kernel!");
    }
	START_TIMER
    err = clEnqueueReadBuffer(commands, d_OutputKey, CL_TRUE, 0, sizeof (unsigned int) *N, h_OutputKeyGPU, 0, NULL, &myEvent);
    END_TIMER
	COUNT_D2H
    CHKERR(err, "Failed to read output array!");

    err = clEnqueueReadBuffer(commands, d_OutputVal, CL_TRUE, 0, sizeof (unsigned int) *N, h_OutputValGPU, 0, NULL, &myEvent);
	CL_FINISH(commands)
    CHKERR(err, "Failed to read output array!");
	END_TIMER
	COUNT_D2H
    qsort(h_InputKey, N, sizeof (unsigned int), compare);
   #ifndef ENABLE_TIMER
	 for (i = 0; i < N; i++)
        printf("%u  ", h_OutputKeyGPU[i]);
	#endif
    for (i = 0; i < N; i++)
        if (h_OutputKeyGPU[i] != h_InputKey[i]){
            printf("The result is invalid");
            exit(1);
        }

    printf("\n");
    free(h_OutputValGPU);
    free(h_OutputKeyGPU);
    free(h_InputVal);
    free(h_InputKey);
    clReleaseMemObject(d_InputKey);
    clReleaseMemObject(d_InputVal);
    clReleaseMemObject(d_OutputKey);
    clReleaseMemObject(d_OutputVal);
    clReleaseProgram(program);
    clReleaseKernel(kernel);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);
	PRINT_COUNT
    return (EXIT_SUCCESS);
}

