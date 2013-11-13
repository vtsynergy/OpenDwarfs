#include <config.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif
#include "../include/rdtsc.h"
#include "../include/common_args.h"

#define EPSILON 0.0001

#define DATA_SIZE 1048576

int main(int argc, char** argv)
{
	ocd_init(&argc, &argv, NULL);
	cl_int err,dev_type;
    float *input = (float*)malloc(DATA_SIZE*sizeof(float));
    float *output = (float*)malloc(DATA_SIZE*sizeof(float));
    unsigned int correct;

    size_t global_size;
    size_t local_size;

    cl_platform_id platform_id;
    cl_device_id device_id;
    cl_context context;
    cl_command_queue commands;
    cl_program program;
    cl_kernel kernel;

    cl_mem dev_input;
    cl_mem dev_output;

    FILE *kernelFile;
    char *kernelSource;
    size_t kernelLength;
    size_t lengthRead;

    /* Fill input set with random float values */
    int i;
    unsigned int count = DATA_SIZE;
    for (i = 0; i < count; i++)
        input[i] = rand() / (float)RAND_MAX;

	ocd_options opts = ocd_get_options();
	int n_platform_id = opts.platform_id;
	int n_device_id = opts.device_id;


	#ifdef USEGPU
		 dev_type = CL_DEVICE_TYPE_GPU;
	#elif defined(USE_AFPGA)
		 dev_type = CL_DEVICE_TYPE_ACCELERATOR;
	#else
		dev_type = CL_DEVICE_TYPE_CPU;
	#endif

	device_id = _ocd_get_device(n_platform_id, n_device_id,dev_type);
    /* Create a compute context */
    context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
    CHKERR(err, "Failed to create a compute context!");

    /* Create a command queue */
    commands = clCreateCommandQueue(context, device_id, CL_QUEUE_PROFILING_ENABLE, &err);
    CHKERR(err, "Failed to create a command queue!");

    /* Load kernel source */
    kernelFile = fopen("samplecl_kernel.cl", "r");
    fseek(kernelFile, 0, SEEK_END);
    kernelLength = (size_t) ftell(kernelFile);
    kernelSource = (char *) malloc(sizeof(char)*kernelLength);
    rewind(kernelFile);
    lengthRead = fread((void *) kernelSource, kernelLength, 1, kernelFile);
    fclose(kernelFile);

    /* Create the compute program from the source buffer */
    program = clCreateProgramWithSource(context, 1, (const char **) &kernelSource, &kernelLength, &err);
    CHKERR(err, "Failed to create a compute program!");

    /* Free kernel source */
    free(kernelSource);

    /* Build the program executable */
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (err == CL_BUILD_PROGRAM_FAILURE)                                                                                                                                       
    {
        char *buildLog;
        size_t logLen;
        err = clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &logLen);
        buildLog = (char *) malloc(sizeof(char)*logLen);
        err = clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, logLen, (void *) buildLog, NULL);
        fprintf(stderr, "CL Error %d: Failed to build program! Log:\n%s", err, buildLog);
        free(buildLog);
        exit(1);
    }
    CHKERR(err, "Failed to build program!");

    /* Create the compute kernel in the program we wish to run */
    kernel = clCreateKernel(program, "copy", &err);
    CHKERR(err, "Failed to create a compute kernel!");

    /* Create the input and output arrays in device memory for our calculation */
    dev_input = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float)*count, NULL, &err);
    CHKERR(err, "Failed to allocate device memory!");
    dev_output = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(float)*count, NULL, &err);
    CHKERR(err, "Failed to allocate device memory!");

    cl_event writeEvent;
    /* Write our data set into the input array in device memory */
    err = clEnqueueWriteBuffer(commands, dev_input, CL_TRUE, 0, sizeof(float)*count, input, 0, NULL, &writeEvent);
    CHKERR(err, "Failed to write to source array!");
    clFinish(commands);
    struct ocdTimer * writeTimer;
    //The event must have actually started or else START_TIMER will fail
    //This is why we use clFinish
    START_TIMER(writeEvent, OCD_TIMER_H2D, NULL, writeTimer)

    /* Set the arguments to our compute kernel */
    err = 0;
    err = clSetKernelArg(kernel, 0, sizeof(cl_mem), &dev_input);
    err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &dev_output);
    err |= clSetKernelArg(kernel, 2, sizeof(unsigned int), &count);
    CHKERR(err, "Failed to set kernel arguments!");

    /* Get the maximum work group size for executing the kernel on the device */
    err = clGetKernelWorkGroupInfo(kernel, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), (void *) &local_size, NULL);
    CHKERR(err, "Failed to retrieve kernel work group info!");

    /* Execute the kernel over the entire range of our 1d input data set */
    /* using the maximum number of work group items for this device */
    cl_event kernEvent;
    global_size = count;
    err = clEnqueueNDRangeKernel(commands, kernel, 1, NULL, &global_size, &local_size, 0, NULL, &kernEvent);
    CHKERR(err, "Failed to execute kernel!");

    /* Wait for the command commands to get serviced before reading back results */
    clFinish(commands);
    struct ocdTimer * kernTimer;
    START_TIMER(kernEvent, OCD_TIMER_KERNEL, NULL, kernTimer)

    /* Read back the results from the device to verify the output */
    cl_event readEvent;
    err = clEnqueueReadBuffer(commands, dev_output, CL_TRUE, 0, sizeof(float)*count, output, 0, NULL, &readEvent);
    CHKERR(err, "Failed to read output array!");
    clFinish(commands);
    struct ocdTimer * readTimer;
    START_TIMER(readEvent, OCD_TIMER_D2H, NULL, readTimer);

    /* Validate our results */
    correct = 0;
    for (i = 0; i < count; i++)
    {
        if (fabs(output[i] - input[i]) < EPSILON)
            correct++;
    }
    //Clock in end times for each event
    START_DUAL_TIMER(writeEvent, readEvent, NULL, ocdTempDualTimer)
    //The event must have actually finished before END_TIMER
    //or else the call will fail
    END_DUAL_TIMER(ocdTempDualTimer)
    END_TIMER(writeTimer)
    END_TIMER(kernTimer)
    END_TIMER(readTimer)

    /* Print a brief summary detailing the results */
    printf("Computed '%d/%d' correct values!\n", correct, count);
    
    clReleaseEvent(writeEvent);
    clReleaseEvent(kernEvent);
    clReleaseEvent(readEvent);

    /* Shutdown and cleanup */
    clReleaseMemObject(dev_input);
    clReleaseMemObject(dev_output);
    clReleaseProgram(program);
    clReleaseKernel(kernel);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);
    
    free(input);
    free(output);
    
    //Finalize the timer suite
    ocd_finalize();

    return 0;
}
