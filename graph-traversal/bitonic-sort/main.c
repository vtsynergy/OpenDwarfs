#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <unistd.h>
#include <CL/cl.h>

#define CHKERR(err, str) \
    if (err != CL_SUCCESS) \
    { \
        fprintf(stderr, "CL Error %d: %s\n", err, str); \
        exit(1); \
    }



#define USEGPU 1
#define MAX_VALUE 65536

int compare(const void * a, const void * b) {
    return ( *(int*) a - *(int*) b);
}

void usage() {
    printf("bsort [ns]\n");
    printf("n <number> - number of items to be sorted( must be power of 2\n");
    printf("s <seed>  - set the seed for the random number\n");
}

int main(int argc, char** argv) {
    cl_int err;
    unsigned int i, j, k;
    int c;

    const unsigned int numValues = MAX_VALUE;

    size_t global_size;
    size_t local_size, max_local, i_, j_;

    cl_platform_id platform_id;
    cl_device_id device_id;
    cl_context context;
    cl_command_queue commands;
    cl_program program;
    cl_kernel kernel;
    cl_event events[2];


    cl_mem input_mem;

    FILE *kernelFile;
    char *kernelSource;
    size_t kernelLength;
    size_t lengthRead;

    unsigned int* input;
    unsigned int* output;
    unsigned int size = 512;
    int seed = 2010;
    while ((c = getopt(argc, argv, "h:n:s")) != -1) {
        switch (c) {
            case 'h':
                usage();
                exit(0);
                break;
            case 'n':
                size = atoi(optarg);
                break;
            case 's':
                seed = atoi(optarg);
                break;
            default:
                abort();
        }
    }
    if (size <= 1) {
        usage();
        abort();
    }
    for (i = size; i > 1; i /= 2) {
        if (i % 2 ==1) {
            usage();
            abort();
        }

    }

    input = malloc(sizeof (unsigned int) * size);
    output = malloc(sizeof (unsigned int) * size);
    /* Fill input set  */
    srand(seed);
    for (i = 0; i < size; i++) {
        input [i] = rand() % numValues;
        output[i] = input[i];
    }
    for (i = 0; i < size; i++)
        printf("%u  ", input[i]);
    printf("\nRunning GPU bitonic sort...\n\n");



    /* Retrieve an OpenCL platform */
    err = clGetPlatformIDs(1, &platform_id, NULL);
    CHKERR(err, "Failed to get a platform!");

    /* Connect to a compute device */
    err = clGetDeviceIDs(platform_id, USEGPU ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU, 1, &device_id, NULL);
    CHKERR(err, "Failed to create a device group!");

    /* Create a compute context */
    context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
    CHKERR(err, "Failed to create a compute context!");

    /* Create a command queue */
    commands = clCreateCommandQueue(context, device_id, 0, &err);
    CHKERR(err, "Failed to create a command queue!");

    /* Load kernel source */
    kernelFile = fopen("sort.cl", "r");
    fseek(kernelFile, 0, SEEK_END);
    kernelLength = (size_t) ftell(kernelFile);
    kernelSource = (char *) malloc(sizeof (char) *kernelLength);
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
    if (err == CL_BUILD_PROGRAM_FAILURE) {
        char *buildLog;
        size_t logLen;
        err = clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &logLen);
        buildLog = (char *) malloc(sizeof (char) *logLen);
        err = clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, logLen, (void *) buildLog, NULL);
        fprintf(stderr, "CL Error %d: Failed to build program! Log:\n%s", err, buildLog);
        free(buildLog);
        exit(1);
    }
    CHKERR(err, "Failed to build program!");

    /* Create the compute kernel in the program we wish to run */
    kernel = clCreateKernel(program, "sort", &err);
    CHKERR(err, "Failed to create a compute kernel!");


    /* Create the input and output arrays in device memory for our calculation */
    input_mem = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof (unsigned int) *size, NULL, &err);
    CHKERR(err, "Failed to allocate device memory!");


    /* Write our data set into the input array in device memory */
    err = clEnqueueWriteBuffer(commands, input_mem, CL_TRUE, 0, sizeof (unsigned int) *size, input, 0, NULL, NULL);
    CHKERR(err, "Failed to write to source array!");


    /* Execute the kernel over the entire range of our 1d input data set */
    /* Get the maximum work group size for executing the kernel on the device */
    err = clGetKernelWorkGroupInfo(kernel, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof (size_t), (void *) &max_local, NULL);
    CHKERR(err, "Failed to retrieve kernel work group info!");
    for (local_size = 1024; local_size >= 1; local_size /= 2) {
        if (max_local >= local_size) {
            break;
        }
    }
    global_size = size;
    if (global_size < local_size)
        local_size = global_size;
    /* Timer starts here*/
    /* Set the arguments to our compute kernel */
    err = clSetKernelArg(kernel, 0, sizeof (cl_mem), &input_mem);
    CHKERR(err, "Failed to set kernel arguments!");
    for (i = 2; i <= size; i *= 2) {
        for (j = i / 2; j > 0; j /= 2) {
            err = clSetKernelArg(kernel, 1, sizeof (unsigned int), (void *) &j);
            err |= clSetKernelArg(kernel, 2, sizeof (unsigned int), (void *) &i);
            err = clEnqueueNDRangeKernel(commands, kernel, 1, NULL, &global_size, &local_size, 0, NULL, &events [0]);
            CHKERR(err, "Failed to execute kernel!");
            err = clWaitForEvents(1, &events [0]);
            CHKERR(err, "Failed to execute event");

        }
    }

    /* Wait for the command commands to get serviced before reading back results */
    clFinish(commands);

    /* Read back the results from the device to verify the output */
    err = clEnqueueReadBuffer(commands, input_mem, CL_TRUE, 0, sizeof (unsigned int) * size, input, 0, NULL, NULL);
    CHKERR(err, "Failed to read output array!");
    /* timer ends here */

    /* Validate our results */
    qsort(output, size, sizeof (unsigned int), compare);
    for (i = 0; i < size; i++)
        if (input[i] != output[i]) {
            printf("The result is invalid");
            exit(1);
        }
    /* Print a brief summary detailing the results */
    for (i = 0; i < size; i++)
        printf("%u  ", input[i]);

    printf("The result is valid\n");

    /* Shutdown and cleanup */
    clReleaseMemObject(input_mem);
    clReleaseProgram(program);
    clReleaseKernel(kernel);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);

    return 0;
}
