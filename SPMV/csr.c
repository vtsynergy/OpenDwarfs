#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

#define CHKERR(err, str) \
    if (err != CL_SUCCESS) \
    { \
        fprintf(stderr, "CL Error %d: %s\n", err, str); \
        exit(1); \
    }

#define USEGPU 1
#include "common.h"
#include "common.c"
#include "sparse_formats.h"

int main(int argc, char** argv)
{
    cl_int err;

    unsigned int correct;

    size_t global_size;
    size_t local_size;

    cl_platform_id platform_id;
    cl_device_id device_id;
    cl_context context;
    cl_command_queue commands;
    cl_program program;
    cl_kernel kernel;

    stopwatch sw;

    cl_mem csr_ap;
    cl_mem csr_aj;
    cl_mem csr_ax;
    cl_mem x_loc;
    cl_mem y_loc;

    FILE *kernelFile;
    char *kernelSource;
    size_t kernelLength;
    size_t lengthRead;

    /* Fill input set with random float values */
    int i;

    csr_matrix csr;
    csr = laplacian_5pt(512);
    int k = 0;
      for(k = 0; k < csr.num_nonzeros+1; k++){
         csr.Ax[k] = 1.0 - 2.0 * (rand() / (RAND_MAX + 1.0));
      }

    //The other arrays
    float * x_host = float_new_array(csr.num_cols+1);
    float * y_host = float_new_array(csr.num_rows+1);

    unsigned int ii;
    for(ii = 0; ii < csr.num_cols+1; ii++){
        x_host[ii] = rand() / (RAND_MAX + 1.0);
    }
    for(ii = 0; ii < csr.num_rows+1; ii++){
        y_host[ii] = rand() / (RAND_MAX + 2.0);
    }

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
    kernelFile = fopen("spmv_csr_kernel.cl", "r");
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
    kernel = clCreateKernel(program, "csr", &err);
    CHKERR(err, "Failed to create a compute kernel!");

    /* Create the input and output arrays in device memory for our calculation */
    csr_ap = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned int)*csr.num_rows+1, NULL, &err);
    CHKERR(err, "Failed to allocate device memory!");
    csr_aj = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(unsigned int)*csr.num_nonzeros+1, NULL, &err);
    CHKERR(err, "Failed to allocate device memory!");
    csr_ax = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float)*csr.num_nonzeros+1, NULL, &err);
    CHKERR(err, "Failed to allocate device memory!");
    x_loc = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float)*csr.num_cols+1, NULL, &err);
    CHKERR(err, "Failed to allocate device memory!");
    y_loc = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(float)*csr.num_rows+1, NULL, &err);
    CHKERR(err, "Failed to allocate device memory!");

    /* beginning of timing point */
    stopwatch_start(&sw); 
   
    /* Write our data set into the input array in device memory */
    err = clEnqueueWriteBuffer(commands, csr_ap, CL_TRUE, 0, sizeof(unsigned int)*csr.num_rows+1, csr.Ap, 0, NULL, NULL);
    CHKERR(err, "Failed to write to source array!");
    err = clEnqueueWriteBuffer(commands, csr_aj, CL_TRUE, 0, sizeof(unsigned int)*csr.num_nonzeros+1, csr.Aj, 0, NULL, NULL);
    CHKERR(err, "Failed to write to source array!");
    err = clEnqueueWriteBuffer(commands, csr_ax, CL_TRUE, 0, sizeof(float)*csr.num_nonzeros+1, csr.Ax, 0, NULL, NULL);
    CHKERR(err, "Failed to write to source array!");
    err = clEnqueueWriteBuffer(commands, x_loc, CL_TRUE, 0, sizeof(float)*csr.num_cols+1, x_host, 0, NULL, NULL);
    CHKERR(err, "Failed to write to source array!");
    err = clEnqueueWriteBuffer(commands, y_loc, CL_TRUE, 0, sizeof(float)*csr.num_rows+1, y_host, 0, NULL, NULL);
    CHKERR(err, "Failed to write to source array!");

    /* Set the arguments to our compute kernel */
    err = 0;
    err = clSetKernelArg(kernel, 0, sizeof(unsigned int), &csr.num_rows);
    err |= clSetKernelArg(kernel, 1, sizeof(cl_mem), &csr_ap);
    err |= clSetKernelArg(kernel, 2, sizeof(cl_mem), &csr_aj);
    err |= clSetKernelArg(kernel, 3, sizeof(cl_mem), &csr_ax);
    err |= clSetKernelArg(kernel, 4, sizeof(cl_mem), &x_loc);
    err |= clSetKernelArg(kernel, 5, sizeof(cl_mem), &y_loc);
    CHKERR(err, "Failed to set kernel arguments!");

    /* Get the maximum work group size for executing the kernel on the device */
    err = clGetKernelWorkGroupInfo(kernel, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), (void *) &local_size, NULL);
    CHKERR(err, "Failed to retrieve kernel work group info!");

    /* Execute the kernel over the entire range of our 1d input data set */
    /* using the maximum number of work group items for this device */
    global_size = csr.num_rows;
    err = clEnqueueNDRangeKernel(commands, kernel, 1, NULL, &global_size, &local_size, 0, NULL, NULL);
    CHKERR(err, "Failed to execute kernel!");

    /* Wait for the command commands to get serviced before reading back results */
    clFinish(commands);
    float output[csr.num_rows+1];
    
    /* Read back the results from the device to verify the output */
    err = clEnqueueReadBuffer(commands, y_loc, CL_TRUE, 0, sizeof(float)*csr.num_rows+1, output, 0, NULL, NULL);
    CHKERR(err, "Failed to read output array!");

    /* end of timing point */
    stopwatch_stop(&sw);
    printf("Time consumed(ms): %lf Gflops: %f \n", 1000*get_interval_by_sec(&sw), (2.0 * (double) csr.num_nonzeros / get_interval_by_sec(&sw)) / 1e9);


    /* Validate our results */
    //for (i = 0; i < csr.num_rows; i++)
    //{
    //    printf("row: %d	output: %f \n", i, output[i]);  
    //}

    /* Print a brief summary detailing the results */

    /* Shutdown and cleanup */
    clReleaseMemObject(csr_ap);
    clReleaseMemObject(csr_aj);
    clReleaseMemObject(csr_ax);
    clReleaseMemObject(x_loc);
    clReleaseMemObject(y_loc);
    clReleaseProgram(program);
    clReleaseKernel(kernel);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);

    return 0;
}
