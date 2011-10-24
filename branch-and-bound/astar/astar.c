#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "../../include/rdtsc.h"

#define CHKERR(err, str) \
    if (err != CL_SUCCESS) \
    { \
        fprintf(stderr, "CL Error %d: %s\n", err, str); \
        exit(1); \
    }



//#define USEGPU 1
#define CITIES 14

int platform_id = PLATFORM_ID, n_device = DEVICE_ID;
void CPUsearch(
        int *h,
        int *city,
        int start,
        int end,
        int *result
        ) {
    int i, j, counter, counter2, tmp, x, y, tmp2, tentative_score, tent_is_better, pos;
    int closedSet[CITIES];
    int openSet[CITIES];
    int f_score[CITIES];
    int g_score[CITIES];
    //make openSet contain the starting node
    counter = 1;
    counter2 = 0;
    pos = -1;
    for (i = 0; i < CITIES; i++) {
        f_score[i] = 0;
        g_score[i] = 0;
    }

    openSet[0] = start;
    g_score[0] = 0;
    f_score[0] = h[0];
    while (counter != 0) {
        tmp = 0;
        tmp2 = f_score[openSet[0]];
        for (i = 0; i < counter; i++) {
            if (tmp2 > f_score[openSet[i]]) {
                tmp = i;
                tmp2 = f_score[openSet[i]];
            }
        }

        x = openSet[tmp];
        if (x == end) {
            result[0] = g_score[x];
            return;
        }

        closedSet[counter2] = x;
        counter2++;
        openSet[tmp] = openSet[counter - 1];
        counter--;
        for (i = 0; i < CITIES; i++) {
            if (i == x) {
                continue;
            }
            if (city[x * CITIES + i] == -1) {
                continue;
            }
            y = i;
            tmp2 = 0;

            for (j = 0; j < counter2; j++) {
                if (y == closedSet[j]) {
                    tmp2 = 1;
                }
            }
            if (tmp2) {
                continue;
            }
            tentative_score = g_score[x] + city[x * CITIES + y];
            tent_is_better = 0;
            tmp2 = 0;

            for (j = 0; j < counter; j++) {
                if (y == openSet[j]) {
                    tmp2 = 1;
                }
            }

            if (tmp2 == 0) {
                openSet[counter] = i;
                counter++;
                tent_is_better = 1;
            } else if (tentative_score < g_score[y]) {
                tent_is_better = 1;
            } else {
                tent_is_better = 0;
            }

            if (tent_is_better == 1) {
                pos++;
                g_score[y] = tentative_score;
                f_score[y] = g_score[y] + h[end * CITIES + y];
            }
        }

    }
    result[0] = -1;
    return;
}

int main(int argc, char** argv) {
    INI_TIMER
	cl_int err;
    int i, j, k;

    unsigned int correct;

    size_t global_size;
    size_t local_size;

    cl_device_id device_id;
    cl_context context;
    cl_command_queue commands;
    cl_program program;
    cl_kernel kernel;

    cl_mem h_mem;
    cl_mem city_mem, result_mem, traverse_mem;

    FILE *kernelFile;
    char *kernelSource;
    size_t kernelLength;
    size_t lengthRead;
    int start = 0, end = 1, result[CITIES * CITIES], traverse[CITIES * CITIES * CITIES], CPU_result[1];
	if(argc == 3)
	{
		platform_id = atoi(argv[1]);
		n_device = atoi(argv[2]);
	}
    /* Fill input set with distance and heuristic value */
    int h[] = {0, 366, 272, 219, 139, 229, 389, 176, 90, 269, 166, 116, 79, 36,
        366, 0, 160, 242, 244, 178, 77, 241, 380, 98, 193, 253, 329, 374,
        272, 160, 0, 87, 123, 179, 119, 93, 332, 103, 129, 186, 202, 299,
        219, 242, 87, 0, 76, 209, 212, 33, 296, 169, 134, 173, 133, 256,
        139, 244, 123, 76, 0, 156, 242, 33, 206, 156, 69, 83, 69, 169,
        229, 178, 179, 209, 156, 0, 199, 179, 219, 83, 73, 93, 216, 219,
        389, 77, 119, 212, 242, 199, 0, 219, 412, 106, 199, 259, 322, 376,
        176, 241, 93, 33, 33, 179, 219, 0, 256, 153, 89, 126, 96, 212,
        90, 380, 332, 296, 206, 219, 412, 256, 0, 298, 199, 136, 173, 39,
        269, 98, 103, 169, 156, 83, 106, 153, 298, 0, 93, 143, 232, 276,
        166, 193, 129, 134, 69, 73, 199, 89, 199, 93, 0, 46, 136, 163,
        116, 253, 186, 173, 83, 93, 259, 126, 136, 143, 46, 0, 116, 119,
        79, 329, 202, 133, 69, 216, 322, 96, 173, 232, 136, 116, 0, 129,
        36, 374, 299, 256, 169, 219, 376, 212, 39, 276, 163, 119, 129, 0};
    int city[] = {-1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 140, 118, 75,
        -1, -1, -1, -1, -1, 211, 90, -1, -1, 101, -1, -1, -1, -1,
        -1, -1, -1, 120, -1, -1, -1, -1, -1, 138, 146, -1, -1, -1,
        -1, -1, 120, -1, -1, -1, -1, 75, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, 70, -1, -1, -1, -1, 111, -1,
        -1, 211, -1, -1, -1, -1, -1, -1, -1, -1, -1, 99, -1, -1,
        -1, 90, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, 75, 70, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 151, -1, 71,
        -1, 101, 138, -1, -1, -1, -1, -1, -1, -1, 97, -1, -1, -1,
        -1, -1, 146, -1, -1, -1, -1, -1, -1, 97, -1, 80, -1, -1,
        140, -1, -1, -1, -1, 99, -1, -1, 151, -1, 80, -1, -1, -1,
        118, -1, -1, -1, 111, -1, -1, -1, -1, -1, -1, -1, -1, -1,
        75, -1, -1, -1, -1, -1, -1, -1, 71, -1, -1, -1, -1, -1};

	device_id = GetDevice(platform_id, n_device);

    /* Create a compute context */
    context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
    CHKERR(err, "Failed to create a compute context!");

    /* Create a command queue */
    commands = clCreateCommandQueue(context, device_id, TIMER_ENABLE, &err);
    CHKERR(err, "Failed to create a command queue!");

    /* Load kernel source */
    kernelFile = fopen("astar.cl", "r");
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
    kernel = clCreateKernel(program, "search", &err);
    CHKERR(err, "Failed to create a compute kernel!");

    /* Create the input and output arrays in device memory for our calculation */
    h_mem = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof (int) *CITIES*CITIES, NULL, &err);
    CHKERR(err, "Failed to allocate device memory!");
    city_mem = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof (int) *CITIES*CITIES, NULL, &err);
    CHKERR(err, "Failed to allocate device memory!");
    result_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof (int) *CITIES*CITIES, NULL, &err);
    CHKERR(err, "Failed to allocate device memory!");
    traverse_mem = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof (int) *CITIES * CITIES*CITIES, NULL, &err);
    CHKERR(err, "Failed to allocate device memory!");

    /* Write our data set into the input array in device memory */
   START_TIMER 
	err = clEnqueueWriteBuffer(commands, h_mem, CL_TRUE, 0, sizeof (int) *CITIES*CITIES, h, 0, NULL, &myEvent);
	END_TIMER
	COUNT_H2D
    CHKERR(err, "Failed to write to source array!");
    err = clEnqueueWriteBuffer(commands, city_mem, CL_TRUE, 0, sizeof (int) *CITIES*CITIES, city, 0, NULL, &myEvent);
    CHKERR(err, "Failed to write to source array!");
	CL_FINISH(commands)
	END_TIMER
	COUNT_H2D

    /* Set the arguments to our compute kernel */
    err = 0;
    err = clSetKernelArg(kernel, 0, sizeof (cl_mem), &h_mem);
    err |= clSetKernelArg(kernel, 1, sizeof (cl_mem), &city_mem);
    //err |= clSetKernelArg(kernel, 2, sizeof (int), (void *) &start);
   // err |= clSetKernelArg(kernel, 3, sizeof (int), (void *) &end);
    err |= clSetKernelArg(kernel, 2, sizeof (cl_mem), &result_mem);
    err |= clSetKernelArg(kernel, 3, sizeof (cl_mem), &traverse_mem);
    CHKERR(err, "Failed to set kernel arguments!");

    /* Get the maximum work group size for executing the kernel on the device */
    err = clGetKernelWorkGroupInfo(kernel, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof (size_t), (void *) &local_size, NULL);
    CHKERR(err, "Failed to retrieve kernel work group info!");

    /* Execute the kernel over the entire range of our 1d input data set */
    /* Get the maximum work group size for executing the kernel on the device */
    err = clGetKernelWorkGroupInfo(kernel, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof (size_t), (void *) &local_size, NULL);
    CHKERR(err, "Failed to retrieve kernel work group info!");

    /* using the maximum number of work group items for this device */
    global_size = CITIES*CITIES;
    local_size = CITIES;
	START_TIMER
    err = clEnqueueNDRangeKernel(commands, kernel, 1, NULL, &global_size, &local_size, 0, NULL, &myEvent);
    CHKERR(err, "Failed to execute kernel!");
    /* Wait for the command commands to get serviced before reading back results */
    clFinish(commands);

	END_TIMER
	COUNT_K
    /* Read back the results from the device to verify the output */
    START_TIMER
	err = clEnqueueReadBuffer(commands, result_mem, CL_TRUE, 0, sizeof (int) *CITIES*CITIES, result, 0, NULL, &myEvent);
	END_TIMER
	COUNT_D2H
    CHKERR(err, "Failed to read output array!");
    err = clEnqueueReadBuffer(commands, traverse_mem, CL_TRUE, 0, sizeof (int) *CITIES * CITIES*CITIES, traverse, 0, NULL, &myEvent);
	CL_FINISH(commands)
    CHKERR(err, "Failed to read output array!");
	END_TIMER
	COUNT_D2H
    /* Validate our results */
    for (i = 0; i < CITIES; i++)
        for (j = 0; j < CITIES; j++) {
            CPUsearch(h, city, i, j, CPU_result);
            if (CPU_result[0] != result[CITIES * i + j]) {
                printf("Validation fail");
                return -1;
            }
        }
    /* Print a brief summary detailing the results */
   #ifndef ENABLE_TIMER 
    for (i = 0; i < CITIES; i++) {
        for (k = 0; k < CITIES; k++) {
            printf("The shortest path from node %d to node %d is %d", i, k, i);
            for (j = 0; traverse[CITIES * CITIES * i + CITIES * k + j] != -1; j++)
        printf("  %d", traverse[CITIES * CITIES * i + CITIES * k + j]);
            printf("  %d\n", k);
            printf("The shortest distance from %d to node %d is %d\n", i, k, result[CITIES * i + k]);
        }
    }
#endif

    /* Shutdown and cleanup */
    clReleaseMemObject(h_mem);
    clReleaseMemObject(city_mem);
    clReleaseMemObject(result_mem);
    clReleaseMemObject(traverse_mem);
    clReleaseProgram(program);
    clReleaseKernel(kernel);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);
	PRINT_COUNT
    return 0;
}


