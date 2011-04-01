#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <string.h>

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
#define DATA_SIZE 32768
#define MIN(a,b) a < b ? a : b


const char *KernelSourceFile = "crc_kernel.cl";

void usage()
{
	printf("graphCreator [hsivp]\n");
	printf("h         - print this help message\n");
	printf("s <seed>  - set the seed for the vertex\n");
	printf("i <file>  - take input from file instead of randomly generating code\n");
	printf("v         - verify parallel code with serial implementation of crc\n");
	printf("p <int>   - change the last 8 bits of the crc polynomial\n");
}

unsigned char serialCrc(unsigned char* h_num, size_t size, unsigned char crc)
{
    unsigned int i;
    unsigned char num = h_num[0];
    for(i = 1; i < size + 1; i++)
    {
        unsigned char crcCalc = h_num[i];
        unsigned int k;
        if(i == size)
            crcCalc = 0;
        for(k = 0; k < 8; k++)
        {
            //If the k-th bit is 1
            if((num >> (7-k)) % 2 == 1)
            {
                num ^= crc >> (k + 1);
                crcCalc ^= crc << (7-k);
            }
        }
        num = crcCalc;
    }

    return num;
}

int main(int argc, char** argv)
{
    cl_int err;

	unsigned char* h_num;
    unsigned char* h_answer;
    unsigned char* data = NULL;
    unsigned char h_table[256];
    unsigned char crc = 0x9B;
    unsigned char finalCRC;
	unsigned int run_serial = 0;
	char* file = NULL;	
    srand(time(NULL));
		
	int c;
	while((c = getopt (argc, argv, "vs:i:p:h")) != -1)
	{
		switch(c)
		{
			case 'h':
				usage();
				exit(0);
				break;
			case 'p':
				crc = atoi(optarg);				
				break;
			case 'v':
				run_serial = 1;
				break;
			case 'i':
				file = malloc(sizeof(*file) * strlen(optarg));
				strncpy(file, optarg, sizeof(*file) * strlen(optarg));
				break;
			case 's':
				srand(atoi(optarg));
				break;
			default:
				abort();
		}	
	}
    

    size_t global_size;
    size_t local_size;

    cl_platform_id platform_id;
    cl_device_id device_id;
    cl_context context;
    cl_command_queue commands;
    cl_program program;
    cl_kernel kernel_precompute;
    cl_kernel kernel_compute;
    cl_kernel kernel_reduce;

    cl_mem dev_input;
    cl_mem dev_table;
    cl_mem dev_output;
	
	unsigned int count = DATA_SIZE;
	
    //Initialize input
	if(file == NULL)
	{
		h_num = malloc(sizeof(*h_num) * DATA_SIZE);
		h_answer = malloc(sizeof(*h_num) * DATA_SIZE);
	    int i;
	    for(i = 0; i < count; i++)
	        h_num[i] = rand();
	}
	else
	{
		FILE* fp = fopen(file, "rb");
		printf("%s\n", file);
		if(!fp)
		{
			printf("Error reading file\n");
			exit(1);
		}
		h_num = malloc(sizeof(*h_num) * DATA_SIZE + 256);
		h_answer = malloc(sizeof(*h_num) * DATA_SIZE);
		size_t read = fread(h_num + 256, 1, DATA_SIZE, fp);
		printf("%zd\n", read);
		size_t pad = 256 - read % 256;
		count = read + pad;
		data = h_num;
		h_num = h_num + pad;
	}
	
    // Retrieve an OpenCL platform
    err = clGetPlatformIDs(1, &platform_id, NULL);
    CHKERR(err, "Failed to get a platform!");

    // Connect to a compute device
    err = clGetDeviceIDs(platform_id, USEGPU ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU, 1, &device_id, NULL);
    CHKERR(err, "Failed to create a device group!");

    // Create a compute context
    context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
    CHKERR(err, "Failed to create a compute context!");

    // Create a command queue
    commands = clCreateCommandQueue(context, device_id, 0, &err);
    CHKERR(err, "Failed to create a command queue!");

    // Get program source.
	FILE* kernelFile = NULL;
	kernelFile = fopen(KernelSourceFile, "r");
	if(!kernelFile)
		printf("Error reading file.\n"), exit(0);
	fseek(kernelFile, 0, SEEK_END);
	size_t kernelLength = (size_t) ftell(kernelFile);
	char* kernelSource = (char *) malloc(sizeof(char)*kernelLength);
	rewind(kernelFile);
	fread((void *) kernelSource, kernelLength, 1, kernelFile);
	fclose(kernelFile);
    
	// Create the compute program from the source buffer
    program = clCreateProgramWithSource(context, 1, (const char **) &kernelSource, NULL, &err);
    CHKERR(err, "Failed to create a compute program!");

    free(kernelSource);

    // Build the program executable
    err = clBuildProgram(program, 0, NULL, NULL, NULL, NULL);
    if (err == CL_BUILD_PROGRAM_FAILURE)
    {
        char *log;
        size_t logLen;
        err = clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &logLen);
        log = (char *) malloc(sizeof(char)*logLen);
        err = clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, logLen, (void *) log, NULL);
        fprintf(stderr, "CL Error %d: Failed to build program! Log:\n%s", err, log);
        free(log);
        exit(1);
    }
    CHKERR(err, "Failed to build program!");

    // Create the pre-compute kernel in the program we wish to run
    kernel_precompute = clCreateKernel(program, "precompute", &err);
    CHKERR(err, "Failed to create a compute kernel!");

    // Create the compute kernel in the program we wish to run
    kernel_compute = clCreateKernel(program, "compute", &err);
    CHKERR(err, "Failed to create a compute kernel!");

    // Create the reduce kernel in the program we wish to run
    kernel_reduce = clCreateKernel(program, "reduce", &err);
    CHKERR(err, "Failed to create a reduce kernel!");

    // Create the input and output arrays in device memory for our calculation
    dev_input = clCreateBuffer(context, CL_MEM_READ_ONLY, sizeof(char)*count, NULL, &err);
    CHKERR(err, "Failed to allocate device memory!");
    dev_table = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(char)*256, NULL, &err);
    CHKERR(err, "Failed to allocate device memory!");
    dev_output = clCreateBuffer(context, CL_MEM_WRITE_ONLY, sizeof(char)*count, NULL, &err);
    CHKERR(err, "Failed to allocate device memory!");

    // Write our data set into the input array in device memory
    err = clEnqueueWriteBuffer(commands, dev_input, CL_TRUE, 0, sizeof(char)*count, h_num, 0, NULL, NULL);
    CHKERR(err, "Failed to write to source array!");

    // Set the arguments to our precompute kernel
    err = 0;
    err = clSetKernelArg(kernel_precompute, 0, sizeof(cl_mem), &dev_table);
    err |= clSetKernelArg(kernel_precompute, 1, sizeof(unsigned char), &crc);
    CHKERR(err, "Failed to set precompute kernel arguments!");

    // Execute the kernel over the entire range of our 1d input data set
    // using the maximum number of work group items for this device
    global_size = 256;
    local_size = 256;
    err = clEnqueueNDRangeKernel(commands, kernel_precompute, 1, NULL, &global_size, &local_size, 0, NULL, NULL);
    CHKERR(err, "Failed to execute precompute kernel!");

    // Set the arguments to our compute kernel
    err = 0;
    err = clSetKernelArg(kernel_compute, 0, sizeof(cl_mem), &dev_input);
    err |= clSetKernelArg(kernel_compute, 1, sizeof(cl_mem), &dev_table);
    err |= clSetKernelArg(kernel_compute, 2, sizeof(cl_mem), &dev_output);
    err |= clSetKernelArg(kernel_compute, 3, sizeof(unsigned int), &count);
    CHKERR(err, "Failed to set compute kernel arguments!");

    // Get the maximum work group size for executing the kernel on the device
    err = clGetKernelWorkGroupInfo(kernel_compute, device_id, CL_KERNEL_WORK_GROUP_SIZE, sizeof(size_t), (void *) &local_size, NULL);
    CHKERR(err, "Failed to retrieve kernel_compute work group info!");

    // Wait for the command commands to get serviced before reading back results
    clFinish(commands);

    // Read back the results from the device to verify the output
    err = clEnqueueReadBuffer(commands, dev_table, CL_TRUE, 0, sizeof(char)*256, h_table, 0, NULL, NULL);
    CHKERR(err, "Failed to read output array!");

    // Execute the kernel over the entire range of our 1d input data set
    // using the maximum number of work group items for this device
    global_size = count;
    local_size = MIN(local_size, count);
    err = clEnqueueNDRangeKernel(commands, kernel_compute, 1, NULL, &global_size, &local_size, 0, NULL, NULL);
    CHKERR(err, "Failed to execute compute kernel!");

    // Wait for the command commands to get serviced before reading back results
    clFinish(commands);

    // Read back the results from the device to verify the output
    err = clEnqueueReadBuffer(commands, dev_output, CL_TRUE, 0, sizeof(char)*count, h_answer, 0, NULL, NULL);
    CHKERR(err, "Failed to read output array!");

    // Reduce our results
    finalCRC = 0;
    int i;
	for (i = 0; i < count; i++)
    {
        finalCRC ^= h_answer[i];
    }

    // Print a brief summary detailing the results
    printf("GPU Computation: '%X'\n", (int)finalCRC);
    
	// Calculate the result if done in serial to verify that we have the correct answer.
	if(run_serial)
	{
    	unsigned char serialCRC = serialCrc(h_num, count, crc);
    	printf("Serial Computation: '%X'\n", serialCRC);
	}	

    // Shutdown and cleanup
    clReleaseMemObject(dev_input);
    clReleaseMemObject(dev_table);
    clReleaseMemObject(dev_output);
    clReleaseProgram(program);
    clReleaseKernel(kernel_precompute);
    clReleaseKernel(kernel_compute);
    clReleaseKernel(kernel_reduce);
    clReleaseCommandQueue(commands);
    clReleaseContext(context);

    return 0;
}
