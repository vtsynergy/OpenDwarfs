#ifndef __RDTSC_H__
#define __RDTSC_H__
#include <time.h>
#define CHECK_ERROR(err) \
	if (err != CL_SUCCESS) \
{ \
	fprintf(stderr, "Error1: %d\n", err);\
	exit(1); \
}
#ifdef __APPLE__
#include <OpenCL/opencl.h>
#else
#include <CL/opencl.h>
#endif

#define USEGPU 1
#define PLATFORM_ID 0
#define DEVICE_ID 0
#include <stdint.h>

cl_event myEvent; 
cl_ulong startTime, endTime;
unsigned long long int h2d_t=0, k_t=0, d2h_t=0;
#ifdef ENABLE_TIMER
#define START_TIMER
#define INI_TIMER printf("TIMER ENABLE\n");
#define TIMER_ENABLE CL_QUEUE_PROFILING_ENABLE
#define END_TIMER clGetEventProfilingInfo(myEvent, CL_PROFILING_COMMAND_START, sizeof(cl_ulong), &startTime, NULL);\
		clGetEventProfilingInfo(myEvent, CL_PROFILING_COMMAND_END, sizeof(cl_ulong), &endTime, NULL);
#define COUNT_H2D h2d_t += (endTime-startTime);
#define COUNT_K   k_t   += (endTime-startTime);
#define COUNT_D2H d2h_t += (endTime-startTime);

#define PRINT_COUNT printf("H2D Transfer Time: %f (ms)\n",h2d_t*1e-6f);\
		    printf("Kernel Time: %f (ms)\n",k_t*1e-6f);\
		    printf("D2H Transfer Time: %f (ms)\n",d2h_t*1e-6f);

#define CL_FINISH(c) clFinish(c);
#else
#define INI_TIMER
#define START_TIMER
#define TIMER_ENABLE 0
#define END_TIMER
#define COUNT_H2D
#define COUNT_K
#define COUNT_D2H
#define PRINT_COUNT 
#define CL_FINISH(c)
#endif

#ifdef START_POWER
#define START_KERNEL printf("Kernel Start\n");
#define END_KERNEL printf("Kernel END\n");
#else
#define START_KERNEL
#define END_KERNEL
#endif

	cl_device_id 
GetDevice(int platform, int device)
{
	cl_int err;
	cl_uint nPlatforms = 1;
	err = clGetPlatformIDs(0, NULL, &nPlatforms);
	CHECK_ERROR(err);

	if (nPlatforms <= 0)
	{
		printf("No OpenCL platforms found. Exiting.\n");
		exit(0);
	}
	if (platform<0 || platform>=nPlatforms)  // platform ID out of range
	{
		printf("Platform index %d is out of range. \n", platform);
		exit(-4);
	}
	cl_platform_id *platforms = (cl_platform_id *)malloc(sizeof(cl_platform_id)*nPlatforms);
	err = clGetPlatformIDs(nPlatforms, platforms, NULL);
	CHECK_ERROR(err);

	cl_uint nDevices = 1;
	char platformName[100];
	err = clGetPlatformInfo(platforms[0], CL_PLATFORM_VENDOR, sizeof(platformName), platformName, NULL);
	CHECK_ERROR(err);
	printf("Platform Chosen : %s\n", platformName);
	// query devices
	err = clGetDeviceIDs(platforms[platform], USEGPU ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU, 0, NULL, &nDevices);
	CHECK_ERROR(err);
	if (nDevices <= 0)
	{
		printf("No OpenCL Device found. Exiting.\n");
		exit(0);
	}
	if (device<0 || device>=nDevices)  // platform ID out of range
	{
		printf("Device index %d is out of range. \n", device);
		exit(-4);
	}
	cl_device_id* devices = (cl_device_id *) malloc (sizeof(cl_device_id) * nDevices);
	err = clGetDeviceIDs(platforms[platform], USEGPU ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU, nDevices, devices, NULL);
	char DeviceName[100];
	err = clGetDeviceInfo(devices[device], CL_DEVICE_NAME, sizeof(DeviceName), DeviceName, NULL);
	CHECK_ERROR(err);
	printf("Device Chosen : %s\n", DeviceName);

	return devices[device];
}

#endif

