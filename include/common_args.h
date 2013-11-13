#ifndef __COMMON_ARGS_H__
#define __COMMON_ARGS_H__

#ifdef __cplusplus
extern "C" {
#endif


#include <opts/opts.h>
#ifdef OPENCL_HEADER_CL_CL
#include <CL/cl.h>
#endif
#ifdef OPENCL_HEADER_LONG
#include <OpenCL/opencl.h>
#endif
#include <string.h>
#include "rdtsc.h"
#include <stdio.h>

#define USEGPU 1

#define MINIMUM(i,j) ((i)<(j) ? (i) : (j))
#define ACL_ALIGNMENT 64 // Minimum alignment for DMA transfer to Altera FPGA board

extern int _deviceType;

#define CHKERR(err, str) \
    if (err != CL_SUCCESS) \
    { \
        fprintf(stdout, "CL Error %d: %s\n", err, str); \
        exit(1); \
    }

typedef struct ocd_options
{
	int platform_id;
	int device_id;
	int device_type;
} ocd_options;
extern ocd_options _settings;

typedef struct ocd_requirements
{
	cl_ulong local_mem_size;
	cl_ulong global_mem_size;
	size_t workgroup_size;
} ocd_requirements;
extern ocd_requirements _requirements;

extern option* _options;
extern int _options_length;
extern int _options_size;

extern int n_platform;
extern int n_device;
extern cl_device_id device_id;
extern cl_context context;
extern cl_command_queue commands;


extern void _ocd_create_arguments();
extern ocd_options ocd_get_options();
extern int ocd_parse(int* argc, char*** argv);
extern cl_device_id _ocd_get_device(int platform, int device, cl_int dev_type);
extern int ocd_check_requirements(ocd_requirements* reqs);
extern void _ocd_expand_list();
extern void _ocd_add_arg(option o, int size);
extern int ocd_register_arg(int type, char abbr, char* name, char* desc, void* value, optsverify verify, optssettor settor);
extern void ocd_usage();
extern void ocd_init(int* argc, char*** argv, ocd_requirements* reqs);
extern void ocd_finalize();
extern void checkDeviceChoice(int devType);
extern void ocd_initCL();

//From nz-ocl
extern void check(int b,const char* msg);
extern cl_program ocdBuildProgramFromFile(cl_context context,cl_device_id device_id,const char* kernel_file_name);
extern void* char_new_array(const size_t N,const char* error_msg);
extern void* int_new_array(const size_t N,const char* error_msg);
extern void* long_new_array(const size_t N,const char* error_msg);
extern void* float_new_array(const size_t N,const char* error_msg);
extern void* float_array_realloc(void* ptr,const size_t N,const char* error_msg);


#ifdef __cplusplus
}
#endif

#endif //__COMMONARGS_H__
