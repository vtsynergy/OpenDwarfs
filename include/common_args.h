#ifndef __COMMON_ARGS_H__
#define __COMMON_ARGS_H__

#ifdef __cplusplus
extern "C" {
#endif


#include <opts/opts.h>
#include <CL/cl.h>
#include <string.h>
#include "rdtsc.h"
#include <stdio.h>

#define USEGPU 1

typedef struct ocd_options
{
	int platform_id;
	int device_id;
	int use_cpu;
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

extern void _ocd_create_arguments();
extern ocd_options ocd_get_options();
extern int ocd_parse(int* argc, char*** argv);
extern cl_device_id _ocd_get_device(int platform, int device);
extern int ocd_check_requirements(ocd_requirements* reqs);
extern void _ocd_expand_list();
extern void _ocd_add_arg(option o, int size);
extern int ocd_register_arg(int type, char abbr, char* name, char* desc, void* value, optsverify verify, optssettor settor);
extern void ocd_usage();
extern void ocd_init(int* argc, char*** argv, ocd_requirements* reqs);
extern void ocd_finalize();

#ifdef __cplusplus
}
#endif

#endif //__COMMONARGS_H__
