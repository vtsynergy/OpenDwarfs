#ifndef __COMMON_ARGS_H__
#define __COMMON_ARGS_H__

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

#include <opts.h>

typedef struct ocd_options
{
	int platform_id;
	int device_id;
	int use_cpu;
} ocd_options;

ocd_options _settings = {0, 0, 0};

option* _options = NULL;

void _ocd_create_arguments()
{
	free(_options);
	_options = malloc(sizeof(option) * 5);
	//_options = {{OTYPE_INT, 'p', "platform", "OpenCL Platform ID",
	//		OFLAG_NONE, &_settings.platform_id, NULL, NULL, NULL, NULL}}
}

int ocd_parse(int argc, char** argv)
{
	if(!_options)
		_ocd_create_arguments();
	option* op;
	if (optsgets(argc, argv, _options)) {
		fprintf(stderr, "optsgets() errored.\nUsage:");
		for (op = &_options[0]; op->type; ++op)
			if (!((op->type == OTYPE_NUL) &&
			  (op->flags & OFLAG_ARG)))
				fprintf(stderr, "%s\n", optsusage(op));
	}
}


int ocd_register_arg(int type, char abbr, char* name, char* desc, void* value, optsverify verify, optssettor settor)
{
	option o;
	
	if(!_options)
		_ocd_create_arguments();
}

#endif
