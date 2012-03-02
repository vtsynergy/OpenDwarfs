#ifndef __COMMON_ARGS_H__
#define __COMMON_ARGS_H__

#include <opts.h>

typedef struct ocd_options
{
	int platform_id;
	int device_id;
	int use_cpu;
} ocd_options;

ocd_options _settings = {0, 0, 0};

option* _options = NULL;
int _options_length = 0;
int _options_size = 0;

void _ocd_create_arguments()
{
	free(_options);
	_options = malloc(sizeof(option) * 5);
	option ops[3] = {{OTYPE_INT, 'p', "platform", "OpenCL Platform ID",
                     OFLAG_NONE, &_settings.platform_id, NULL, NULL, NULL, NULL},
		{OTYPE_INT, 'd', "device", "OpenCL Device ID",
                     OFLAG_NONE, &_settings.device_id, NULL, NULL, NULL, NULL},
		{OTYPE_END, '\0', "", 0,
                     OFLAG_NONE, NULL, NULL, NULL, NULL, NULL}};
	
	_options[0] = ops[0];
	_options[1] = ops[1];
	_options[2] = ops[2];
	_options_length = 5;
	_options_size = 3;
}

ocd_options ocd_get_options()
{
	return _settings;
}

int ocd_parse(int* argc, char*** argv)
{
	if(!_options)
		_ocd_create_arguments();
	
	int largc = *argc;
	char** largv = *argv;

	//Determine if there are any arguments to parse.
	int i;
	int origargc = largc;
	for(i = 0; i < largc; i++)
	{
		if(strlen(largv[i]) == 2 && strncmp(largv[i], "--", 2) == 0)
		{
			largc = i;
			break;
		}
	}

	if(largc == origargc) // No double dash
	{
		//Assume we failed because they were actually
		//Program specific arguments.
		return 0;
	}

	if (optsgets(largc, largv, _options)) {

		fprintf(stderr, "optsgets() errored.\nUsage:");
		option* op;
		for (op = &_options[0]; op->type; ++op)
			if (!((op->type == OTYPE_NUL) &&
			  (op->flags & OFLAG_ARG)))
				fprintf(stderr, "%s\n", optsusage(op));
	}

	*argv += largc;
	*argc = *argc - largc;
	return largc;
}

void _ocd_expand_list()
{
	option* new_options = malloc(sizeof(option) * (_options_length + 5));
	memcpy(new_options, _options, sizeof(option) * _options_length);
	free(_options);
	_options = new_options;
	_options_length += 5;
}

void _ocd_add_arg(option o, int size)
{
	if(_options_size >= _options_length)
		_ocd_expand_list();
	
	option end = _options[size-1];
	_options[size-1] = o;
	_options[size] = end;	
	
	_options_size++;
}

int ocd_register_arg(int type, char abbr, char* name, char* desc, void* value, optsverify verify, optssettor settor)
{
	option o = {type, abbr, name, desc, OFLAG_NONE, value, 0, verify, settor, 0};
	
	if(!_options)
		_ocd_create_arguments();
	_ocd_add_arg(o, _options_size);	
}

void ocd_usage()
{
	option* op;
	for (op = &_options[0]; op->type; ++op)
		printf("%s\n", optsusage(op));
}

#endif
