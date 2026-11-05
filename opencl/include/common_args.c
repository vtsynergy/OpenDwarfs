#include "common_args.h"

ocd_options _settings = {0, -1, 0, 0};
ocd_requirements _requirements = {0,0,0};
option* _options = NULL;

int _options_length = 0;
int _options_size = 0;

int n_platform;
int n_device;
int optimized;  //If otimized kernel shall be used
cl_device_id device_id;
cl_context context;
cl_command_queue commands;

int _deviceType=0;//default for CPU, if no device option given

void _ocd_create_arguments()
{
	free(_options);
	_options = (option*)malloc(sizeof(option) * 6);
	option ops[5] = {{OTYPE_INT, 'p', (char*)"platform", (char*)"OpenCL Platform ID",
                     OFLAG_NONE, &_settings.platform_id, NULL, NULL, NULL, NULL},
		{OTYPE_INT, 'd', (char*)"device", (char*)"OpenCL Device ID",
                     OFLAG_NONE, &_settings.device_id, NULL, NULL, NULL, NULL},
        {OTYPE_INT, 't', (char*)"device type", (char*)"OpenCL Device type",
                     OFLAG_NONE, &_settings.device_type, NULL, NULL, NULL, NULL},
        {OTYPE_BOL, 'o', (char*)"Optimized Kernel", (char*)"Use Optimized kernel for the given platform",
                     OFLAG_SET, &_settings.optimized, NULL, NULL, NULL, NULL},
		{OTYPE_END, '\0', (char*)"", NULL,
                     OFLAG_NONE, NULL, NULL, NULL, NULL, NULL}};
	
	_options[0] = ops[0];
	_options[1] = ops[1];
	_options[2] = ops[2];
	_options[3] = ops[3];
	_options[4] = ops[4];
	_options_length = 6; // why?
	_options_size = 5;
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

cl_device_id _ocd_get_device(int platform, int device, cl_int dev_type)
{
    cl_int err;
    cl_uint nPlatforms = 1;
    char DeviceName[100];
    cl_device_id* devices;
    err = clGetPlatformIDs(0, NULL, &nPlatforms);
    CHECK_ERROR(err);

    if (nPlatforms <= 0) {
        printf("No OpenCL platforms found. Exiting.\n");
        exit(0);
    }
    if (platform < 0 || platform >= nPlatforms) // platform ID out of range
    {
        printf("Platform index %d is out of range. \n", platform);
        exit(-4);
    }
    cl_platform_id *platforms = (cl_platform_id *) malloc(sizeof (cl_platform_id) * nPlatforms);
    err = clGetPlatformIDs(nPlatforms, platforms, NULL);
    CHECK_ERROR(err);

    cl_uint nDevices = 1;
    char platformName[100];
    err = clGetPlatformInfo(platforms[platform], CL_PLATFORM_VENDOR, sizeof (platformName), platformName, NULL);
    CHECK_ERROR(err);
    printf("Platform Chosen : %s\n", platformName);


	//IF given device ID, use this, and disregard -t parameter if given
	if(device!=-1){
		err = clGetDeviceIDs(platforms[platform], CL_DEVICE_TYPE_ALL, 0, NULL, &nDevices);
		printf("Number of available devices: %d\n", nDevices);
    	if (nDevices <= 0) {
        	printf("No OpenCL Device found. Exiting.\n");
        	exit(0);
    	}
		if (device < 0 || device >= nDevices) // platform ID out of range
    	{
        	printf("Device index %d is out of range. \n", device);
        	exit(-4);
    	}
    	devices = (cl_device_id *) malloc(sizeof (cl_device_id) * nDevices);
		err = clGetDeviceIDs(platforms[platform], CL_DEVICE_TYPE_ALL, nDevices, devices, NULL);
    	err = clGetDeviceInfo(devices[device], CL_DEVICE_NAME, sizeof (DeviceName), DeviceName, NULL);
    	CHECK_ERROR(err);
	}
	//OTHERWISE, check at the device type parameter
	else{
		// query devices
		err = clGetDeviceIDs(platforms[platform], dev_type, 0, NULL, &nDevices);
		if(err == CL_DEVICE_NOT_FOUND)
		{
			fprintf(stderr,"No supported device of requested type found. Falling back to CPU.\n");
			dev_type = CL_DEVICE_TYPE_CPU;
			err = clGetDeviceIDs(platforms[platform], dev_type, 0, NULL, &nDevices);
			if(err == CL_DEVICE_NOT_FOUND){
				fprintf(stderr, "No CPU device available in this platform. Please, check your available OpenCL devices.\n"); 
				exit(-4);
			}
		}
		CHECK_ERROR(err);
		printf("Number of available devices: %d\n", nDevices);
    	if (nDevices <= 0) {
        	printf("No OpenCL Device found. Exiting.\n");
        	exit(0);
    	}
		//if (device < 0 || device >= nDevices) // platform ID out of range
    	//{
       	//	printf("Device index %d is out of range. \n", device);
        //	exit(-4);
    	//}
    	devices = (cl_device_id *) malloc(sizeof (cl_device_id) * nDevices);
    	err = clGetDeviceIDs(platforms[platform], dev_type, nDevices, devices, NULL);
    	//Get the first available device of requested type
    	err = clGetDeviceInfo(devices[0], CL_DEVICE_NAME, sizeof (DeviceName), DeviceName, NULL);
    	device=0;
    	CHECK_ERROR(err);	
	}
	    
    //Return
    printf("Device Chosen : %s\n", DeviceName);
    return devices[device];
}

int ocd_check_requirements(ocd_requirements* reqs)
{
	if(reqs == NULL)
		return 1;

	int pass = 1;

	ocd_options opts = ocd_get_options();
	cl_device_id d_id = _ocd_get_device(opts.platform_id, opts.device_id, opts.device_type);

	cl_ulong local_mem;
	clGetDeviceInfo(d_id, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(cl_ulong), &local_mem, NULL);
	if(local_mem < reqs->local_mem_size)
		pass = 0;
	reqs->local_mem_size = local_mem;
	
	cl_ulong global_mem;
	clGetDeviceInfo(d_id, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(cl_ulong), &global_mem, NULL);
	if(global_mem < reqs->global_mem_size)
		pass = 0;
	reqs->global_mem_size = global_mem;
	
	size_t workgroup_size;
	clGetDeviceInfo(d_id, CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(size_t), &workgroup_size, NULL);
	if(workgroup_size < reqs->workgroup_size)
		pass = 0;
	reqs->workgroup_size = workgroup_size;

	return pass;
}
	

void _ocd_expand_list()
{
	option* new_options = (option*)malloc(sizeof(option) * (_options_length + 5));
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

void ocd_init(int* argc, char*** argv, ocd_requirements* reqs)
{
	ocd_parse(argc, argv);
	ocd_check_requirements(reqs);
	#ifdef ENABLE_TIMER
	TIMER_INIT;
	#endif
}

void ocd_finalize()
{
	#ifdef ENABLE_TIMER
	//TIMER_FINISH;//TIMER_FINISH was broken into TIMER_STOP, TIMER_PRINT, TIMER_DEST, per nz-ocl to allow for extra functionalities
	TIMER_STOP
	TIMER_PRINT
	TIMER_DEST
	#endif
}

void checkDeviceChoice(int devType)
{

	if ( devType == 0)
		printf("CPU was selected\n");
	else if ( devType == 1)
		printf("GPU was selected\n");
	else if ( devType == 2)
		printf("MIC was selected\n");
	else if ( devType == 3)
		printf("FPGA was selected\n");
	else
		printf("Selection based on platform and device ID\n");

}

void ocd_initCL()
{

	cl_int err,dev_type;
	
	ocd_options opts = ocd_get_options();
    n_platform = opts.platform_id;
    n_device = opts.device_id;
    optimized = opts.optimized;
    _deviceType = opts.device_type;
    printf("Opt.optimized = %0d \n", opts.optimized);

	if (_deviceType == 1){
		 dev_type = CL_DEVICE_TYPE_GPU;
	}
	else if (_deviceType == 2){
		 dev_type = CL_DEVICE_TYPE_ACCELERATOR;
	}
	else if (_deviceType == 3){
		dev_type = CL_DEVICE_TYPE_ACCELERATOR;
	}
	else{
		dev_type = CL_DEVICE_TYPE_CPU;
	}

	checkDeviceChoice(_deviceType);//for debugging

	device_id = _ocd_get_device(n_platform, n_device, dev_type);
	
    // Create a compute context
    context = clCreateContext(0, 1, &device_id, NULL, NULL, &err);
    CHKERR(err, "Failed to create a compute context!");

    // Create a command queue
    commands = clCreateCommandQueue(context, device_id, CL_QUEUE_PROFILING_ENABLE, &err);
    CHKERR(err, "Failed to create a command queue!");
    
}

//Below, from nz-ocl

void check(int b,const char* msg)
{
	if(!b)
	{
		fprintf(stderr,"error: %s\n\n",msg);
		exit(-1);
	}
}

cl_program ocdBuildProgramFromFile(cl_context   context
                                 , cl_device_id device_id
                                 , const  char* kernel_file_name
                                 , const  char* args)
{
	cl_int err;
	cl_program program;
	size_t kernelLength;
	char* kernelSource;
	FILE* kernel_fp;
	size_t items_read;
	const char* kernel_file_mode;
    char kernel_file[80];
    char define_argument[100];
    cl_device_type device_type; //Reconfirming that it was correct

    clGetDeviceInfo(device_id
                   ,CL_DEVICE_TYPE
                   ,sizeof(cl_device_type)
                   ,&device_type
                   ,NULL); 

    printf("Device type is ");
    strcpy(kernel_file, kernel_file_name);
    strcpy(define_argument, " -DOPENCL -I. ");

    if(optimized)
        strcat(kernel_file,"_opt");

	kernel_file_mode = "r";
    if(device_type == CL_DEVICE_TYPE_CPU) {
        printf("CPU\n");
        strcat(kernel_file,optimized ? "_cpu.cl" : ".cl");
        if(args != NULL) strcat(define_argument,args);
    }
    else if(device_type == CL_DEVICE_TYPE_GPU) {
        printf("GPU\n");
        strcat(kernel_file,optimized ? "_gpu.cl" : ".cl");
        if(args != NULL) strcat(define_argument,args);
    }
    else if(device_type == CL_DEVICE_TYPE_ACCELERATOR) {
        printf("FPGA\n");
        strcat(kernel_file,optimized ? "_fpga.aocx" : ".aocx");
		kernel_file_mode = "rb";
    }
    else { 
        printf("Default\n");
        strcat(kernel_file,".cl");
        if(args != NULL) strcat(define_argument,args);
    }

    printf("Kernel file %s: Defines = %s\n", kernel_file, define_argument);

	kernel_fp = fopen(kernel_file, kernel_file_mode);
	if(kernel_fp == NULL){
		fprintf(stderr,"common_ocl.ocdBuildProgramFromFile() - Cannot open kernel file! %s", kernel_file);
		exit(-1);
	}
	fseek(kernel_fp, 0, SEEK_END);
	kernelLength = (size_t) ftell(kernel_fp);
	kernelSource = malloc(sizeof(char)*kernelLength);
	if(kernelSource == NULL){
		fprintf(stderr,"common_ocl.ocdBuildProgramFromFile() - Heap Overflow! Cannot allocate space for kernelSource.");
		exit(-1);
	}
	rewind(kernel_fp);
	items_read = fread((void *) kernelSource, kernelLength, 1, kernel_fp);
	if(items_read != 1){
		fprintf(stderr,"common_ocl.ocdBuildProgramFromFile() - Error reading from kernelFile");
		exit(-1);
	}
	fclose(kernel_fp);

	/* Create the compute program from the source buffer */
	if (_deviceType == 3) { //use Altera FPGA
        printf("Creating program from binary %s\n", kernel_file);
		program = clCreateProgramWithBinary(context,1,&device_id,&kernelLength,(const unsigned char**)&kernelSource,NULL,&err);
    }
	else //CPU or GPU or MIC
		program = clCreateProgramWithSource(context, 1, (const char **) &kernelSource, &kernelLength, &err);
	CHKERR(err, "common_ocl.ocdBuildProgramFromFile() - Failed to create a compute program!");

	/* Build the program executable */
	if (_deviceType == 3) //use Altera FPGA
		err = clBuildProgram(program,1,&device_id,define_argument,NULL,NULL);
	else
		err = clBuildProgram(program, 0, NULL,define_argument, NULL, NULL);
	
	if (err == CL_BUILD_PROGRAM_FAILURE)
	{
		char *buildLog;
		size_t logLen;
		err = clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, 0, NULL, &logLen);
		buildLog = (char *) malloc(sizeof(char)*logLen);
		if(buildLog == NULL){
			fprintf(stderr,"common_ocl.ocdBuildProgramFromFile() - Heap Overflow! Cannot allocate space for buildLog.");
			exit(-1);
		}
		err = clGetProgramBuildInfo(program, device_id, CL_PROGRAM_BUILD_LOG, logLen, (void *) buildLog, NULL);
		fprintf(stderr, "CL Error %d: Failed to build program! Log:\n%s", err, buildLog);
		free(buildLog);
		exit(1);
	}
	CHKERR(err,"common_ocl.ocdBuildProgramFromFile() - Failed to build program!");

	free(kernelSource); /* Free kernel source */
	return program;
}


void* char_new_array(const size_t N,const char* error_msg)
{
	void* ptr;
	int err;
	if (_deviceType == 3){
		err = posix_memalign(&ptr,ACL_ALIGNMENT,N * sizeof(char));
		check(err == 0,error_msg);
	}
	else{
		ptr = malloc(N * sizeof(char));
		check(ptr != NULL,error_msg);
	}
	return ptr;
}

void* int_new_array(const size_t N,const char* error_msg)
{
	void* ptr;
	int err;
	if (_deviceType == 3){
		err = posix_memalign(&ptr,ACL_ALIGNMENT,N * sizeof(int));
		check(err == 0,error_msg);
	}
	else{
		ptr = malloc(N * sizeof(int));
		check(ptr != NULL,error_msg);
	}
	return ptr;
}

void* long_new_array(const size_t N,const char* error_msg)
{
	void* ptr;
	int err;
	if (_deviceType == 3){
		err = posix_memalign(&ptr,ACL_ALIGNMENT,N * sizeof(long));
		check(err == 0,error_msg);
	}
	else{
		ptr = malloc(N * sizeof(long));
		check(ptr != NULL,error_msg);
	}
	return ptr;
}

void* float_new_array(const size_t N,const char* error_msg)
{
	void* ptr;
	int err;
	if (_deviceType == 3){
		err = posix_memalign(&ptr,ACL_ALIGNMENT,N * sizeof(float));
		check(!err,error_msg);
	}
	else{
		ptr = malloc(N * sizeof(float));
		check(ptr != NULL,error_msg);
	}
	return ptr;
}

void* float_array_realloc(void* ptr,const size_t N,const char* error_msg)
{
	int err;
	if (_deviceType == 3){
		if(ptr != NULL) free(ptr);
		err = posix_memalign(&ptr,ACL_ALIGNMENT,N * sizeof(float));
		check(!err,error_msg);
	}
	else{
		ptr = realloc(ptr,N * sizeof(float));
		check(ptr != NULL,error_msg);
	}
	return ptr;
}
