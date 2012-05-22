#include <stdio.h>
#include <unistd.h>
#include <getopt.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "../../include/rdtsc.h"
#include "../../include/common_args.h"

#include "common.h"

#define CHECKERR(err) \
    if (err != CL_SUCCESS) \
    { \
        fprintf(stderr, "Error: %d\n", err);\
        exit(1); \
    }

//#define USEGPU 1
int BLOCK_SIZE = 16;
int platform_id=PLATFORM_ID, device_id=DEVICE_ID;
static int do_verify = 0;

static struct option long_options[] = {
      /* name, has_arg, flag, val */
      {"input", 1, NULL, 'i'},
      {"platform", 1, NULL, 'p'},
      {"device", 1, NULL, 'd'},
      {"size", 1, NULL, 's'},
      {"verify", 0, NULL, 'v'},
      {0,0,0,0}
};

int
main ( int argc, char *argv[] )
{
    OCD_INIT
  int matrix_dim = 32; /* default matrix_dim */
  int opt, option_index=0;
  func_ret_t ret;
  const char *input_file = NULL;
  float *m, *mm;
  stopwatch sw;

  cl_device_id clDevice;
  cl_context clContext;
  cl_command_queue clCommands;
  cl_program clProgram;
  cl_kernel clKernel_diagonal;
  cl_kernel clKernel_perimeter;
  cl_kernel clKernel_internal;

  cl_int errcode;

  FILE *kernelFile;
  char *kernelSource;
  size_t kernelLength;

  cl_mem d_m;

  ocd_init(&argc, &argv, NULL);
  ocd_options opts = ocd_get_options();
  platform_id = opts.platform_id;
  device_id = opts.device_id;


  while ((opt = getopt_long(argc, argv, "::vs:i:", 
                            long_options, &option_index)) != -1 ) {
      switch(opt){
        case 'i':
          input_file = optarg;
          break;
        case 'v':
          do_verify = 1;
          break;
        case 's':
          matrix_dim = atoi(optarg);
          fprintf(stderr, "Currently not supported, use -i instead\n");
          fprintf(stderr, "Usage: %s [-v] [-s matrix_size|-i input_file|-p platform|-d device]\n", argv[0]);
          exit(EXIT_FAILURE);
        case '?':
          fprintf(stderr, "invalid option\n");
          break;
        case ':':
          fprintf(stderr, "missing argument\n");
          break;
        default:
          fprintf(stderr, "Usage: %s [-v] [-s matrix_size|-i input_file||-p platform|-d device]\n",
                  argv[0]);
          exit(EXIT_FAILURE);
      }
  }
  
  if ( (optind < argc) || (optind == 1)) {
      fprintf(stderr, "Usage: %s [-v] [-s matrix_size|-i input_file|-p platform|-d device]\n", argv[0]);
      exit(EXIT_FAILURE);
  }

  if (input_file) {
      printf("Reading matrix from file %s\n", input_file);
      ret = create_matrix_from_file(&m, input_file, &matrix_dim);
      if (ret != RET_SUCCESS) {
          m = NULL;
          fprintf(stderr, "error create matrix from file %s\n", input_file);
          exit(EXIT_FAILURE);
      }
  } else {
    printf("No input file specified!\n");
    exit(EXIT_FAILURE);
  }

  if (do_verify){
    printf("Before LUD\n");
    print_matrix(m, matrix_dim);
    matrix_duplicate(m, &mm, matrix_dim);
  }

//  errcode = clGetPlatformIDs(NUM_PLATFORM, clPlatform, NULL);
//  CHECKERR(errcode);
//
//  errcode = clGetDeviceIDs(clPlatform[PLATFORM_ID], USEGPU ? CL_DEVICE_TYPE_GPU : CL_DEVICE_TYPE_CPU, 1, &clDevice, NULL);
//  CHECKERR(errcode);
clDevice = GetDevice(platform_id, device_id);
 size_t max_worksize[3];
 errcode = clGetDeviceInfo(clDevice, CL_DEVICE_MAX_WORK_ITEM_SIZES,sizeof(size_t)*3, &max_worksize, NULL);
 CHECKERR(errcode);
	while(BLOCK_SIZE*BLOCK_SIZE>max_worksize[0])
		BLOCK_SIZE = BLOCK_SIZE/2;
	
  clContext = clCreateContext(NULL, 1, &clDevice, NULL, NULL, &errcode);
  CHECKERR(errcode);

  clCommands = clCreateCommandQueue(clContext, clDevice, CL_QUEUE_PROFILING_ENABLE, &errcode);
  CHECKERR(errcode);

  kernelFile = fopen("lud_kernel.cl", "r");
  fseek(kernelFile, 0, SEEK_END);
  kernelLength = (size_t) ftell(kernelFile);
  kernelSource = (char *) malloc(sizeof(char)*kernelLength);
  rewind(kernelFile);
  fread((void *) kernelSource, kernelLength, 1, kernelFile);
  fclose(kernelFile);

  clProgram = clCreateProgramWithSource(clContext, 1, (const char **) &kernelSource, &kernelLength, &errcode);
  CHECKERR(errcode);

  free(kernelSource);
	char arg[100];
	sprintf(arg,"-D BLOCK_SIZE=%d", (int)BLOCK_SIZE);
  errcode = clBuildProgram(clProgram, 1, &clDevice, arg, NULL, NULL);
  if (errcode == CL_BUILD_PROGRAM_FAILURE)                                                                                                                                       
  {
    char *log;
    size_t logLength;
    errcode = clGetProgramBuildInfo(clProgram, clDevice, CL_PROGRAM_BUILD_LOG, 0, NULL, &logLength);
    log = (char *) malloc(sizeof(char)*logLength);
    errcode = clGetProgramBuildInfo(clProgram, clDevice, CL_PROGRAM_BUILD_LOG, logLength, (void *) log, NULL);
    fprintf(stderr, "Kernel build error! Log:\n%s", log);
    free(log);
    return;
  }
  CHECKERR(errcode);

  clKernel_diagonal = clCreateKernel(clProgram, "lud_diagonal", &errcode);
  CHECKERR(errcode);
  clKernel_perimeter = clCreateKernel(clProgram, "lud_perimeter", &errcode);
  CHECKERR(errcode);
  clKernel_internal = clCreateKernel(clProgram, "lud_internal", &errcode);
  CHECKERR(errcode);

  d_m = clCreateBuffer(clContext, CL_MEM_READ_WRITE, matrix_dim*matrix_dim*sizeof(float), NULL, &errcode);
  CHECKERR(errcode);

  /* beginning of timing point */
  stopwatch_start(&sw);
	 
  errcode = clEnqueueWriteBuffer(clCommands, d_m, CL_TRUE, 0, matrix_dim*matrix_dim*sizeof(float), (void *) m, 0, NULL, &ocdTempEvent);

  clFinish(clCommands);
  	START_TIMER(ocdTempEvent, OCD_TIMER_H2D, "Matrix Copy", ocdTempTimer)
	END_TIMER(ocdTempTimer)
	CHECKERR(errcode);

  int i=0;
  size_t localWorkSize[2];
  size_t globalWorkSize[2];
	//printf("BLOCK_SIZE: %d\n",BLOCK_SIZE);	
//	printf("max Work-item Size: %d\n",(int)max_worksize[0]);	
  	#ifdef START_POWER
	for( int iter = 0; iter < 1000; iter++)
	#endif
	for (i=0; i < matrix_dim-BLOCK_SIZE; i += BLOCK_SIZE) {
      errcode = clSetKernelArg(clKernel_diagonal, 0, sizeof(cl_mem), (void *) &d_m);
      errcode |= clSetKernelArg(clKernel_diagonal, 1, sizeof(int), (void *) &matrix_dim);
      errcode |= clSetKernelArg(clKernel_diagonal, 2, sizeof(int), (void *) &i);
      CHECKERR(errcode);
      
      localWorkSize[0] = BLOCK_SIZE;
      globalWorkSize[0] = BLOCK_SIZE;
      	 
	errcode = clEnqueueNDRangeKernel(clCommands, clKernel_diagonal, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
        clFinish(clCommands);
      	START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "Diagonal Kernels", ocdTempTimer)
	END_TIMER(ocdTempTimer)
	CHECKERR(errcode);
      errcode = clSetKernelArg(clKernel_perimeter, 0, sizeof(cl_mem), (void *) &d_m);
      errcode |= clSetKernelArg(clKernel_perimeter, 1, sizeof(int), (void *) &matrix_dim);
      errcode |= clSetKernelArg(clKernel_perimeter, 2, sizeof(int), (void *) &i);
      CHECKERR(errcode);
      localWorkSize[0] = BLOCK_SIZE*2;
      globalWorkSize[0] = ((matrix_dim-i)/BLOCK_SIZE-1)*localWorkSize[0];
      	 
	errcode = clEnqueueNDRangeKernel(clCommands, clKernel_perimeter, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
        clFinish(clCommands);
      START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "Perimeter Kernel", ocdTempTimer)
	CHECKERR(errcode);
	END_TIMER(ocdTempTimer)
      errcode = clSetKernelArg(clKernel_internal, 0, sizeof(cl_mem), (void *) &d_m);
      errcode |= clSetKernelArg(clKernel_internal, 1, sizeof(int), (void *) &matrix_dim);
      errcode |= clSetKernelArg(clKernel_internal, 2, sizeof(int), (void *) &i);
      CHECKERR(errcode);
      localWorkSize[0] = BLOCK_SIZE;
      localWorkSize[1] = BLOCK_SIZE;
      globalWorkSize[0] = ((matrix_dim-i)/BLOCK_SIZE-1)*localWorkSize[0];
      globalWorkSize[1] = ((matrix_dim-i)/BLOCK_SIZE-1)*localWorkSize[1];
      	 
	errcode = clEnqueueNDRangeKernel(clCommands, clKernel_internal, 2, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
        clFinish(clCommands);
      	START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "Internal Kernel", ocdTempTimer)
	END_TIMER(ocdTempTimer)
	CHECKERR(errcode);
  }
  errcode = clSetKernelArg(clKernel_diagonal, 0, sizeof(cl_mem), (void *) &d_m);
  errcode |= clSetKernelArg(clKernel_diagonal, 1, sizeof(int), (void *) &matrix_dim);
  errcode |= clSetKernelArg(clKernel_diagonal, 2, sizeof(int), (void *) &i);
  CHECKERR(errcode);
  localWorkSize[0] = BLOCK_SIZE;
  globalWorkSize[0] = BLOCK_SIZE;
  	 
	errcode = clEnqueueNDRangeKernel(clCommands, clKernel_diagonal, 1, NULL, globalWorkSize, localWorkSize, 0, NULL, &ocdTempEvent);
        clFinish(clCommands);
 	START_TIMER(ocdTempEvent, OCD_TIMER_KERNEL, "Diagonal Kernels", ocdTempTimer)
	 CHECKERR(errcode);

  END_TIMER(ocdTempTimer)
	 
  errcode = clEnqueueReadBuffer(clCommands, d_m, CL_TRUE, 0, matrix_dim*matrix_dim*sizeof(float), (void *) m, 0, NULL, &ocdTempEvent);
        clFinish(clCommands);
	START_TIMER(ocdTempEvent, OCD_TIMER_D2H, "Matrix copy", ocdTempTimer)
	END_TIMER(ocdTempTimer)
  /* end of timing point */
  stopwatch_stop(&sw);
  printf("Time consumed(ms): %lf\n", 1000*get_interval_by_sec(&sw));

  clReleaseMemObject(d_m);

  if (do_verify){
    printf("After LUD\n");
    print_matrix(m, matrix_dim);
    printf(">>>Verify<<<<\n");
	printf("matrix_dim: %d\n",matrix_dim);
    lud_verify(mm, m, matrix_dim); 
    free(mm);
  }

  clReleaseKernel(clKernel_diagonal);
  clReleaseKernel(clKernel_perimeter);
  clReleaseKernel(clKernel_internal);
  clReleaseProgram(clProgram);
  clReleaseCommandQueue(clCommands);
  clReleaseContext(clContext);

  free(m);
  OCD_FINISH
  return EXIT_SUCCESS;
}				/* ----------  end of function main  ---------- */
